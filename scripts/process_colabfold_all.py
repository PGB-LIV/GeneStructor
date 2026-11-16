import argparse
import concurrent.futures as cf
import json
import math
import os
import re
from pathlib import Path

import pandas as pd
from tqdm import tqdm


ALLOWED_EXTENSIONS = {".json", ".pdb", ".cif"}


def find_key_recursive(obj, target_key):
    if isinstance(obj, dict):
        if target_key in obj:
            return obj[target_key]
        for value in obj.values():
            found = find_key_recursive(value, target_key)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for value in obj:
            found = find_key_recursive(value, target_key)
            if found is not None:
                return found
    return None


def to_float_or_none(value):
    try:
        if value is None:
            return None
        num = float(value)
        if math.isnan(num):
            return None
        return num
    except (TypeError, ValueError):
        return None


def summarize_plddt(plddt_values):
    numeric = [to_float_or_none(x) for x in plddt_values]
    numeric = [x for x in numeric if x is not None]
    if not numeric:
        return {
            "mean_plddt": None,
            "plddt_50": None,
            "plddt_60": None,
            "plddt_70": None,
            "plddt_80": None,
            "plddt_90": None,
            "total_residues": None,
        }

    return {
        "mean_plddt": sum(numeric) / len(numeric),
        "plddt_50": sum(x > 50 for x in numeric),
        "plddt_60": sum(x > 60 for x in numeric),
        "plddt_70": sum(x > 70 for x in numeric),
        "plddt_80": sum(x > 80 for x in numeric),
        "plddt_90": sum(x > 90 for x in numeric),
        "total_residues": len(numeric),
    }


def detect_source_from_text(path_obj, text_head):
    haystack = f"{path_obj.as_posix()}\n{text_head}".lower()
    if "colabfold" in haystack:
        return "colabfold"
    if "alphafold" in haystack or "af2" in haystack:
        return "alphafold2"
    return "other"


def extract_ptm_from_text(text):
    patterns = [
        r"\"ptm\"\s*[:=]\s*([0-9]*\.?[0-9]+)",
        r"\bptm\b\s*[:=]\s*([0-9]*\.?[0-9]+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return to_float_or_none(match.group(1))
    return None


def parse_json_file(path_obj):
    with path_obj.open("r", encoding="utf-8") as handle:
        obj = json.load(handle)

    ptm_val = find_key_recursive(obj, "ptm")
    plddt_val = find_key_recursive(obj, "plddt")
    if plddt_val is None:
        plddt_val = []

    if not isinstance(plddt_val, list):
        plddt_val = [plddt_val]

    metrics = summarize_plddt(plddt_val)
    text_head = json.dumps(obj)[:5000]

    return {
        "idx": path_obj.stem,
        "file_name": path_obj.name,
        "file_path": str(path_obj),
        "file_type": "json",
        "source": detect_source_from_text(path_obj, text_head),
        "ptm": to_float_or_none(ptm_val),
        **metrics,
        "error": None,
    }


def parse_pdb_file(path_obj):
    residue_scores = {}
    ptm = None

    with path_obj.open("r", encoding="utf-8", errors="ignore") as handle:
        head_lines = []
        for i, line in enumerate(handle):
            if i < 200:
                head_lines.append(line)

            if ptm is None:
                ptm = extract_ptm_from_text(line)

            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if len(line) < 66:
                continue

            chain_id = line[21:22].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26:27].strip() or "_"
            key = (chain_id, resseq, icode)

            try:
                b_factor = float(line[60:66].strip())
            except ValueError:
                continue

            residue_scores.setdefault(key, []).append(b_factor)

    per_residue = [sum(vals) / len(vals) for vals in residue_scores.values() if vals]
    metrics = summarize_plddt(per_residue)
    text_head = "".join(head_lines)

    return {
        "idx": path_obj.stem,
        "file_name": path_obj.name,
        "file_path": str(path_obj),
        "file_type": "pdb",
        "source": detect_source_from_text(path_obj, text_head),
        "ptm": ptm,
        **metrics,
        "error": None,
    }


def parse_cif_file(path_obj):
    residue_scores = {}
    ptm = None
    atom_count = 0
    in_atom_site = False
    col_indices = {}

    with path_obj.open("r", encoding="utf-8", errors="ignore") as handle:
        text_head_lines = []
        for i, line in enumerate(handle):
            if i < 250:
                text_head_lines.append(line)

            if ptm is None:
                ptm = extract_ptm_from_text(line)

            # Detect mmCIF _atom_site loop
            if line.strip().startswith("_atom_site."):
                in_atom_site = True
                col_name = line.strip().split()[0]
                if col_name == "_atom_site.B_iso_or_equiv":
                    col_indices["b_factor"] = len(col_indices)
                elif col_name == "_atom_site.label_asym_id":
                    col_indices["chain"] = len(col_indices)
                elif col_name == "_atom_site.label_seq_id":
                    col_indices["resseq"] = len(col_indices)
                continue

            # Try PDB-style ATOM/HETATM lines (hybrid format)
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) < 2:
                    continue

                b_factor = None
                # If we detected column indices, use them
                if col_indices and len(parts) > max(col_indices.values(), default=-1):
                    b_idx = col_indices.get("b_factor")
                    if b_idx is not None and b_idx < len(parts):
                        b_factor = to_float_or_none(parts[b_idx])
                
                # Fallback: find last numeric value (likely B-factor)
                if b_factor is None:
                    for token in reversed(parts):
                        value = to_float_or_none(token)
                        if value is not None:
                            b_factor = value
                            break
                
                if b_factor is None:
                    continue

                # Extract chain and residue identifiers
                if col_indices:
                    chain_idx = col_indices.get("chain", 6)
                    resseq_idx = col_indices.get("resseq", 8)
                    chain_id = parts[chain_idx] if chain_idx < len(parts) else "_"
                    resseq = parts[resseq_idx] if resseq_idx < len(parts) else "_"
                else:
                    # Fallback positions for hybrid format
                    chain_id = parts[6] if len(parts) > 6 else "_"
                    resseq = parts[8] if len(parts) > 8 else "_"
                
                key = (chain_id, resseq)
                residue_scores.setdefault(key, []).append(b_factor)
                atom_count += 1

    per_residue = [sum(vals) / len(vals) for vals in residue_scores.values() if vals]
    metrics = summarize_plddt(per_residue)
    text_head = "".join(text_head_lines)

    error_msg = None
    if atom_count == 0:
        error_msg = "No ATOM/HETATM records found in CIF file"
    elif len(per_residue) == 0:
        error_msg = "No valid B-factors extracted from CIF file"

    return {
        "idx": path_obj.stem,
        "file_name": path_obj.name,
        "file_path": str(path_obj),
        "file_type": "cif",
        "source": detect_source_from_text(path_obj, text_head),
        "ptm": ptm,
        **metrics,
        "error": error_msg,
    }


def parse_file(path_str):
    path_obj = Path(path_str)
    suffix = path_obj.suffix.lower()

    try:
        if suffix == ".json":
            return parse_json_file(path_obj)
        if suffix == ".pdb":
            return parse_pdb_file(path_obj)
        if suffix == ".cif":
            return parse_cif_file(path_obj)

        return {
            "idx": path_obj.stem,
            "file_name": path_obj.name,
            "file_path": str(path_obj),
            "file_type": suffix.replace(".", "") or "unknown",
            "source": "other",
            "ptm": None,
            "mean_plddt": None,
            "plddt_50": None,
            "plddt_60": None,
            "plddt_70": None,
            "plddt_80": None,
            "plddt_90": None,
            "total_residues": None,
            "error": f"Unsupported extension: {suffix}",
        }
    except Exception as exc:
        return {
            "idx": path_obj.stem,
            "file_name": path_obj.name,
            "file_path": str(path_obj),
            "file_type": suffix.replace(".", "") or "unknown",
            "source": "other",
            "ptm": None,
            "mean_plddt": None,
            "plddt_50": None,
            "plddt_60": None,
            "plddt_70": None,
            "plddt_80": None,
            "plddt_90": None,
            "total_residues": None,
            "error": str(exc),
        }


def chunked(items, chunk_size):
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Process JSON/CIF/PDB files and extract PTM and pLDDT metrics "
            "with laptop-friendly batching."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Directory containing .json/.pdb/.cif files (recursive search).",
    )
    parser.add_argument(
        "--output-csv",
        default="colabfold_structure_metrics.csv",
        help="Path to output CSV file.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Maximum CPU workers (default: 4 for laptop-friendly use).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="Files processed per chunk (default: 20).",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Show sample parsed data for manual verification.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit processing to first N files (for testing).",
    )
    return parser


def main():
    args = build_parser().parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    all_files = [
        p
        for p in input_dir.rglob("*")
        if p.is_file() and p.suffix.lower() in ALLOWED_EXTENSIONS
    ]

    if not all_files:
        print("No .json/.pdb/.cif files found. Nothing to process.")
        pd.DataFrame().to_csv(args.output_csv, index=False)
        return

    # Apply limit if specified (for testing)
    if args.limit and args.limit > 0:
        all_files = all_files[:args.limit]
        print(f"Limiting to first {len(all_files)} files (--limit {args.limit})")

    max_workers = max(1, min(int(args.max_workers), 4, os.cpu_count() or 1))
    batch_size = max(1, min(int(args.batch_size), 20))

    print(f"Found {len(all_files)} files to process")
    print(f"Using max_workers={max_workers}, batch_size={batch_size}")

    results = []
    batches = list(chunked(all_files, batch_size))

    for batch in tqdm(batches, desc="Processing batches"):
        with cf.ProcessPoolExecutor(max_workers=max_workers) as executor:
            batch_results = list(executor.map(parse_file, [str(x) for x in batch]))
        results.extend(batch_results)

    df = pd.DataFrame(results)

    preferred_columns = [
        "idx",
        "file_name",
        "file_path",
        "file_type",
        "source",
        "ptm",
        "mean_plddt",
        "total_residues",
        "plddt_50",
        "plddt_60",
        "plddt_70",
        "plddt_80",
        "plddt_90",
        "error",
    ]

    for col in preferred_columns:
        if col not in df.columns:
            df[col] = None
    df = df[preferred_columns]

    df.to_csv(args.output_csv, index=False)

    print(f"\nSaved results to: {args.output_csv}")
    print("\n=== Summary ===")
    print("\nFile types:")
    print(df["file_type"].value_counts(dropna=False).to_string())
    print("\nSources detected:")
    print(df["source"].value_counts(dropna=False).to_string())
    print(f"\nRows with errors: {df['error'].notna().sum()}")
    
    if df['error'].notna().sum() > 0:
        print("\nError breakdown:")
        print(df[df['error'].notna()]['error'].value_counts().to_string())
    
    # Show pLDDT statistics
    if df['mean_plddt'].notna().sum() > 0:
        print("\n=== pLDDT Statistics ===")
        print(f"Mean pLDDT: {df['mean_plddt'].mean():.2f}")
        print(f"Median pLDDT: {df['mean_plddt'].median():.2f}")
        print(f"Min pLDDT: {df['mean_plddt'].min():.2f}")
        print(f"Max pLDDT: {df['mean_plddt'].max():.2f}")
    
    # Validation mode: show sample data
    if args.validate:
        print("\n=== Validation Samples ===")
        for file_type in df['file_type'].unique():
            sample = df[df['file_type'] == file_type].head(2)
            print(f"\nSample {file_type} files:")
            for _, row in sample.iterrows():
                print(f"  {row['file_name']}:")
                print(f"    source: {row['source']}, ptm: {row['ptm']}, mean_plddt: {row['mean_plddt']}")
                print(f"    plddt_80: {row['plddt_80']}, total_residues: {row['total_residues']}")
                if row['error']:
                    print(f"    ERROR: {row['error']}")


if __name__ == "__main__":
    main()
