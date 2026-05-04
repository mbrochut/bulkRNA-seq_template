#!/usr/bin/env python3

import os
import json
import argparse
import gseapy as gp


def create_directories(organism):
    dirs = [
        f"data/DB/{organism}",
        "meta",
        "QUARTO",
        "figures",
        "figures/QC",
        "results/contrasts",
        "results/models",
        "results/pathway/concat",
        "results/pathway/GSEA",
        "results/pathway/ORA",
        "results/pathway/GSEA_object",
        "results/QC",
        "results/papermill"
    ]

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    print("✅ Directory structure ready.")


def is_db_empty(path):
    return not any(os.scandir(path))


def download_databases(organism):
    outdir = f"data/DB/{organism}"

    libs = [
        "GO_Biological_Process_2025",
        "MSigDB_Hallmark_2020",
    ]

    if organism.lower() == "mouse":
        libs.append("KEGG_2019_Mouse")
    else:
        libs.append("KEGG_2026")

    downloaded = {}

    for lib in libs:
        try:
            print(f"Downloading {lib}...")
            gene_sets = gp.get_library(name=lib, organism=organism)

            outfile = os.path.join(outdir, f"{lib}.json")
            with open(outfile, "w") as f:
                json.dump(gene_sets, f)

            downloaded[lib] = "OK"

        except Exception as e:
            print(f"❌ Failed: {lib} → {e}")
            downloaded[lib] = str(e)

    # Save log
    with open(os.path.join(outdir, "download_log.json"), "w") as f:
        json.dump(downloaded, f, indent=2)

    print("✅ Download complete.")


def main():
    parser = argparse.ArgumentParser(description="Initialize project structure")

    parser.add_argument(
        "--organism",
        nargs="?",
        default="Human",
        help="Organism (e.g. Human, Mouse)",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if DB exists",
    )

    args = parser.parse_args()
    organism = args.organism

    print(f"Initializing repo (organism: {organism})")

    create_directories(organism)

    db_path = f"data/DB/{organism}"

    if is_db_empty(db_path) or args.force:
        print(f"DB for {organism} is empty or force=True → downloading...")
        download_databases(organism)
    else:
        print(f"DB for {organism} already exists → skipping download.")

    print("✅ Init complete!")


if __name__ == "__main__":
    main()