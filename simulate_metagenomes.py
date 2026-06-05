#!/usr/bin/env python3

import random
import zipfile
from pathlib import Path

import numpy as np

###############################################################################
# CONFIG
###############################################################################

BASE_DIR = Path("/home/andrew/Documents/Telost_simulated_microbiomes")

FISHES = [
    "Oncorhynchus_mykiss",
    "Perca_fluviatilis",
    "Salmo_salar"
]

N_CONTROL = 20
N_INFECTED = 20

HOST_RANGE = (0.005, 0.02)
PATHOGEN_RANGE = (0.13, 0.22, 0.37)

OUT_DIR = BASE_DIR / "simulated_output"

random.seed(42)
np.random.seed(42)

###############################################################################
# GENOME DISCOVERY (ROBUST FOR NCBI STRUCTURE)
###############################################################################

def discover_genomes(fish_dir):

    genomes = []

    for zip_path in Path(fish_dir).rglob("*.zip"):

        try:
            with zipfile.ZipFile(zip_path) as z:

                for f in z.namelist():

                    if not f.endswith((".fna", ".fasta")):
                        continue

                    try:
                        with z.open(f) as fh:

                            # scan until FASTA header
                            for _ in range(100):
                                line = fh.readline()
                                if not line:
                                    break

                                line = line.decode(errors="ignore").strip()

                                if line.startswith(">"):
                                    accession = line[1:].split()[0]

                                    genomes.append({
                                        "zip": zip_path,
                                        "inner": f,
                                        "accession": accession
                                    })
                                    break

                    except Exception:
                        continue

        except Exception:
            continue

    print(f"\n[DISCOVERY] {fish_dir.name}: {len(genomes)} genomes found")

    return genomes


###############################################################################
# CORE SELECTION (SAFE)
###############################################################################

def select_core(genomes, fish, n_core=10):

    if len(genomes) == 0:
        print(f"[WARNING] {fish}: no genomes found")
        return []

    random.shuffle(genomes)

    core = genomes[:min(n_core, len(genomes))]

    print(f"[CORE] {fish}: using {len(core)} genomes")

    return core


###############################################################################
# DIRICHLET
###############################################################################

def dirichlet(n, alpha=5):
    return np.random.dirichlet(np.repeat(alpha, n))


###############################################################################
# SAMPLE GENERATION
###############################################################################

def generate_sample(core, bg, fish, infected=False):

    abund = {}

    host = random.uniform(*HOST_RANGE)
    abund["HOST"] = host

    remaining = 1 - host

    pathogen = 0.0

    if infected:
        pathogen = np.random.triangular(*PATHOGEN_RANGE)
        abund["Aeromonas_salmonicida"] = pathogen
        remaining -= pathogen
        core_frac = remaining * 0.55
    else:
        core_frac = remaining * 0.80

    bg_frac = remaining - core_frac

    # CORE
    if len(core) > 0:

        w = dirichlet(len(core), alpha=10)

        for i, g in enumerate(core):
            abund[g["accession"]] = core_frac * w[i]

    # BACKGROUND
    if len(bg) > 0:

        wbg = dirichlet(len(bg), alpha=3)

        for i, g in enumerate(bg):
            abund[g["accession"]] = bg_frac * wbg[i]

    total = sum(abund.values())

    if total == 0:
        return {}

    return {k: v / total for k, v in abund.items()}


###############################################################################
# REFERENCE FASTA BUILDER (FIXED PER FISH)
###############################################################################

def build_reference(fish, genomes):

    ref_dir = OUT_DIR / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    ref_path = ref_dir / f"{fish}_combined.fasta"

    with open(ref_path, "w") as out:

        for g in genomes:

            try:
                with zipfile.ZipFile(g["zip"]) as z:
                    with z.open(g["inner"]) as fh:

                        seq = []
                        header = None

                        for line in fh:
                            line = line.decode(errors="ignore").strip()

                            if line.startswith(">"):

                                # flush previous genome
                                if header and seq:
                                    out.write(header + "\n")
                                    out.write("".join(seq) + "\n")

                                header = f">{g['accession']}"
                                seq = []

                            else:
                                seq.append(line)

                        # flush last record
                        if header and seq:
                            out.write(header + "\n")
                            out.write("".join(seq) + "\n")

            except Exception:
                continue

    print(f"[REFERENCE FIXED] {fish}: single-record-per-genome FASTA")

    return ref_path


###############################################################################
# RUN SCRIPT BUILDER (FIXED PER FISH)
###############################################################################

def write_run_script(fish, ref_path, sample_dir):

    run_dir = OUT_DIR / "runs"
    run_dir.mkdir(parents=True, exist_ok=True)

    reads_dir = OUT_DIR / "reads" / fish
    reads_dir.mkdir(parents=True, exist_ok=True)

    script_path = run_dir / f"run_iss_{fish}.sh"

    with open(script_path, "w") as f:

        f.write("#!/bin/bash\n\n")
        f.write(f"REF={ref_path}\n")
        f.write(f"OUT={reads_dir}\n\n")
        f.write("mkdir -p $OUT\n\n")

        for file in sample_dir.glob("*.txt"):

            name = file.stem

            f.write(
                "iss generate "
                f"--genomes $REF "
                f"--abundance_file {file} "
                f"--model hiseq "
                f"--output $OUT/{name}\n"
            )

    print(f"[SCRIPT] {fish}: {script_path}")
    
###############################################################################
# MAIN PIPELINE
###############################################################################

def main():

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for fish in FISHES:

        print(f"\n================ {fish} ================")

        fish_dir = BASE_DIR / "downloaded_genomes" / fish

        genomes = discover_genomes(fish_dir)

        print(f"[DEBUG] {fish}: genomes parsed = {len(genomes)}")

        core = select_core(genomes, fish)

        if len(core) == 0:
            print(f"[SKIP] {fish}: no usable genomes")
            continue

        core_ids = {c["accession"] for c in core}

        bg = [g for g in genomes if g["accession"] not in core_ids]

        # OUTPUT DIRS
        sample_dir = OUT_DIR / "abundance_files" / fish
        sample_dir.mkdir(parents=True, exist_ok=True)

        # STEP 1: reference FASTA (FIXED)
        ref_path = build_reference(fish, genomes)

        metadata = []

        # STEP 2: abundance generation
        for i in range(N_CONTROL + N_INFECTED):

            infected = i >= N_CONTROL

            n_bg = random.randint(10, 40) if infected else random.randint(30, 60)

            bg_sample = random.sample(bg, min(n_bg, len(bg))) if len(bg) > 0 else []

            abund = generate_sample(core, bg_sample, fish, infected)

            label = "infected" if infected else "control"

            out_file = sample_dir / f"{label}_{i+1:03d}.txt"

            with open(out_file, "w") as f:
                for k, v in abund.items():
                    f.write(f"{k}\t{v}\n")

            metadata.append((out_file.name, fish, label))

        # STEP 3: metadata
        with open(sample_dir / "metadata.tsv", "w") as f:
            f.write("sample\tfish\tstatus\n")
            for row in metadata:
                f.write("\t".join(row) + "\n")

        # STEP 4: run script (FIXED PER FISH)
        write_run_script(fish, ref_path, sample_dir)

    print("\nDONE.")


if __name__ == "__main__":
    main()
