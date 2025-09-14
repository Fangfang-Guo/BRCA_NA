import pandas as pd
import os
from Bio.Seq import Seq
from collections import defaultdict
import requests

try:
    from pyfaidx import Fasta
except Exception:
    Fasta = None

from SigProfilerExtractor import sigpro as sig

# -----------------------------------------------------------
# Step 1. Initialize empty mutation matrix from COSMIC catalog
# -----------------------------------------------------------
file_path = "COSMIC_v3.4_SBS_GRCh37.txt"  # Reference COSMIC file
data = pd.read_csv(file_path, sep="\t")

# Use mutation type labels as row indices
mutation_types = data['Type'].tolist()
empty_matrix = pd.DataFrame(index=mutation_types)

# Save as initial empty matrix
empty_matrix.to_csv("mutation_matrix.csv")
print("Empty mutation matrix saved as mutation_matrix.csv")


# -----------------------------------------------------------
# Step 2. Utility functions for mutation context extraction
# -----------------------------------------------------------
UCSC_BASE_URL = "https://api.genome.ucsc.edu/getData/sequence"

def get_context_from_fasta(fasta, chrom, pos_1based, ref):
    """Extract trinucleotide context using local FASTA file."""
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    start = pos_1based - 1
    end = pos_1based + 1
    seq = fasta[chrom][start - 1:end].seq.upper()
    if len(seq) == 3 and seq[1] == ref:
        return seq
    return None

def get_context_from_ucsc(genome, chrom, pos_1based, ref):
    """Extract trinucleotide context using UCSC API."""
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    start0 = pos_1based - 2
    end0 = pos_1based + 1
    url = f"{UCSC_BASE_URL}?genome={genome};chrom={chrom};start={start0};end={end0}"
    r = requests.get(url, timeout=15)
    if r.status_code == 200:
        ctx = r.json().get("dna", "").upper()
        if len(ctx) == 3 and ctx[1] == ref:
            return ctx
    return None

def complementary_mutation_key(mutation_key):
    """Return complementary mutation representation."""
    left, ref, alt, right = mutation_key[0], mutation_key[2], mutation_key[4], mutation_key[6]
    return f"{str(Seq(left).complement())}[{str(Seq(ref).complement())}>{str(Seq(alt).complement())}]{str(Seq(right).complement())}"

def read_maf(maf_path):
    """Read MAF file into DataFrame."""
    compression = "gzip" if maf_path.endswith(".gz") else None
    return pd.read_csv(maf_path, sep="\t", comment="#", compression=compression, low_memory=False)

def detect_barcode_column(df):
    """Detect sample barcode column in MAF file."""
    for col in ["Tumor_Sample_Barcode", "Tumor_Sample", "Tumor_Sample_UUID"]:
        if col in df.columns:
            return col
    raise ValueError("Sample ID column not found in MAF.")


# -----------------------------------------------------------
# Step 3. Update mutation matrix with variants from MAF files
# -----------------------------------------------------------
def update_mutation_matrix_from_combined_maf(
    mutation_matrix_path,
    combined_maf_path,
    genome="hg38",
    fasta_path=None,
):
    mm = pd.read_csv(mutation_matrix_path, index_col=0)
    mm.index = mm.index.astype(str)

    maf = read_maf(combined_maf_path)
    barcode_col = detect_barcode_column(maf)

    # Keep only SNPs
    maf = maf[(maf["Reference_Allele"].astype(str).str.len() == 1) &
              (maf["Tumor_Seq_Allele2"].astype(str).str.len() == 1)]

    # Add missing sample columns
    all_samples = maf[barcode_col].dropna().astype(str).unique().tolist()
    missing = [s for s in all_samples if s not in mm.columns]
    if missing:
        mm = pd.concat([mm, pd.DataFrame(0, index=mm.index, columns=missing, dtype=int)], axis=1)

    # Load FASTA if provided
    fasta = None
    if fasta_path:
        if Fasta is None:
            raise RuntimeError("pyfaidx not installed. Please run `pip install pyfaidx`.")
        fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)

    ctx_cache = {}
    groups = list(maf.groupby(barcode_col))
    total = len(groups)

    for i, (sample, df_s) in enumerate(groups, start=1):
        if pd.isna(sample):
            continue
        sample = str(sample)
        counts = defaultdict(int)

        for _, row in df_s.iterrows():
            chrom = str(row["Chromosome"])
            pos = int(row["Start_Position"])
            ref = str(row["Reference_Allele"]).upper()
            alt = str(row["Tumor_Seq_Allele2"]).upper()
            if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                continue

            key = (chrom, pos, ref)
            if key in ctx_cache:
                ctx = ctx_cache[key]
            else:
                ctx = get_context_from_fasta(fasta, chrom, pos, ref) if fasta else \
                      get_context_from_ucsc(genome, chrom, pos, ref)
                ctx_cache[key] = ctx

            if not ctx or len(ctx) != 3:
                continue

            mut_key = f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"
            comp_key = complementary_mutation_key(mut_key)

            if mut_key in mm.index:
                counts[mut_key] += 1
            elif comp_key in mm.index:
                counts[comp_key] += 1

        if counts:
            s = pd.Series(counts, dtype=int)
            mm.loc[s.index, sample] = mm.loc[s.index, sample].to_numpy() + s.values

        print(f"[{i}/{total}] {sample}: {sum(counts.values())} mutations")

    mm = mm.copy()
    mm.to_csv(mutation_matrix_path)
    print(f"Mutation matrix updated: {mutation_matrix_path}")


# -----------------------------------------------------------
# Step 4. Build matrix and run SigProfilerExtractor
# -----------------------------------------------------------

# Update mutation matrix for two cohorts
update_mutation_matrix_from_combined_maf(
    mutation_matrix_path="mutation_matrix.csv",
    combined_maf_path="filtered_CA_final.maf",
    genome="hg38"
)

update_mutation_matrix_from_combined_maf(
    mutation_matrix_path="mutation_matrix.csv",
    combined_maf_path="filtered_NA_final.maf",
    genome="hg19"
)

# Load final matrix
df = pd.read_csv("mutation_matrix.csv", index_col=0)

# Identify zero-mutation samples
zero_cols = df.columns[df.sum(axis=0) == 0]
print("Samples with zero mutations:", list(zero_cols))
print(f"Total {len(zero_cols)} samples with zero mutations")

# Save as tab-delimited
df.to_csv("mutation_matrix.txt", sep="\t")

# Run SigProfilerExtractor
sig.sigProfilerExtractor(
    input_type="matrix",
    output="mutational_signature_results",  # Output folder
    input_data="mutation_matrix.txt",       # Input matrix
    minimum_signatures=1,
    maximum_signatures=10,
    nmf_replicates=100,
    cpu=-1
)
