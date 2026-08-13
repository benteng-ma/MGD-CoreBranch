import os
import pandas as pd
import gseapy as gp

# -----------------------------
# Input and output paths
# -----------------------------
score_file = r"D:\MGD_cNMF\results\GSE274498_allMG_pilot\GSE274498_allMG_pilot.gene_spectra_score.k_8.dt_0_1.txt"

outdir = r"D:\MGD_cNMF\results\GSE274498_allMG_pilot\K8_enrichment_top100"
os.makedirs(outdir, exist_ok=True)

# -----------------------------
# Read cNMF gene spectra scores
# -----------------------------
df = pd.read_csv(score_file, sep="\t", index_col=0)

print("Score matrix shape:", df.shape)

# Experimental gene universe
background = list(dict.fromkeys(
    str(g).upper() for g in df.columns
))

print("Background genes:", len(background))

# -----------------------------
# Download gene-set libraries
# -----------------------------
print("Downloading GO Biological Process 2026...")
go_bp = gp.get_library(
    name="GO_Biological_Process_2026",
    organism="Mouse"
)

print("GO terms:", len(go_bp))

print("Downloading Reactome Pathways 2024...")
reactome = gp.get_library(
    name="Reactome_Pathways_2024",
    organism="Mouse"
)

print("Reactome terms:", len(reactome))

libraries = {
    "GO_Biological_Process_2026": go_bp,
    "Reactome_Pathways_2024": reactome
}

summary_tables = []

# -----------------------------
# Run enrichment for P1-P8
# -----------------------------
for program in df.index:

    program_name = f"P{program}"

    genes = list(dict.fromkeys(
        str(g).upper()
        for g in df.loc[program].nlargest(100).index
    ))

    # Save Top100 genes
    gene_file = os.path.join(
        outdir,
        f"{program_name}_Top100_genes.txt"
    )

    pd.Series(genes).to_csv(
        gene_file,
        index=False,
        header=False
    )

    print()
    print("=" * 60)
    print(program_name)
    print("=" * 60)

    for library_name, library in libraries.items():

        enr = gp.enrich(
            gene_list=genes,
            gene_sets=library,
            background=background,
            outdir=None,
            no_plot=True,
            verbose=False
        )

        results = enr.results.copy()

        results.insert(0, "Program", program_name)
        results.insert(1, "Library", library_name)

        result_file = os.path.join(
            outdir,
            f"{program_name}_{library_name}.csv"
        )

        results.to_csv(
            result_file,
            index=False
        )

        n_sig = int(
            (results["Adjusted P-value"] < 0.05).sum()
        )

        print(
            library_name,
            "significant FDR<0.05:",
            n_sig
        )

        top10 = (
            results
            .sort_values("Adjusted P-value")
            .head(10)
            .copy()
        )

        summary_tables.append(top10)

# -----------------------------
# Combined Top10 summary
# -----------------------------
summary = pd.concat(
    summary_tables,
    ignore_index=True
)

summary_file = os.path.join(
    outdir,
    "K8_all_programs_Top10_enrichment.csv"
)

summary.to_csv(
    summary_file,
    index=False
)

print()
print("=" * 60)
print("DONE")
print("=" * 60)
print("Results saved to:")
print(outdir)
print()
print("Summary file:")
print(summary_file)