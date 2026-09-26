import sys
from pathlib import Path
import pandas as pd
import pyhmmer

def run_dbcan_hmm(faa_path, hmm_db_path, output_tsv, evalue_cutoff=1e-15, coverage_cutoff=0.35):
    """
    Scan protein sequences in faa_path against dbCAN HMM database using pyhmmer.
    Filters by E-value and domain coverage standard for dbCAN.
    """
    faa_path = Path(faa_path)
    hmm_db_path = Path(hmm_db_path)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    results = []

    with pyhmmer.easel.SequenceFile(str(faa_path), digital=True) as seq_file:
        sequences = list(seq_file)

    with pyhmmer.plan7.HMMFile(str(hmm_db_path)) as hmm_file:
        hmms = list(hmm_file)

    hmm_lengths = {}
    for h in hmms:
        h_name = h.name.decode('utf-8') if isinstance(h.name, bytes) else str(h.name)
        hmm_lengths[h_name] = h.M

    # Run hmmsearch
    for hmm, top_hits in zip(hmms, pyhmmer.hmmsearch(hmms, sequences)):
        h_name = hmm.name.decode('utf-8') if isinstance(hmm.name, bytes) else str(hmm.name)
        cazyme_family = h_name.replace('.hmm', '')
        hmm_len = hmm.M

        for hit in top_hits:
            if hit.evalue <= evalue_cutoff:
                gene_id = hit.name.decode('utf-8') if isinstance(hit.name, bytes) else str(hit.name)
                for domain in hit.domains:
                    if domain.i_evalue <= evalue_cutoff:
                        # calculate domain coverage
                        domain_cov = (domain.alignment.hmm_to - domain.alignment.hmm_from + 1) / hmm_len
                        if domain_cov >= coverage_cutoff:
                            
                            # Determine CAZyme class
                            cazyme_class = "Unknown"
                            if cazyme_family.startswith("GH"):
                                cazyme_class = "Glycoside Hydrolase (GH)"
                            elif cazyme_family.startswith("GT"):
                                cazyme_class = "GlycosylTransferase (GT)"
                            elif cazyme_family.startswith("PL"):
                                cazyme_class = "Polysaccharide Lyase (PL)"
                            elif cazyme_family.startswith("CE"):
                                cazyme_class = "Carbohydrate Esterase (CE)"
                            elif cazyme_family.startswith("AA"):
                                cazyme_class = "Auxiliary Activity (AA)"
                            elif cazyme_family.startswith("CBM"):
                                cazyme_class = "Carbohydrate-Binding Module (CBM)"

                            results.append({
                                "gene_id": gene_id,
                                "cazyme_family": cazyme_family,
                                "cazyme_class": cazyme_class,
                                "evalue": domain.i_evalue,
                                "score": domain.score,
                                "coverage": domain_cov,
                                "hmm_from": domain.alignment.hmm_from,
                                "hmm_to": domain.alignment.hmm_to,
                            })

    df = pd.DataFrame(results)
    if not df.empty:
        # Keep best hit per gene
        df = df.sort_values(by=["gene_id", "evalue"]).drop_duplicates(subset=["gene_id", "cazyme_family"], keep="first")
    else:
        df = pd.DataFrame(columns=["gene_id", "cazyme_family", "cazyme_class", "evalue", "score", "coverage", "hmm_from", "hmm_to"])

    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Annotated {len(df)} CAZyme domain hits in {faa_path.name}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python run_dbcan_pyhmmer.py <faa_path> <hmm_db_path> <output_tsv>")
        sys.exit(1)
    run_dbcan_hmm(sys.argv[1], sys.argv[2], sys.argv[3])
