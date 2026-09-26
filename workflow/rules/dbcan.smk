dbcan_db_dir = Path("resources/dbcan_db")

rule install_dbcan_db:
    output:
        hmm = dbcan_db_dir / "dbCAN-HMMdb.txt",
        pressed = dbcan_db_dir / "dbCAN-HMMdb.txt.h3i"
    conda:
        "../envs/dbcan.yaml"
    log:
        "logs/dbcan/install_dbcan_db.log"
    shell:
        """
        mkdir -p {dbcan_db_dir}
        curl -sSL "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm" -o {output.hmm} 2>> {log}
        hmmpress -f {output.hmm} &>> {log}
        """


rule dbcan_annotate:
    input:
        faa = "data/interim/prokka/{strains}/{strains}.faa",
        hmm = rules.install_dbcan_db.output.hmm
    output:
        tsv = "data/interim/dbcan/{strains}/{strains}_cazyme.tsv"
    conda:
        "../envs/dbcan.yaml"
    log:
        "logs/dbcan/dbcan_annotate_{strains}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/run_dbcan_pyhmmer.py {input.faa} {input.hmm} {output.tsv} 2>> {log}
        """

rule dbcan_summary:
    input:
        tsv = lambda wildcards: [f"data/interim/dbcan/{s}/{s}_cazyme.tsv" for s in PEP_PROJECTS[wildcards.name].sample_table["genome_id"]]
    output:
        csv = "data/processed/{name}/tables/df_cazyme.csv"

    conda:
        "../envs/dbcan.yaml"
    log:
        "logs/dbcan/dbcan_summary_{name}.log"
    shell:
        """
        python -c "
import pandas as pd
from pathlib import Path
inputs = [Path(p) for p in '{input.tsv}'.split()]
dfs = [pd.read_csv(p, sep='\t') for p in inputs if p.exists() and p.stat().st_size > 0]
if dfs:
    df_all = pd.concat(dfs, ignore_index=True).drop_duplicates()
else:
    df_all = pd.DataFrame(columns=['gene_id', 'cazyme_family', 'cazyme_class', 'evalue', 'score', 'coverage', 'hmm_from', 'hmm_to'])
df_all.to_csv('{output.csv}', index=False)
" 2>> {log}
        """

rule enrich_ppanggolin_gexf:
    input:
        gexf = "data/processed/{name}/ppanggolin/genome_roary/gexf/pangenomeGraph.gexf",
        pangenome_csv = "data/processed/{name}/ppanggolin/genome_roary/gene_pres_abs",
        cazyme_csv = "data/processed/{name}/tables/df_cazyme.csv",
        bgc_dir = f"data/interim/bgcs/{{name}}/{dependency_version['antismash']}"
    output:
        gexf_annotated = "data/processed/{name}/ppanggolin/genome_roary/gexf/pangenomeGraph_annotated.gexf",
        overlap_csv = "data/processed/{name}/tables/df_bgc_cazyme_overlap.csv"
    conda:
        "../envs/dbcan.yaml"
    log:
        "logs/ppanggolin/enrich_gexf_{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/inject_gexf_metadata.py {input.gexf} {input.pangenome_csv} {input.cazyme_csv} {output.gexf_annotated} {input.bgc_dir} {output.overlap_csv} 2>> {log}
        """

