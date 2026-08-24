antismash_db_path = Path(f"resources/antismash_db")

if antismash_major_version <= 6:
    rule antismash_db_setup:
        output:
            clusterblast=directory(antismash_db_path / "clusterblast"),
            clustercompare_mibig=directory(antismash_db_path / "clustercompare/mibig"),
            pfam=directory(antismash_db_path / "pfam/34.0"),
            resfam=directory(antismash_db_path / "resfam"),
            tigrfam=directory(antismash_db_path / "tigrfam"),
        conda:
            "../envs/antismash_v6.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        params:
            antismash_db_path=antismash_db_path
        shell:
            """
            download-antismash-databases --database-dir {params.antismash_db_path} &>> {log}
            antismash --version >> {log}
            antismash --check-prereqs --databases {params.antismash_db_path} &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            clusterblast=rules.antismash_db_setup.output.clusterblast,
            clustercompare_mibig=rules.antismash_db_setup.output.clustercompare_mibig,
            pfam=rules.antismash_db_setup.output.pfam,
            resfam=rules.antismash_db_setup.output.resfam,
            tigrfam=rules.antismash_db_setup.output.tigrfam
        output:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
            json="data/interim/antismash/{version}/{strains}/{strains}.json",
            zip="data/interim/antismash/{version}/{strains}/{strains}.zip",
        conda:
            "../envs/antismash_v6.yaml"
        threads: 4
        log:
            "logs/antismash/antismash/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            genefinding="none",
            antismash_db_path=antismash_db_path,
        shell:
            """
            antismash \
                --genefinding-tool {params.genefinding} \
                --database {params.antismash_db_path} \
                --output-dir {params.folder} \
                --cb-general \
                --cb-subclusters \
                --cb-knownclusters \
                -c {threads} {input.gbk} --logfile {log} 2>> {log}
            """

elif antismash_major_version >= 7:
    rule antismash_db_setup:
        output:
            as_js=directory(antismash_db_path / "as-js/0.16"),
            clusterblast=directory(antismash_db_path / "clusterblast"),
            clustercompare_mibig=directory(antismash_db_path / "clustercompare/mibig/4.0"),
            comparippson_asdb=directory(antismash_db_path / "comparippson/asdb/4.0"),
            comparippson_mibig=directory(antismash_db_path / "comparippson/mibig/4.0"),
            knownclusterblast=directory(antismash_db_path / "knownclusterblast/4.0"),
            nrps_pks_stachelhaus=directory(antismash_db_path / "nrps_pks/stachelhaus/1.1"),
            nrps_pks_svm=directory(antismash_db_path / "nrps_pks/svm/2.0"),
            nrps_pks_transATor=directory(antismash_db_path / "nrps_pks/transATor/2023.02.23"),
            pfam=directory(antismash_db_path / "pfam/35.0"),
            resfam=directory(antismash_db_path / "resfam"),
            tigrfam=directory(antismash_db_path / "tigrfam"),
            mite=directory(antismash_db_path / "mite/1.3"),
        conda:
            "../envs/antismash.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        params:
            antismash_db_path=antismash_db_path
        shell:
            """
            download-antismash-databases --database-dir {params.antismash_db_path} &>> {log}
            antismash --version >> {log}
            antismash --database {params.antismash_db_path} --prepare-data &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            as_js=rules.antismash_db_setup.output.as_js,
            clusterblast=rules.antismash_db_setup.output.clusterblast,
            clustercompare_mibig=rules.antismash_db_setup.output.clustercompare_mibig,
            comparippson_asdb=rules.antismash_db_setup.output.comparippson_asdb,
            comparippson_mibig=rules.antismash_db_setup.output.comparippson_mibig,
            knownclusterblast=rules.antismash_db_setup.output.knownclusterblast,
            nrps_pks_stachelhaus=rules.antismash_db_setup.output.nrps_pks_stachelhaus,
            nrps_pks_svm=rules.antismash_db_setup.output.nrps_pks_svm,
            nrps_pks_transATor=rules.antismash_db_setup.output.nrps_pks_transATor,
            pfam=rules.antismash_db_setup.output.pfam,
            resfam=rules.antismash_db_setup.output.resfam,
            tigrfam=rules.antismash_db_setup.output.tigrfam
        output:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
            json="data/interim/antismash/{version}/{strains}/{strains}.json",
            zip="data/interim/antismash/{version}/{strains}/{strains}.zip",
        conda:
            "../envs/antismash.yaml"
        threads: 4
        log:
            "logs/antismash/antismash/{version}/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            antismash_db_path=antismash_db_path,
            genefinding="none",
            taxon=os.getenv("BGCFLOW_ANTISMASH_MODE", "bacteria")
        shell:
            """
            set +e

            # Find the latest existing JSON output for this strain
            latest_json=$(find data/interim/antismash/*/* -name "{wildcards.strains}.json" | sort -V | tail -n 1) 2>> {log}
            if [ -n "$latest_json" ]; then
                # Use existing JSON result as starting point
                echo "Using existing JSON from $latest_json as starting point..." >> {log}
                antismash_input="--reuse-result $latest_json"
            else
                # No existing JSON result found, use genbank input
                echo "No existing JSON result found, starting AntiSMASH from scratch..." >> {log}
                antismash_input="{input.gbk}"
            fi

            # Store common parameters in a variable
            antismash_params="--genefinding-tool {params.genefinding} \
                --output-dir {params.folder} \
                --database {params.antismash_db_path} \
                --taxon {params.taxon} \
                --cb-knownclusters \
                --cb-subclusters \
                --cc-mibig \
                --clusterhmmer \
                --tigrfam \
                --pfam2go \
                --rre \
                --asf \
                --tfbs \
                -c {threads} \
                --logfile {log}"

            # Run AntiSMASH
            echo "Running AntiSMASH {params.taxon} mode..." >> {log}
            antismash $antismash_params $antismash_input 2>> {log}

            # Check if the run failed due to changed detection results or changed protocluster types
            if grep -q -e "ValueError: Detection results have changed. No results can be reused" \
                    -e "RuntimeError: Protocluster types supported by HMM detection have changed, all results invalid" {log}
            then
                # Use genbank input instead
                echo "Previous JSON result is invalid, starting AntiSMASH from scratch..." >> {log}
                antismash $antismash_params {input.gbk} 2>> {log}
            fi
            """

rule copy_antismash:
    input:
        dir="data/interim/antismash/{version}/{strains}",
    output:
        dir=directory("data/processed/{name}/antismash/{version}/{strains}"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/antismash/copy_antismash/copy_antismash_{version}-{strains}-{name}.log",
    shell:
        """
        base_dir=$PWD
        mkdir {output.dir}
        (cd {output.dir} && for item in $(ls $base_dir/{input.dir}); do ln -s $base_dir/{input.dir}/$item $(basename $item); done) 2>> {log}
        """

rule antismash_json_extract:
    input:
        json = "data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        cdss = temp("data/interim/database/as_{version}/{strains}/{strains}_cdss.json"),
        dna_sequences = temp("data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json"),
        regions = temp("data/interim/database/as_{version}/{strains}/{strains}_regions.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/scatter/as_{version}_json_extract_{strains}.log"
    params:
        outdir = "data/interim/database/as_{version}/{strains}",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/bgc_meta.py {input.json} {params.outdir} {wildcards.strains} 2>> {log}
        """

rule build_dna_sequences_table:
    input:
        dna_sequences = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/gather/as_{version}_dna_sequences_gather_{name}.log"
    params:
        index_key = "sequence_id",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.dna_sequences}' {params.index_key} {output.dna_sequences} 2>> {log}
        """

rule build_regions_table:
    input:
        regions = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_regions.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/gather/as_{version}_regions_gather_{name}.log"
    params:
        index_key = "region_id",
        exclude = "_regions.json"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet_with_correction.py '{input.regions}' {input.mapping_dir} {output.regions} {params.exclude} {params.index_key} 2>> {log}
 2>> {log}
        """

rule build_cdss_table:
    input:
        cdss = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/gather/as_{version}_cdss_gather_{name}.log"
    params:
        index_key = "cds_id",
        exclude = "_cdss.json"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet_with_correction.py '{input.cdss}' {input.mapping_dir} {output.cdss} {params.exclude} {params.index_key} 2>> {log}
 2>> {log}
        """

rule build_warehouse:
    input:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        log = "data/processed/{name}/data_warehouse/{version}/database.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/report/database_{version}_{name}.log"
    shell:
        """
        echo {input.cdss} >> {output.log}
        echo {input.regions} >> {output.log}
        echo {input.dna_sequences} >> {output.log}
        """
