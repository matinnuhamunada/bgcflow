



rule get_dbt_template:
    input:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        profile = "data/processed/{name}/dbt/antiSMASH_{version}/models/sources.yml"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "logs/database/report/get_dbt_template_{version}_{name}.log"
    threads: 4
    params:
        dbt = "data/processed/{name}/dbt/antiSMASH_{version}",
        dbt_repo = "https://github.com/NBChub/bgcflow_dbt-duckdb",
        release = "0.4.0",
        cutoff = "0.30",
        as_version = "{version}"
    shell:
        """
        # clone dbt
        if [ -f "{params.dbt}/profiles.yml" ]
        then
            echo "{params.dbt} already exists!" >> {log}
        else
            rm -rf {params.dbt} 2>> {log}
            mkdir -p data/processed/{wildcards.name}/dbt
            (cd data/processed/{wildcards.name}/dbt \
                && wget {params.dbt_repo}/archive/refs/tags/v{params.release}.zip && unzip v{params.release}.zip
            ) &>> {log}
            mv data/processed/{wildcards.name}/dbt/bgcflow_dbt-duckdb-{params.release} {params.dbt}
            rm data/processed/{wildcards.name}/dbt/v{params.release}.zip
        fi

        python {params.dbt}/scripts/source_template.py {params.dbt}/templates/_sources.yml {output.profile} {params.as_version} {params.cutoff} &>> {log}
        """

def exclude_model_dbt(model_to_ignore):
    """
    Returns a string containing the `--exclude` option followed by the models to ignore in a dbt project.

    Args:
        model_to_ignore (list of str): A list of model names to ignore.

    Returns:
        str: A string containing the `--exclude` option followed by the model names, or an empty string if the list is empty.

    Example:
        >>> exclude_model_dbt(['model1', 'model2'])
        '--exclude model1 model2'

    Note:
        The returned string can be used as an argument to the `dbt build` command to exclude the specified models from the build process.
    """
    if len(model_to_ignore) == 0:
        return ""
    else:
        return " ".join(["--exclude"] + model_to_ignore)

rule build_database:
    input:
        database="data/processed/{name}/antismash_database/antiSMASH_database_{version}",
        profile = "data/processed/{name}/dbt/antiSMASH_{version}/models/sources.yml"
    output:
        duckdb = "data/processed/{name}/dbt/antiSMASH_{version}/dbt_bgcflow.duckdb"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "logs/database/report/database_{version}_{name}.log"
    threads: 16
    params:
        dbt = "data/processed/{name}/dbt/antiSMASH_{version}",
        exclude = lambda wildcards: exclude_model_dbt(models_to_ignore[wildcards.name])
    shell:
        """
        cp {input.database}/antismash_db.duckdb {output.duckdb} 2>> {log}
        command="dbt build --threads {threads} {params.exclude} -x"
        echo $command >> {log}
        (cd {params.dbt} \
            && dbt debug \
            && $command \
        ) &>> {log}
        """
