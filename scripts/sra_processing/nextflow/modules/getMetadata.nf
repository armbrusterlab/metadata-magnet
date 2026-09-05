/*
 * Retrieves metadata for the SRA query specified.
 */
process getMetadata {
    conda "${workflow.projectDir}/envs/env_messy.yml"

    input:
    val sra_query

    output:
    path "metadata_esearch.csv", emit: metadata_esearch
    path "metadata_pysradb.tsv", emit: metadata_pysradb
    path "taxonomy_analysis/", emit: taxonomy_analysis

    script:
    """
    projDir="${workflow.projectDir}" # currently, scripts are located in the parent dir of projDir
    bash "\$projDir/../get_sra_metadata.sh" "${sra_query}"
    """
}