process processMetadata {
    conda "${workflow.projectDir}/envs/metadata-magnet-env.yaml"

    input:
    path input_files
    path categoryFile
    path subcategoryFile
    path benchmark_file

    output:  
    path "*_processedMetadata.blast", emit: metadata
    path "benchmarking_afterMetadataProcessing.txt", emit: benchmark

    script:
    // inputBasename = input_file.simpleName

    // projDir="${workflow.projectDir}"
    // outname="${inputBasename}_processedMetadata.blast"

    // cp ${input_file} "\$outname"
    // bash "\$projDir/../scripts/metadata_processing_wrapper.sh" -f "\$outname" -c ${categoryFile} -s ${subcategoryFile}

// 

    // projDir="${workflow.projectDir}"
    // outname="${inputBasename}_processedMetadata.blast"

    // newBenchmarkFile="${inputBasename}_benchmarking_afterMetadataProcessing.txt"
    // cat ${benchmark_file} > \$newBenchmarkFile
    // echo "" >> \$newBenchmarkFile

    // cp ${input_file} "\$outname"
    // bash "\$projDir/../scripts/metadata_processing_wrapper.sh" -f "\$outname" -c ${categoryFile} -s ${subcategoryFile} >> \$ newBenchmarkFile
    """
    projDir="${workflow.projectDir}"

    newBenchmarkFile="benchmarking_afterMetadataProcessing.txt"
    cat ${benchmark_file} > \$newBenchmarkFile
    echo "" >> \$newBenchmarkFile

    # Process each metadata file
    for input_file in ${input_files}; do
        inputBasename=\$(basename "\$input_file" .blast)
        outname="\${inputBasename}_processedMetadata.blast"
        
        cp "\$input_file" "\$outname"
        bash "\$projDir/../scripts/metadata_processing_wrapper.sh" -f "\$outname" -c ${categoryFile} -s ${subcategoryFile} >> "\$newBenchmarkFile"
        echo "" >> "\$newBenchmarkFile"
    done
    """

    // stub:
    // inputBasename = input_file.simpleName
    // """
    // cat "${input_file}" > "${inputBasename}_processedMetadata.blast"
    // echo "metadata processing process" >> "${inputBasename}_processedMetadata.blast"
    // wc -l "${inputBasename}_processedMetadata.blast" >> "${inputBasename}_processedMetadata.blast" 
    // """
}

// process collectMetadataBenchmarking {
//     input:
//     path benchmark_file
//     path metadataBenchmarking

//     output:
//     path "benchmarking_metadata.txt"

//     script:
//     """
//     newBenchmarkFile="benchmarking_metadata.txt"
//     cat "${benchmark_file}" > \$newBenchmarkFile
//     echo "" >> \$newBenchmarkFile
//     cat ${metadataBenchmarking} >> \$newBenchmarkFile 
//     """
// }