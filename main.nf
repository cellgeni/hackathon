include { HACKATON } from './subworkflows/local/hackaton'


workflow {
    // Read sample sheet and extract sample IDs and file paths
    files = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                   .splitCsv(header: true, sep: ',')
                   .map { row ->
                       tuple(row.sample_id, file(row.filename))
                   }

    // Read metadata and embedding files
    metadata = file(params.metadata)
    embedding = file(params.embedding)
    
    HACKATON(
        files,
        metadata,
        embedding,
        params.delimiter ? params.delimiter : "",
        params.barcode_column ? params.barcode_column : "",
        params.target_sum ? params.target_sum : "",
        params.hvg_flavor ? params.hvg_flavor : "",
        params.n_top_genes ? params.n_top_genes : "",
        params.method,
        params.vis_label ? params.vis_label : "",
        params.embedding_name ? params.embedding_name : "",
    )
}