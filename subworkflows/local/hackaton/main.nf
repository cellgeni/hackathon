include { TOH5AD } from '../../../modules/sanger/toh5ad'

process AttachCellMetadata {
    tag "Attaching cell metadata to .obs section of AnnData object for ${sample}"
    container "quay.io/cellgeni/metacells-python:latest"
    publishDir "${params.output_dir}/${sample}", mode: 'copy', overwrite: true
    input:
        tuple val(sample), path(h5ad, name: 'input/*')
        path(metadata)
        val(barcode_column)
    output:
        tuple val(sample), path("${sample}_meta.h5ad")
    script:
        """
        attach_annotation.py \
            --h5ad_file ${h5ad} \
            --sample_id ${sample} \
            --metadata ${metadata} \
            --barcode_column ${barcode_column} \
            --output ${sample}_meta.h5ad
        """
}

process ProcessGex {
    tag "Processing GEX data for ${sample}"
    container "quay.io/cellgeni/metacells-python:latest"
    publishDir "${params.output_dir}/${sample}", mode: 'copy', overwrite: true
    input:
        tuple val(sample), path(h5ad_file, name: 'input/*')
        val(target_sum)
        val(hvg_flavor)
        val(n_top_genes)
    output:
        tuple val(sample), path("${sample}_processed.h5ad")
    script:
        """
        process_gex.py \
            --h5ad_file ${h5ad_file} \
            ${target_sum ? "--target_sum ${target_sum}" : ""} \
            ${hvg_flavor ? "--hvg_flavor ${hvg_flavor}" : ""} \
            ${n_top_genes ? "--n_top_genes ${n_top_genes}" : ""} \
            --output ${sample}_processed.h5ad
        """
}

process Embedding {
    tag "Embedding GEX data for ${sample}"
    container "quay.io/cellgeni/metacells-python:latest"
    publishDir "${params.output_dir}/${sample}", mode: 'copy', overwrite: true
    input:
        tuple val(sample), path(h5ad_file, name: 'input/*')
        val(method)
        val(vis_label)
    output:
        tuple val(sample), path("${sample}_${method}.h5ad"), emit: h5ad
        path "*.npy", emit: embedding
        tuple path("*pdf"), path("*png"), emit: plots
    script:
        """
        calculate_embedding.py \
            --h5ad_file ${h5ad_file} \
            --sample_id ${sample} \
            --method ${method} \
            ${vis_label ? "--vis_label ${vis_label}" : ""} \
            --output ${sample}_${method}.h5ad \
            --neighors_file neighbors.npy \
            --embedding_file embedding.npy
        """
}

process AttachEmbedding {
    tag "Attaching embedding to .obs section of AnnData object for ${sample}"
    container "quay.io/cellgeni/metacells-python:latest"
    publishDir "${params.output_dir}/${sample}", mode: 'copy', overwrite: true
    input:
        tuple val(sample), path(h5ad_file, name: 'input/*')
        path(embedding)
        val(barcode_column)
        val(embedding_name)
        val(vis_label)
    output:
        tuple val(sample), path("${sample}*.h5ad"), emit: h5ad
        tuple path("*pdf"), path("*png"), emit: plots
    script:
        """
        attach_embedding.py \
            --h5ad_file ${h5ad_file} \
            --sample_id ${sample} \
            --embedding_file ${embedding} \
            --barcode_column "barcode" \
            ${embedding_name ? "--embedding_name ${embedding_name}" : ""} \
            ${vis_label ? "--vis_label ${vis_label}" : ""} \
            ${embedding_name ? "--output ${sample}_${embedding_name}.h5ad" : "--output ${sample}_atemb.h5ad"}
        """
}

process ConcatAnndata {
    tag "Concatenating AnnData objects for ${sample}"
    container "quay.io/cellgeni/metacells-python:latest"
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    input:
        tuple val(sample), path(anndata, name: 'input/*')
    output:
        path("combined.h5ad")
    script:
        """
        concat_adata.py \
            --anndata ${anndata} \
            --output combined.h5ad
        """
}

workflow HACKATON {

    take:
    files // channel: [ val(meta), [ file/dir ] ]
    metadata // channel: [ path(metadata.csv) ]
    embedding // channel: [ path(embedding.csv) ]
    delimiter // channel: [ val(delimiter) ]
    barcode_column // channel: [ val(barcode_column) ]
    target_sum // channel: [ val(target_sum) ]
    hvg_flavor // channel: [ val(hvg_flavor) ]
    n_top_genes // channel: [ val(n_top_genes) ]
    method // channel: [ val(method) ]
    vis_label // channel: [ val(vis_label) ]
    embedding_name // channel: [ val(embedding_name) ]

    main:

    h5ad_files = TOH5AD(
        files,
        delimiter ? delimiter : ""
        ).h5ad
    
    h5ad_files = h5ad_files.map { meta, h5ad ->
        tuple(meta.id, h5ad)
    }
    
    // Attach cell metadata to the .obs section of the AnnData object
    metadata_files = AttachCellMetadata(
        h5ad_files,
        metadata,
        barcode_column
    )

    // Process GEX data
    processed_gex_files = ProcessGex(
        metadata_files,
        target_sum ? target_sum : "",
        hvg_flavor ? hvg_flavor : "",
        n_top_genes ? n_top_genes : ""
    )

    // Embed GEX data
    embedded_files = Embedding(
        processed_gex_files,
        method,
        vis_label ? vis_label : ""
    ).h5ad

    // Attach embedding to the .obs section of the AnnData object
    tsne_files = AttachEmbedding(
        embedded_files,
        embedding,
        barcode_column,
        embedding_name ? embedding_name : "",
        vis_label ? vis_label : ""
    ).h5ad

    // Check if concatenation is needed
    tsne_files.toList()
              .branch {
                it ->
                combine_objects: it.size() > 1
                sample: true
              }
              .set { objects }

    // Concatenate AnnData objects
    concatenated_files = ConcatAnndata(
        objects.combine_objects.transpose().toList()
    )

    emit:
    // TODO nf-core: edit emitted channels
    h5ad      = AttachEmbedding.out.h5ad           // channel: [ val(meta), [ h5ad ] ]
    combined      = ConcatAnndata.out         // channel: [ val(meta), [ h5ad ] ]
}

