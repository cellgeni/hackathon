process EMBEDDING {
    tag "Embedding GEX data for  $meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/cellgeni/metacells-python:latest':
        'quay.io/cellgeni/metacells-python:latest' }"

    input:
    tuple val(meta), path(h5ad_file, name: 'input/*')
    val(method)
    val(vis_label)

    output:
    tuple val(meta),  path("${meta.id}_${method}.h5ad"), emit: h5ad 
    tuple val(meta),  path("*.npy"), emit: npy
    tuple val(meta),  path("*.pdf"), emit: pdf
    tuple val(meta),  path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    script:
    """
    calculate_embedding.py \
        --h5ad_file ${h5ad_file} \
        --sample_id ${meta.id} \
        --method ${method} \
        ${vis_label ? "--vis_label ${vis_label}" : ""} \
        --output ${meta.id}_${method}.h5ad \
        --neighors_file neighbors.npy \
        --embedding_file embedding.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}_${method}.h5ad"
    touch "neighbors.npy"
    touch "embedding.npy"
    touch "${method}_embedding.png"
    touch "${method}_embedding.pdf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
    END_VERSIONS
    """
}