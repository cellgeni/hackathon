#!/bin/bash

sample_sheet="/lustre/scratch127/cellgen/cellgeni/aljes/hackathon/example/sample_sheet.csv"
metadata="/lustre/scratch127/cellgen/cellgeni/aljes/hackathon/example/metadata_splited.csv"
embedding="/lustre/scratch127/cellgen/cellgeni/aljes/hackathon/example/projection.csv"
delimiter="__"

nextflow run main.nf \
    --sample_sheet "$sample_sheet" \
    --metadata "$metadata" \
    --delimiter "__" \
    --method "umap" \
    --vis_label "celltype" \
    --embedding $embedding \
    --embedding_name "tsne" \
    -resume