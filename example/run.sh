#!/bin/zsh

sample_sheet="example/sample_sheet.csv"
metadata="./example/metadata_splited.csv"
embedding="./example/projection.csv"
delimiter="__"

nextflow run main.nf \
    --sample_sheet "$sample_sheet" \
    --metadata "$metadata" \
    --delimiter "__" \
    --method "umap" \
    --vis_label "celltype" \
    --embedding $embedding \
    --embedding_name "tsne"