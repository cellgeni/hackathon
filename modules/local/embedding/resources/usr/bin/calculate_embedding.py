#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import scanpy as sc
import numpy as np


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    stream=sys.stdout,  # Direct output to stdout instead of a file
)


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Script validates sample and annotation tables and splits annotation table into separate celltypes"
    )
    parser.add_argument(
        "--h5ad_file",
        type=str,
        metavar="<path>",
        help="Specify a path to .h5ad input file",
        required=True,
    )
    parser.add_argument(
        "--sample_id",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample name for the file",
    )
    parser.add_argument(
        "--method",
        type=str,
        metavar="<str>",
        default="pca",
        help="Specify the method for dimensionality reduction: pca, umap, tsne. Default: 'pca'",
    )
    parser.add_argument(
        "--neighbors_key",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify the key for neighbors in AnnData object. Default: None (will compute neighbors if not present)",
    )
    parser.add_argument(
        "--output",
        metavar="<path>",
        type=str,
        help="Specify a path to output .h5ad file",
        default=10,
    )
    parser.add_argument(
        "--layer",
        type=str,
        metavar="<str>",
        default=None,
        help="If provided, which element of layers to use for PCA. Default: None (use .X)",
    )
    parser.add_argument(
        "--n_comps",
        type=str,
        metavar="<str>",
        default=50,
        help="Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation.",
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        metavar="<int>",
        default=15,
        help="Specify the number of neighbors for the neighbors graph. Default: 15",
    )
    parser.add_argument(
        "--vis_label",
        type=str,
        metavar="<str>",
        nargs="*",
        help="Specify the label for visualization in the embedding plot",
        default=None,
    )
    parser.add_argument(
        "--neighors_file",
        type=str,
        metavar="<path>",
        default="neighbors.npy",
        help="Specify the path to save neighbors information. Default: 'neighbors.npy'",
    )
    parser.add_argument(
        "--embedding_file",
        type=str,
        metavar="<path>",
        default="embedding.npy",
        help="Specify the path to save embedding information. Default: 'embedding.npy'",
    )

    return parser


def main():
    """
    Main function to process the input file and convert it to .h5ad format
    """
    parser = init_parser()
    args = parser.parse_args()

    # Load the input file
    logging.info(f"Loading input file: {args.h5ad_file}")
    adata = sc.read_h5ad(args.h5ad_file)

    # Calculate the embedding based on the specified method
    match args.method:
        case "pca":
            logging.info("Using PCA for dimensionality reduction")
            sc.tl.pca(adata, n_comps=args.n_comps, layer=args.layer)
        case "umap":
            logging.info("Using UMAP for dimensionality reduction")
            neighbors_key = (
                args.neighbors_key if args.neighbors_key else "pca_neighbors"
            )
            if args.neighbors_key not in adata.uns:
                logging.info("Computing neighbors for UMAP")
                sc.tl.pca(adata, n_comps=args.n_comps, layer=args.layer)
                sc.pp.neighbors(
                    adata, n_neighbors=args.n_neighbors, key_added=neighbors_key
                )
                adata.uns[neighbors_key]["connectivities_key"] = (
                    neighbors_key + "_connectivities"
                )
                adata.uns[neighbors_key]["distances_key"] = neighbors_key + "_distances"
            sc.tl.umap(adata, neighbors_key=neighbors_key, random_state=4)
        case "tsne":
            logging.info("Using t-SNE for dimensionality reduction")
            sc.tl.pca(adata, n_comps=args.n_comps, layer=args.layer)
            sc.tl.tsne(adata, use_rep="X_pca", random_state=4)
        case _:
            raise ValueError(f"Unknown method: {args.method}")

    # Plot the embedding if there is a label for visualization
    if args.vis_label:
        logging.info("Plotting the embedding")
        fig = sc.pl.embedding(
            adata,
            basis=f"X_{args.method}",
            color=args.vis_label,
            title=args.sample_id,
            return_fig=True,
        )
        fig.savefig(f"{args.method}_embedding.png", bbox_inches="tight")
        fig.savefig(f"{args.method}_embedding.pdf", bbox_inches="tight")

    # Save the AnnData object to the specified output file
    logging.info(f"Saving the AnnData object to {args.output}")
    adata.write_h5ad(args.output)

    # Save neighbors (if computed) and embedding information to separate numpy files
    if "pca_neighbors" in adata.uns:
        logging.info(f"Saving neighbors to {args.neighors_file}")
        np.save(args.neighors_file, adata.uns["pca_neighbors"])

    logging.info(f"Saving embedding to {args.embedding_file}")
    np.save(args.embedding_file, adata.obsm[f"X_{args.method}"])


if __name__ == "__main__":
    main()
