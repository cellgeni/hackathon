#!/usr/bin/env python3

import sys
import logging
import argparse
import scanpy as sc


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
        "--output",
        metavar="<path>",
        type=str,
        help="Specify a path to output .h5ad file",
        default=10,
    )
    parser.add_argument(
        "--target_sum",
        type=float,
        metavar="<int>",
        default=1e4,
        help="Specify target sum for normalization. Default: 1e4",
    )
    parser.add_argument(
        "--hvg_flavor",
        type=str,
        metavar="<str>",
        default="seurat",
        help="Specify the flavor for highly variable genes detection. Default: 'seurat_v3'",
    )
    parser.add_argument(
        "--n_top_genes",
        type=int,
        metavar="<int>",
        default=2000,
        help="Specify the number of top variable genes to select. Default: 2000",
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
    adata.layers["counts"] = adata.X.copy()  # Store raw counts in a layer

    # Normalize the data
    logging.info("Normalizing the data")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    # Find highly variable genes
    logging.info("Finding highly variable genes")
    sc.pp.highly_variable_genes(
        adata, flavor=args.hvg_flavor, n_top_genes=args.n_top_genes
    )

    # Save the processed data to .h5ad format
    logging.info(f"Saving processed data to: {args.output}")
    adata.write_h5ad(args.output)
    logging.info("Processing complete.")


if __name__ == "__main__":
    main()
