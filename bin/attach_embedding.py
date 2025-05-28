#!/usr/bin/env python3

import sys
import logging
import argparse
import scanpy as sc
import pandas as pd


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
        help="Specify a path to h5ad input file",
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
        "--embedding_file",
        type=str,
        metavar="<path>",
        help="Specify a path to a .csv file with embedding data and must contain 'barcode', 'dim1' and 'dim2' columns",
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
        "--barcode_column",
        type=str,
        metavar="<str>",
        default="obs_names",
        help='Specify barcode column name in AnnData object to match with metadata. Specify "obs_names" if you wish to match with AnnData object obs_names. Default: "obs_names"',
    )
    parser.add_argument(
        "--embedding_name",
        type=str,
        metavar="<str>",
        default="embedding",
        help="Specify the name for the embedding in AnnData object. Default: 'embedding'",
    )
    parser.add_argument(
        "--vis_label",
        type=str,
        metavar="<str>",
        nargs="*",
        help="Specify the label for visualization in the embedding plot",
        default=None,
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

    # Load the embedding data
    logging.info(f"Loading embedding data from: {args.embedding_file}")
    embedding_df = pd.read_csv(args.embedding_file)

    # Check that all barcodes in AnnData object are present in the embedding data
    if args.barcode_column == "obs_names":
        barcodes = adata.obs_names
    else:
        barcodes = adata.obs[args.barcode_column]
    missing_barcodes = set(barcodes) - set(embedding_df["barcode"])
    if missing_barcodes:
        logging.error(f"Missing barcodes in embedding data: {missing_barcodes}")
        sys.exit(1)

    # Merge the embedding data with the AnnData object
    logging.info("Attaching embedding to AnnData object")
    embedding_df.set_index("barcode", inplace=True)
    embedding = embedding_df.loc[barcodes, ["dim1", "dim2"]].values
    adata.obsm[f"X_{args.embedding_name}"] = embedding

    # Plot the embedding
    logging.info("Plotting the embedding")
    fig = sc.pl.embedding(
        adata,
        basis=f"X_{args.embedding_name}",
        color=args.vis_label,
        title=args.sample_id,
        return_fig=True,
    )
    fig.savefig(f"{args.embedding_name}_embedding.png", bbox_inches="tight")
    fig.savefig(f"{args.embedding_name}_embedding.pdf", bbox_inches="tight")

    # Save the AnnData object to the specified output file
    logging.info(f"Writing AnnData object to {args.output}")
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
