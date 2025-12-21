# Count mapped reads to open reading frames (ORFs)
# -------------------------------------------------
# This script counts the number of reads mapping to each contig or ORF in
# SAM alignment files (e.g., Bowtie2 outputs). It can be used for either
# metagenomic or metatranscriptomic datasets.
# 
# The script reads a list of SAM files, extracts the number of unique reads
# mapped to each reference (contig or ORF), and writes one count table per sample.
# Each output file lists the reference IDs and corresponding read counts.
# -------------------------------------------------

# Import pandas for tabular data manipulation
import pandas as pd


def get_counts(sam_file, sample_name, ref, out):
    """
    Parse a SAM file and count the number of unique reads that map to each
    reference sequence (e.g., contig or ORF).

    Parameters
    ----------
    sam_file : str
        Path to the SAM file containing alignments for one sample.
    sample_name : str
        Sample identifier (used for naming the output column).
    ref : str
        Reference type (e.g., 'orf' or 'contig') used to label the first column.
    out : str
        Output filename for the resulting counts table (.csv).

    Returns
    -------
    pandas.DataFrame
        A dataframe with two columns: <ref>_id and <sample_name>_counts.
    """

    # Initialize dictionary to hold read IDs per contig/ORF
    contig_dict = {}

    # Open and parse the SAM file
    with open(f"{sam_file}") as inFile:
        for line in inFile:
            # Skip SAM header lines (which begin with "@")
            if line[0] != "@":
                toks = line.split("\t")

                # Extract relevant SAM fields
                contig_id = toks[2]   # Reference/contig name
                read_id = toks[0]     # Query/read name
                read = toks[9]        # Read sequence
                qual = toks[5]        # CIGAR string

                # Ignore unmapped reads (where contig_id == "*")
                if contig_id != "*":
                    # Define the expected perfect-match CIGAR string (e.g., "150M" for a 150 bp read)
                    target_qual = f"{len(read)}M"

                    # Only count full-length, perfect matches
                    if qual == target_qual:
                        # If the contig_id is not yet in the dictionary, initialize a list
                        if contig_dict.get(contig_id) is None:
                            contig_dict[contig_id] = []

                        # Append read_id to the list for this contig
                        contig_dict[contig_id].append(read_id)

    # Remove redundant read IDs (keep only unique reads per contig/ORF)
    for contig_id in contig_dict:
        contig_dict[contig_id] = len(list(set(contig_dict[contig_id])))

    # Convert the dictionary to a DataFrame
    df_counts = pd.DataFrame()
    df_counts[f"{ref}_id"] = contig_dict.keys()
    df_counts[f"{sample_name}_counts"] = contig_dict.values()

    # Sort results by descending read count
    df_counts = df_counts.sort_values(by=f"{sample_name}_counts", ascending=False)

    # Write output CSV in the same directory as the input SAM file
    out_wd = "/".join(sam_file.split("/")[0:-1])
    df_counts.to_csv(f"{out_wd}/{out}", index=False)

    return df_counts


# -------------------------------------------------
# Build lists of SAM files, output filenames, and sample names
# -------------------------------------------------
sam_files, out_files, sample_names = [], [], []

# Read the list of SAM files from an external text file
# (each line should contain the relative or absolute path to one SAM file)
with open("sam_files.txt") as inFile:
    for line in inFile:
        sam_file = line.strip()

        # Extract a short sample name from the file path by splitting on "_"
        out_base_name = sam_file.split("_")[-2]

        sample_names.append(out_base_name)
        sam_files.append(sam_file)

        # Define output filename (one CSV per sample)
        out_files.append(f"{out_base_name}_ORF_counts.csv")


# -------------------------------------------------
# Main loop: process each sample in order
# -------------------------------------------------

# Define working directory containing SAM files
wd = 'data/Transcriptomics/MicorbiomeTranscriptCounts'

# Loop through each file and count reads
for idx in range(len(sam_files)):
    get_counts(f"{wd}/{sam_files[idx]}", sample_names[idx], "orf", out_files[idx])
    print(f"sample: {sample_names[idx]} processed...")

# -------------------------------------------------
# End of script
# -------------------------------------------------
