# Emapper output to tsv file format
# -----------------------------------------------------------------------------
# This scripts converts an eggNOG-mapper *.emapper.annotations file into a clean, tab-delimited
# TSV suitable for downstream analysis (R/Python).
# -----------------------------------------------------------------------------

# Import pandas for tabular data manipulation
import pandas as pd  

# -----------------------------------------------------------------------------
# Function
# -----------------------------------------------------------------------------
def orf_to_csv(inputF, out):
    """
    Convert an eggNOG-mapper *.emapper.annotations file to TSV.

    Parameters
    ----------
    inputF : str
        Path to the input *.emapper.annotations file.
    out : str
        Path to the output *.tsv file.

    Processing rules
    ----------------
    - Lines starting with "##" are skipped (metadata).
    - The header line (starting with "#") is split on tabs and the first
      column name is set to "orf_id".
    - In data rows, any field equal to "-" is written as "NA".
    """
    with open(out,"w") as outFile:
        with open(inputF) as inFile:
            for line in inFile:
                # Skip extended metadata
                if line[0:2] != "##":
                    # Header line
                    if line[0] == "#":
                        hs = line.strip().split("\t")
                        hs[0] = "orf_id"
                        my_line = '\t'.join(hs)
                        outFile.write(f"{my_line}\n")
                    else:
                        # Data line: replace "-" with "NA"
                        toks = line.strip().split("\t")
                        toks_no_dash = []
                        for t in toks:
                            if t == "-":
                                toks_no_dash.append("NA")
                            else:
                                toks_no_dash.append(t)
                        my_line = "\t".join(toks_no_dash)
                        outFile.write(f"{my_line}\n")

# -----------------------------------------------------------------------------
# User-defined paths (edit these before running)
# -----------------------------------------------------------------------------
file = "/path/to/input/oyster_transcriptome_annotations.emapper.annotations"
out  = "/path/to/output/oyster_transcriptome_annotations.emapper.annotations.tsv"

# -----------------------------------------------------------------------------
# Execute
# -----------------------------------------------------------------------------
orf_to_csv(file, out)

# -----------------------------------------------------------------------------
# End of script
# -----------------------------------------------------------------------------
