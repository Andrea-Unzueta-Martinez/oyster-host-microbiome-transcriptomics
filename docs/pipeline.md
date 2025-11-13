**#Pipeline: QC → counts (host and microbiome transcriptome)**

## Summary:
This document outlines the workflow used to process raw RNA-seq reads from oysters and their associated microbiome, producing gene-level count tables for downstream differential expression and co-expression analyses. All steps were performed on a high-performance computing cluster (HPC) running a Linux environment.

## Data source
Raw reads are available under NCBI BioProject accession [PRJNA1313129](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1313129).


## 1. Input data
- **Sequencing type:** Paired-end Illumina NovaSeq S4, 2 × 150 bp reads  
- **Flow cell type:** S4 (10B)  
- **Sequencing yield:** Two lanes were used, producing approximately 5 billion paired-end reads in total. This corresponds to an average of ~125–166 million read pairs per library after demultiplexing.  
- **Samples:** 24 biological samples, each representing total RNA from an individual oyster, containing both host and microbiome sequences  
- **Controls:** 4 procedural RNA extraction negatives and 2 mock community standards  
- **Input files:** 30 total paired-end libraries (`*.fastq.gz`), demultiplexed per sample  
- **Reference assemblies:**  
  - **Host:** *Crassostrea virginica* reference genome, NCBI accession [GCF_002022765.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002022765.2/)  
  - **Microbiome:** de novo metagenome assembly generated from pooled reads across all samples, with predicted and annotated open reading frames (ORFs) used as reference for mapping and counting


## 2. Quality control

### Tools
- **Trim Galore! v0.6.10** — adapter and quality trimming (uses Cutadapt + FastQC)  
- **FastQC** — per-sample quality assessment  
- **MultiQC** — aggregate summary of QC metrics  
- **Conda environment:** `metagenomics` (Python 3.x) with MultiQC installed  

### Workflow

#### 1. Environment setup
```bash
ml python                # Load Python module on cluster
conda activate metagenomics  # Activate environment with MultiQC installed
```

#### 2. Adapter trimming and quality filtering
Trim Galore was run on all paired-end libraries to remove adapters, discard short reads, and generate FastQC reports:

```bash
trim_galore \
  --paired \
  --length 150 \
  --fastqc_args "--nogroup --noextract" \
  --output_dir TrimGalore_output/ \
  /path/to/raw_reads/*_R1_001.fastq /path/to/raw_reads/*_R2_001.fastq
```

- Reads shorter than 150 bp after trimming were removed.  
- Adapter sequences were automatically detected and removed.  
- FastQC reports were generated for all libraries (samples, negative controls, and mock communities).

#### 3. Aggregate QC reports
A combined summary report was generated using MultiQC:

```bash
cd TrimGalore_output/
multiqc . -o qc_reports/
```

The resulting `multiqc_report.html` summarizes quality metrics (per-base quality, GC content, sequence duplication, and adapter content) across all libraries.

#### 4. File management
Quality-filtered reads and QC reports were transferred to long-term storage and local machines using secure copy (`scp`) or `rsync` commands:

```bash
scp /path/to/cluster/project/qc_reports/multiqc_report.html /path/to/local/storage/
```

### Notes
- QC was performed on all **24 biological samples**, **4 extraction negatives**, and **2 mock communities**.  
- Quality filtering confirmed consistent base quality and read-length distributions across libraries.  

## 3. Removal of ribosomal reads

### Overview
Ribosomal RNA (rRNA) reads were removed prior to downstream mapping and quantification to retain only mRNA-derived sequences for both host and microbiome analyses.  

### Tools and databases
- **SortMeRNA v4.3.6** — rRNA read removal  
- **Reference databases:** recommended SortMeRNA database set from release 4.3.4  
  - `smr_v4.3_default_db.fasta` (recommended default database which includes  Silva 138 SSURef NR99 (16S, 18S), Silva 132 LSURef (23S, 28S), and RFAM v14.1 (5S, 5.8S))
  - Databases were obtained from the official SortMeRNA release package:  
  [https://github.com/biocore/sortmerna/releases/tag/v4.3.4](https://github.com/biocore/sortmerna/releases/tag/v4.3.4)

### Workflow

#### 1. rRNA removal using SortMeRNA
SortMeRNA was used to identify and remove ribosomal RNA sequences from all paired-end libraries prior to mapping.  
Each sample was processed in parallel on the cluster using a SLURM batch submission script.

- Memory allocation: **500 GB**  
- Threads per job: **8**

Example command:

```bash
sortmerna \
  --idx-dir /path/to/rRNA_databases_v4 \
  -ref /path/to/rRNA_databases_v4/smr_v4.3_default_db.fasta \
  -reads /path/to/quality_filtered_reads/${sample}_R1_001_val_1.fq \
  -reads /path/to/quality_filtered_reads/${sample}_R2_001_val_2.fq \
  --aligned ${sample}_rRNA \
  --other ${sample}_clean \
  --task 4 \
  --threads 8 \
  -v \
  --paired_in \
  --fastx \
  --out2 \
  --workdir /path/to/rRNAremoval_outputs/${sample}_outputdir
```

#### 2. Outputs
- **`--aligned`** output contained rRNA reads.  
- **`--other`** output contained non-rRNA (clean) reads used in downstream host and microbiome quantification.  
- Log files for each sample were stored in the designated output directories for quality verification.

### Notes
- rRNA removal was applied to all **24 biological samples**, **4 extraction negatives**, and **2 mock communities**.  
- Effective rRNA depletion was confirmed based on read retention and SortMeRNA summary reports.  
- The resulting “clean” reads served as input for mapping and quantification.


## 4. Splice-aware mapping to the host genome

### Overview
Clean, non-rRNA reads were aligned to the *Crassostrea virginica* reference genome to separate host-derived sequences from microbial reads and to generate host gene quantification files.  
Mapping was performed with the splice-aware aligner **STAR**, optimized for large eukaryotic transcriptomes, on a high-performance computing cluster (HPC) running a Linux environment.

### Tools and reference
- **STAR v2.7.11a** — splice-aware aligner for RNA-seq data  
- **Reference genome:** *Crassostrea virginica* (NCBI accession [GCF_002022765.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002022765.2/))  
- **Index:** pre-built STAR genome index generated from the reference genome and GTF annotation file using default parameters  
- **Environment:** `metagenomics` conda environment with STAR installed

### Workflow

#### 1. Parallel alignment of clean reads
STAR was executed as a SLURM array job to process multiple samples in parallel.

**Cluster resource settings**
- Memory allocation: **500 GB**  
- Threads per task: **6**  
- Wall time: **48 hours**  
- Scheduler: SLURM array (`--array=1-27`)

**Example SLURM submission script**

```bash
STAR \
  --genomeDir /path/to/cvirginica_genome_indices/ \
  --runThreadN 6 \
  --readFilesIn /path/to/rRNAremoval_outputs_clean/${sample}_clean_fwd.fq \
                 /path/to/rRNAremoval_outputs_clean/${sample}_clean_rev.fq \
  --outFileNamePrefix /path/to/STAR_Alignment/${sample}_alignment_results_ \
  --outSAMunmapped Within KeepPairs \
  --outReadsUnmapped Fastx \
  --outSAMstrandField intronMotif
```

- `--outSAMunmapped Within KeepPairs` retained unmapped read pairs for downstream microbiome analysis.  
- `--outReadsUnmapped Fastx` exported unmapped reads in FASTQ format.  
- `--outSAMstrandField intronMotif` ensured proper intron motif annotation for downstream counting tools.

#### 2. Outputs
- **Mapped host reads:** BAM/SAM files for each sample located in `STAR_Alignment/`  
- **Unmapped reads:** FASTQ files containing reads not aligned to the host genome; these were used as input for microbiome quantification  
- **Log files:** STAR alignment summaries (`Log.final.out`) stored per sample and used to calculate mapping efficiency

### Notes
- Host alignment rates and read distributions were verified using STAR’s `Log.final.out` files.  
- Unmapped read pairs were retained for microbiome mapping to avoid loss of non-host transcripts.


## 5. De novo assembly and filtering of microbiome transcripts

### Overview
Reads that did not align to the *Crassostrea virginica* reference genome (obtained as unmapped reads from STAR) were used to assemble the microbiome transcriptome.  
Assembly was performed de novo using **rnaSPAdes**, followed by **BLASTn** filtering to remove any residual host contamination.  
All steps were executed on a high-performance computing cluster (HPC) running a Linux environment.

### Tools
- **rnaSPAdes v3.15.5** — de novo transcriptome assembler  
- **BLAST+ v2.13.0** — sequence similarity search for host filtering  
- **makeblastdb** — creation of nucleotide BLAST databases  
- **Python 3.10** — custom filtering script for recovering contigs with no BLAST hits  
- **Reference genome for filtering:** *Crassostrea virginica* (NCBI accession [GCF_002022765.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002022765.2/))

### Workflow

#### 1. Prepare unmapped reads for assembly
Unmapped reads from STAR (both mates) were located in the `STAR_Alignment` directory.  
To ensure compatibility with rnaSPAdes, `.fastq` extensions were added to all files:

```bash
for f in *Unmapped.out.mate1; do mv "$f" "$f.fastq"; done
for f in *Unmapped.out.mate2; do mv "$f" "$f.fastq"; done
```

#### 2. De novo assembly with rnaSPAdes
rnaSPAdes was used to assemble transcripts from the unmapped reads, representing the combined microbiome of each oyster sample.

```bash
rnaspades.py \
  --pe1-1 sample_Unmapped.out.mate1.fastq \
  --pe1-2 sample_Unmapped.out.mate2.fastq \
  -o rnaSPAdes_output/sample_directory/ \
  -t 8 \
  -m 500
```

#### 3. Remove host contamination using BLASTn
To ensure that only microbiome-derived transcripts were retained, each assembled transcriptome was compared to the *C. virginica* genome using BLASTn.

**Create host BLAST database**
```bash
makeblastdb \
  -in /path/to/C_virginica_genomic.fna \
  -dbtype nucl \
  -out cvirginica_blastdb
```

**Run BLASTn against the host genome**
```bash
blastn \
  -query rnaSPAdes_output/sample_directory/transcripts.fasta \
  -db cvirginica_blastdb \
  -outfmt 6 \
  -evalue 1e-10 \
  -num_threads 8 \
  -out sample_vs_cvirginica_blast_results.txt
```

#### 4. Filter assembled transcripts with a custom Python script
A Python script was used to recover contigs that produced **no BLAST hits** to the oyster genome and therefore represent putative microbiome transcripts.

```python
# Function to recover hits that didn’t align to the oyster genome
def filter_assemblies_blast(blastFile):
    filt = []
    with open(blastFile) as inFile:
        for line in inFile:
            if "# Query:" in line:
                q1 = line.split()[2]
            if "hits found" in line:
                n = int(line.split()[1])
                if n == 0:
                    filt.append(q1)
    return filt

# Function to extract FASTA sequences from filtered hits
def get_fasta_seqs(hits, input_fasta_file, out_file, wd):
    seq = ""
    read = 0
    with open(f"{wd}/{out_file}", "w") as outFile:
        with open(input_fasta_file) as inFile:
            for line in inFile:
                if line[0] == ">":
                    if seq != "":
                        outFile.write(f">{h}\n")
                        outFile.write(f"{seq}\n")
                        seq = ""
                        read = 0
                    h = line.strip()[1:]
                    if h in hits:
                        read = 1
                else:
                    if read:
                        seq += line.strip()
        if seq != "":
            outFile.write(f">{h}\n")

# Example usage (update file paths for each sample)
wd = "/path/to/BlastMicrobiomeAssembledTranscripts"
blastFile = "/path/to/sample_transcripts_vs_cvirginica.blast"
fastaFile = "/path/to/sample_rnaSPAdes_output/transcripts.fasta"
hits = filter_assemblies_blast(blastFile)
get_fasta_seqs(hits, fastaFile, "sample_transcripts_filtered.fasta", wd)
```

The script first identifies transcripts with **no BLAST hits** to the host genome and then extracts those sequences from the rnaSPAdes output FASTA file into a new filtered file.

#### 5. Outputs
- `rnaSPAdes_output/` — per-sample transcriptome assemblies  
- `blast_results/` — BLASTn output tables (tabular format, outfmt 6)  
- `filtered_transcripts/` — FASTA files of microbiome-only transcripts recovered by the Python script  

### Notes
- The custom filtering script provided a transparent, reproducible way to identify non-host contigs.  
- BLAST summary logs confirmed that most host contamination was removed after filtering.  
- The resulting filtered transcriptomes served as references for downstream microbiome quantification and annotation.
- **negative controls failed de novo assembly** due to extremely low read depth and were excluded from further analysis.


## 6. Quantification of microbiome transcripts

### Overview
To estimate transcript abundance for the microbiome, non-oyster (unmapped) reads were mapped back to the filtered microbiome transcript assemblies.  
This step was performed to quantify transcript-level expression within the microbial community.  
ORFs were first predicted from the filtered assemblies using **Prodigal**, and relative abundance was estimated by mapping reads to those ORFs with **Bowtie2**.

### Tools
- **Prodigal v2.6.3** — identification of open reading frames (ORFs) in assembled transcripts  
- **Bowtie2 v2.5.1** — short-read alignment and quantification  
- **Samtools v1.17** — downstream file handling (conversion, sorting, indexing)  
- **Environment:** `metagenomics` conda environment containing Prodigal, Bowtie2, and Samtools  

### Workflow

#### 1. Predict open reading frames (ORFs)
Prodigal was used in metagenomic mode (`-p meta`) to identify coding sequences from each filtered microbiome transcriptome.

```bash
prodigal \
  -i ${sample}_rnaSPAdes_output_filt.fasta \
  -d ${sample}_non_oyster_transcript_orfs.fa \
  -p meta \
  -o ${sample}_non_oyster_transcripts.prodigal
```

- The output `*.fa` file contains the predicted nucleotide sequences of ORFs for each sample.  
- Prodigal was run on each filtered assembly in parallel using SLURM batch submission.

#### 2. Build Bowtie2 reference indices
Each sample-specific ORF FASTA file was indexed using Bowtie2 to create mapping references.

```bash
bowtie2-build \
  ${sample}_non_oyster_transcript_orfs.fa \
  ${sample}_non_oyster_transcript_orfs.reference
```

- Reference indices were stored in a dedicated directory for subsequent quantification.  
- Indexing was performed once per sample, using the corresponding ORF FASTA file.

#### 3. Map reads to microbiome ORFs
Clean, non-oyster reads (unmapped reads from STAR) were aligned back to each ORF reference to estimate abundance.

```bash
bowtie2 \
  -x /path/to/NonOysterFunctionalAnnotations/COUNTS/${sample}_non_oyster_transcript_orfs.reference \
  -1 ${sample}_alignment_results_Unmapped.out.mate1.fastq \
  -2 ${sample}_alignment_results_Unmapped.out.mate2.fastq \
  -S /path/to/NonOysterFunctionalAnnotations/COUNTS/${sample}.sam \
  --very-sensitive-local \
  -p ${SLURM_CPUS_PER_TASK} \
  --verbose
```

- The `--very-sensitive-local` mode was used to maximize read alignment sensitivity across diverse microbial transcripts.  
- Output SAM files were generated per sample and later converted to BAM for downstream processing.  
- Mapping statistics were reviewed to confirm consistent alignment rates across samples.

#### 4. Outputs
- `*_non_oyster_transcript_orfs.fa` — nucleotide sequences of predicted ORFs  
- `*_non_oyster_transcript_orfs.reference*` — Bowtie2 index files  
- `*.sam` — alignment files of non-oyster reads mapped to sample-specific ORFs  
- `*.bam` — optional, sorted binary alignment files (converted with Samtools)  
- Summary mapping statistics for each sample

### Notes
- This quantification step produced per-sample abundance profiles for microbiome genes and transcripts.  
- The resulting SAM/BAM files were used to generate combined count matrices for downstream analyses (e.g., differential expression and WGCNA).  
- Mapping results confirmed that non-host reads primarily aligned to microbial ORFs, validating the de novo assembly and host-filtering steps.


## 7. Build microbiome count matrix from per-sample alignments

### Overview
Per-sample alignments of non-oyster reads to sample-specific ORFs (Section 6) were summarized into **per-sample count tables** (rows = ORFs/contigs for that sample; one CSV per sample).  
**Important:** Because ORFs were assembled independently per sample, ORF IDs are **not shared across samples**, so a naïve union join on ORF IDs is not appropriate.

### Script
- **Path:** `scripts/make_microbiome_counts.py`  
- **Language:** Python (pandas required)  
- **Behavior:** For each sample, parse its SAM file and output a CSV with:
  - `<ref>_id` (e.g., `orf_id`)
  - `<sample>_counts` (unique full-length matches; CIGAR = `{read_length}M`)

### Inputs
- `sam_files.txt` — one relative path per SAM file, e.g.:
  ```
  STAR_Alignment/S01.sam
  STAR_Alignment/S02.sam
  ...
  ```
- SAM files from Section 6 (Bowtie2 `--very-sensitive-local`).

### Outputs
- Per-sample CSVs written alongside each SAM, e.g.:
  - `S01_ORF_counts.csv`
  - `S02_ORF_counts.csv`
  - …

### Example usage
```bash
# Ensure sam_files.txt contains one path per SAM file (no quotes)
python - <<'PY'
import os
from scripts.make_microbiome_counts import get_counts

sam_files, sample_names = [], []
with open("sam_files.txt") as f:
    for line in f:
        sam = line.strip()
        if not sam: 
            continue
        base = os.path.basename(sam)
        sample = base.split("_")[-2]   # adapt to your naming convention
        sam_files.append(sam)
        sample_names.append(sample)

for sam_file, sample in zip(sam_files, sample_names):
    out_csv = f"{sample}_ORF_counts.csv"
    get_counts(sam_file=sam_file, sample_name=sample, ref="orf", out=out_csv)
    print(f"[OK] {sample}: wrote {out_csv}")
PY
```

### ⚠️ NOTE on cross-sample merging
- ORF IDs are **sample-specific** (each assembly is independent), so feature names **do not match across samples**.
- Therefore, a simple “outer join on ORF IDs” is **not valid**.

### TODO — to be revised later
- **Planned approach:** Merge at the **functional level** in **R** using per-sample **functional annotation tables** (e.g., eggNOG/KEGG/COG) **and** per-sample ORF count tables, then aggregate counts by a shared functional hierarchy (e.g., KO → pathway, COG category, or custom ontology).
- **To add here later:**
  - R code to (1) read per-sample counts and functional tables, (2) map ORFs → function, (3) aggregate counts by function, (4) merge across samples on function identifiers.
  - Clear description of functional levels used (e.g., KO term, pathway/module).
  - Data dictionary for the function-level count matrix.

## 8. Functional annotation of microbiome ORFs (eggNOG-mapper v2)

### Overview
Functional annotations were generated per sample using **eggNOG-mapper v2**.  
For each filtered microbiome transcript FASTA (Section 5.4), annotations were computed using **DIAMOND** search mode.  
Runs were executed **one sample at a time** (parallelization attempts frequently failed on the HPC).

### Tools and databases
- **eggNOG-mapper v2** (`emapper.py`)
- **DIAMOND** for sequence similarity search (`-m diamond`)
- **Gene prediction inside emapper:** `--genepred prodigal` with `--itype metagenome`  
  *(eggNOG-mapper internally predicts ORFs from nucleotide inputs before annotation)*
- **eggNOG database**: pre-downloaded with `--data_dir /path/to/eggnog_data`  
  - See: https://github.com/eggnogdb/eggnog-mapper/wiki

> **Note:** If ORFs are **already** predicted externally (e.g., Prodigal in Section 6.1), you may pass amino-acid FASTA files with `--itype proteins` and **omit** `--genepred`. Here, following the executed workflow, `--genepred prodigal` was used so eggNOG-mapper handles gene prediction per sample.

### Workflow

#### 1. Prepare inputs
- Input per sample: filtered non-oyster transcript FASTA from Section 5.4  
  Example: `filtered_non_oyster_transcripts/SAMPLE_X_rnaSPAdes_output_filt.fasta`

#### 2. Run eggNOG-mapper per sample (serialized)
Run each sample independently to avoid HPC crashes observed with parallel submission.

```bash
emapper.py \
  -i /path/to/filtered_non_oyster_transcripts/SAMPLE_X_rnaSPAdes_output_filt.fasta \
  -o SAMPLE_X_non_oyster_transcripts_ \
  --genepred prodigal \
  --itype metagenome \
  -m diamond \
  --output_dir /path/to/NonOysterFunctionalAnnotations \
  --cpu 0 \
  --data_dir /path/to/NonOysterFunctionalAnnotations/eggnog_data
```

- `--cpu 0` uses **all available cores** on the node.
- `--genepred prodigal` + `--itype metagenome` instructs emapper to predict ORFs from nucleotide input before annotation.
- `-m diamond` provides fast homology search suitable for large metatranscriptomic inputs.

#### 3. Outputs (per sample)
eggNOG-mapper produces several files (prefix = `SAMPLE_X_non_oyster_transcripts_`), commonly including:
- `*.emapper.annotations` — main functional annotation table (orthologs, KOs/COGs/GO, etc.)
- `*.emapper.seed_orthologs` — seed ortholog hits
- `*.emapper.hits` — (optional, if enabled) detailed hit table
- `*.emapper.genepred.fasta` — predicted protein sequences from `--genepred prodigal` (AA FASTA)

### Notes
- Each sample was run **individually** to prevent job failures encountered when attempting parallel runs on the HPC.   
- The output `*.emapper.annotations` files contained predicted orthologs, KEGG ortholog (KO) IDs, COG categories, and GO terms used for downstream functional analyses.  
- Annotation tables were reformatted to tab-delimited files using the make_microbiome_counts.py script (available in the scripts/ directory) prior to import into R for downstream analyses.  
- Annotation tables were reformatted to tab-delimited files using the make_microbiome_counts.py script (available in the scripts/ directory) prior to import into R for downstream analyses.


## 9. Functional annotation of oyster transcripts (eggNOG-mapper v2)

### Overview
The *Crassostrea virginica* reference transcriptome was re-annotated using **eggNOG-mapper v2** to obtain consistent functional annotations across host and microbiome datasets.  
Annotations were performed at the **transcript isoform level** using the NCBI *C. virginica* transcript FASTA (`rna.fna`).

### Tools and databases
- **eggNOG-mapper v2.1.11** (`emapper.py`)  
  - Expected eggNOG DB version: **5.0.2**  
  - Installed eggNOG DB version: **5.0.2**  
  - DIAMOND version: **2.1.8**  
  - MMseqs2 version: **14.7e284**  
- **Gene prediction:** `--genepred prodigal` with `--itype genome`  
- **eggNOG database:** downloaded locally with `download_eggnog_data.py` and indexed using `create_dbs.py` (Metazoa subset)

### Workflow

#### 1. Download reference transcriptome
```bash
genome_updater.sh \
  -d refseq \
  -g invertebrate \
  -s "Crassostrea virginica" \
  -f rna.fna \
  -o c_virginica2
```

#### 2. Download and prepare eggNOG databases
```bash
download_eggnog_data.py --data_dir /path/to/eggnog_databases

create_dbs.py \
  -m diamond \
  --dbname eggnog_db_metazoa \
  --taxa Metazoa \
  --data_dir /path/to/eggnog_databases/
```

#### 3. Annotate oyster transcripts
```bash
emapper.py \
  -i /path/to/c_virginica2/files/rna.fna \
  -o oyster_transcriptome_annotations \
  --genepred prodigal \
  --itype genome \
  -m diamond \
  --output_dir /path/to/eggnog_output/ \
  --cpu 0 \
  --data_dir /path/to/eggnog_databases
```

#### 4. Outputs
- `oyster_transcriptome_annotations.emapper.annotations` — main functional annotation table (orthologs, KOs, COGs, GO terms)  
- `oyster_transcriptome_annotations.emapper.seed_orthologs` — seed ortholog hits  
- `oyster_transcriptome_annotations.emapper.genepred.fasta` — predicted protein sequences (AA FASTA)

### Notes
- Functional annotations were generated for the entire *C. virginica* transcriptome to ensure consistency between host and microbiome datasets.  
- Annotations are used to link host transcript counts (Section 10) to functional categories.  
- Annotation results were downloaded to a local machine for integration and visualization in R.  
- **Mock community samples were not processed** for these oyster steps because they represent microbial community standards, not host tissue.  
- Some **RNA extraction negative controls** were excluded from downstream host analyses because they contained insufficient oyster reads for annotation or quantification.


## 10. Quantification of oyster transcripts (Salmon)

### Overview
Oyster transcript abundance was quantified using **Salmon** in **mapping-based mode**, aligning oyster-mapped reads against the NCBI *C. virginica* reference transcriptome (`rna.fna`).  
Reads were extracted from the STAR alignments produced in Section 4 and converted to paired-end FASTQ files prior to quantification.

### Tools and reference
- **Salmon v1.10.3** — lightweight transcript quantifier  
- **Reference transcriptome:** NCBI *Crassostrea virginica* (`rna.fna`)  
- **Index:** built once with `salmon index`  
- **Input reads:** paired-end FASTQs containing only oyster reads extracted from STAR-mapped BAMs  
- **Environment:** `metagenomics` conda environment with Salmon installed

### Workflow

#### 1. Build Salmon index
```bash
salmon index \
  -t /path/to/c_virginica2/files/rna.fna \
  -i c_virginica_index
```

#### 2. Extract oyster-mapped reads from STAR BAMs
Convert STAR alignments to FASTQ files containing only mapped oyster reads.

```bash
# Sort BAMs and filter mapped reads
for f in *.sam; do
    filename="${f%%.*}"
    samtools sort "$f" -o "${filename}.sorted.bam"
    samtools view -b -F 4 "${filename}.sorted.bam" > "${filename}_mapped.bam"
done

# Collate paired reads and convert to FASTQ
for f in *_mapped.bam; do
    filename="${f%%.*}"
    samtools collate "$f" -o "${filename}.collate.bam"
    bamToFastq \
        -i "${filename}.collate.bam" \
        -fq "${filename}_mapped_R1.fq" \
        -fq2 "${filename}_mapped_R2.fq"
done
```

#### 3. Quantify transcript abundance with Salmon
```bash
salmon quant \
  -i /path/to/c_virginica_index \
  -l A \
  -1 ${sample}_mapped_R1.fq \
  -2 ${sample}_mapped_R2.fq \
  -p 15 \
  -o /path/to/Salmon_output/${sample}_quant
```

- `-l A` autodetects library type.  
- Output directories (`*_quant/`) contain `quant.sf` files for each sample.

#### 4. Outputs
- Per-sample directories: `Salmon_output/${sample}_quant/`  
  - `quant.sf` — transcript-level abundance estimates  
  - `cmd_info.json` — run parameters  
- Combined output will be used to generate host count matrices for downstream statistical analyses.

### Notes 
- **Mock community samples were not included** because they represent microbial standards rather than host transcriptomes.  
- **negative controls failed Salmon quantification** due to extremely low read depth and were excluded from further analysis.  
- Final outputs include transcript-level abundances (TPM and raw counts) matching the transcript identifiers annotated in Section 9.  
- Analyses in R (e.g., DESeq2, WGCNA) were performed at the transcript isoform level for improved biological resolution.
