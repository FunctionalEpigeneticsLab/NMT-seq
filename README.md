# Nextflow pipeline for data processing of scNMT-seq and sciMET
The cDNA and gDNA library processing are independent pipelines located in the folders `rna_pipeline` and `dna_pipeline`. The gDNA pipeline includes processing for both single-cell/single-nucleus NMT-seq and sciMET, as they share major processes.

### Prerequisites
- Ensure that `Nextflow` is installed and accessible in your system's PATH.
- Place `main.nf` and `nextflow.config` in your run folder e.g. a directory for processing data from a specific experiment batch.
- To build Singularity image: use the `def` files in the `singularity` folder.
- Modify `nextflow.config`: adjust the parameters listed to suit your own need; change the reference genome files to their location; set `singularity.cacheDir` to your own scratch folder with sufficient space; replace `accessToken` with your personal Nextflow Tower token for workflow monitoring; change the Singularity image (sif file) locations correspondingly to where your containers have been built and stored.

### Running the pipeline
For NMT-seq (use `rna_pipeline` and `dna_pipeline` for cDNA and gDNA libraries, respectively):
```shell
nohup nextflow run main.nf > pipeline.log 2>&1 &
```

For sciMET (`dna_pipeline`):
```shell
nextflow run main.nf -entry pipeline_sciMET > pipeline.log 2>&1 &
```

### Key steps and parameters of data processing
Overview of workflows:
![nmt_flowchart](https://github.com/user-attachments/assets/41414d4f-7b92-4596-930a-cd9af254e69e)
#### 1. cDNA (Smart-seq2) libraries
- **Raw reads QC**: generate `FastQC` & `MultiQC` report. Set `run_fastqc=false` if you do not need this.

- **Trimming**: trim off adapters and low-quality bases from the 3' ends by `TrimGalore`, and filter by read length of 36 bp.

- **Alignment**: generate genome and transcriptome alignment BAM files by `STAR`.

- **Filtering**: filter aligned reads by `Samtools`: unmapped reads and secondary alignments are filtered out, and only uniquely mapped reads are kept.

- **Gene expression quantification**:
  - for single-cell data, specify `rsem = true` in `nextflow.config` to use RSEM for quantification on the gene level based on the transcriptome alignment BAM files; this only counts exonic mapped reads, and accounts for bias in the read start position distribution.
  - for single-nucleus data, specify `featurecount_gene = true` in `nextflow.config` to use `featureCounts` for quantification on the gene level.

- **To generate a comprehensive summary file**:
  ```shell
  python rna_pipeline/bin/collect.RNAalign.nf.py
  ```

#### 2. gDNA libraries (DNA methylation & chromatin accessibility)
- **Raw reads QC**: generate `FastQC` & `MultiQC` report. Set `run_fastqc=false` if you do not need this.

- **Trimming of scNMT-seq**: use `TrimGalore` to perform adapter and base quality trimming, hard-clip 20 bp from the 5' end and 7 bp from the 3' end, and filter by read length of 36 bp. 
  
- **Trimming of sciMET**: use `TrimGalore` to perform adapter and base quality trimming, and filter by read length. 

- **Alignment**: alignment with `Bismark` that calls `Bowtie2` in `--local` (allowing for soft-clipping) `--non_directional` (essential for scNMT-seq PBAT libraries) mode.
  
- **Deduplication**: `deduplicate_bismark` is used to remove duplicates.

- **Filtering**: filter aligned reads by `Samtools`: only uniquely mapped reads are kept; reads with mapping quality (MAPQ) score below 30 are filtered out, and only uniquely mapped reads are kept.

- **Downsampling reads from BAM files**: randomly sample 10%, 20%, 40%, 60%, and 80% of reads from the Bismark output BAM files by `Samtools`. Afterwards, deduplication, filtering, and computing coverage are performed on the downsampled BAM files. This is for estimating a sequencing saturation curve for each single-cell library.

- **Computing coverage**: compute genome-wide coverage of individual BAM files using `Samtools`.

- **Generating BED files**: convert BAM files to BED files that indicate read coverage positions using `BEDTools`. The BED files can be used for copy-number calling on the single-cell level with the software of your choice. Set `generate_bedfile=false` if this is not needed. 

- **Context-specific methylation calling**: 
  - Option 1: Enable `bismark_extract=true` to use `bismark_methylation_extractor` to assess the methylation status of each cytosine, and `coverage2cytosine` to create a genome-wide cytosine report. The `--nome-seq` option is used for NMT-seq contexts. 
  - Option 2: Enable `allc_extract=true` to use `Methylpy` for methylation calling to convert BAM files to ALLC format. We include the upstream 1 base to take GpC context into account for NMT-seq. Extraction of different cytosine contexts was then carried out by `extract-allc` module of `ALLCools`. The single-cell ALLC files can be directly used as input for downsteam analysis. 
    - for scNMT-seq: we use `--mc_contexts WCGN` to extract CpG sites (reflecting DNA methylation), and `--mc_contexts GCHN` for GpC sites (reflecting chromatin accessibility due to GpC methyltransferase treatment);
    - for sciMET: we use `--mc_contexts NCGN` to extract CG sites, and `--mc_contexts NCHN` to extract CH sites.

- **To generate a comprehensive summary file**:
  - For NMT-seq:
    ```shell
    python dna_pipeline/bin/generate_allc_summary.py --input_dir allc --out summary/allc_summary.tsv
    python dna_pipeline/bin/generate_mapping_summary.NMT.py
    ```
  - For sciMET:
    ```shell
    python dna_pipeline/bin/count_barcodes_from_fq.sciMET.py --input_dir trim --out summary/trimmed_reads_summary.all_barcodes.tsv &
    sh dna_pipeline/bin/concat_mapping_summary.sciMET.sh bam_sc summary/uniq_mapped_reads_summary.all_barcodes.tsv &
    python dna_pipeline/bin/generate_allc_summary.py --input_dir allc --out summary/allc_summary.tsv &
    wait
    python dna_pipeline/bin/generate_metadata.allCells.py --input_dir summary --out summary/metadata.allCells.csv
    ```
