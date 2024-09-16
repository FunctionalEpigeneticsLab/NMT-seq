# scNMT-seq pipeline

## QC and STAR alignment attributes of scRNA-seq data

### Output from STAR summary log file:

**Total_readpairs**: The total number of read pairs obtained after paired-end sequencing (two reads for each fragment). 

Extracted from STAR alignment output Log.final.out: "Number of input reads"

**Total_reads**: The total number of individual reads, both from the first and second ends of the fragments. This should typically be a bit higher than twice the number of Total_readpairs.

**Avg_length**: The average length (insert size) of the reads. 
*Extracted from STAR alignment output Log.final.out: "Average input read length"*

**Uniq_map_readpairs**: The number of read pairs that uniquely aligned to the reference genome. This means that these read pairs did not align to multiple locations but had a unique match.
Extracted from STAR alignment output Log.final.out: "Uniquely mapped reads number"*

**Uniq_map_rate**: The percentage of read pairs that uniquely aligned to the reference genome out of the total read pairs.
Extracted from STAR alignment output Log.final.out: "Uniquely mapped reads %"*
 
**Avg_map_length**: The average length of the aligned portion of the reads. This can be useful to assess the quality of the alignment and the fragmentation of the RNA.
Extracted from STAR alignment output Log.final.out: "Average mapped length"*

**Annotate_splices**: The number of spliced reads that have been annotated or matched to known splice junctions.
*Extracted from STAR alignment output Log.final.out***GT/AG_splices, GC/AG_splices, AT/AC_splices**: These are counts or percentages of observed splice junctions with specific dinucleotide sequences at the exon-intron boundaries. GT/AG is the canonical splice junction in eukaryotes, while GC/AG and AT/AC are considered non-canonical but still observed in a minor fraction of splicing events.

**Non-canonical_splices**: The count or percentage of spliced reads that do not match any of the canonical or known non-canonical splice junctions.
*Extracted from STAR alignment output Log.final.out*

**Mismatch_rate**: The percentage of bases in aligned reads that do not match the reference sequence. A high mismatch rate might indicate poor-quality reads or a divergent reference.
*Extracted from STAR alignment output Log.final.out*

**Multi_map_rate**: The percentage of reads that map to multiple locations in the genome or transcriptome.
*Extracted from STAR alignment output Log.final.out*

**Toomanyloci_rate**: The percentage of reads that map to an excessive number of locations, making them hard to assign to a specific location.
*Extracted from STAR alignment output Log.final.out*

**Unmap_mismatch_rate**: The percentage of reads that couldn't be mapped due to excessive mismatches.
*Extracted from STAR alignment output Log.final.out*

**Unmap_tooshort_rate**: The percentage of reads that couldn't be mapped because they were "too short". "Too short" relates to the best alignment mapped length, not the read length itself. Typically, most of the unmapped reads are classified as "too short". This could either happen for normal-length read which were not mapped well, or for a read that was actually too short (over-trimmed/low quality) before mapping.
*Extracted from STAR alignment output Log.final.out*

**Unmap_other_rate**: The percentage of reads that couldn't be mapped for reasons other than mismatches or "too short".
*Extracted from STAR alignment output Log.final.out*

**Chimeric_rate**: The percentage of reads that appear to be chimeric. Chimeric reads occur when one sequencing read aligns to two distinct portions of the genome with little or no overlap. Chimeric reads are indicative of structural variation, and they are also called split reads. Chimeric reads can also arise due to issues in library preparation or sequencing.
*Extracted from STAR alignment output Log.final.out*

### Output from featureCount using the default setting (count exonic reads)
**Exonic_map**: The number or percentage of reads that aligned to exonic regions of the genome. Exons are the coding regions of genes and are thus of particular interest in RNA-seq as they represent the expressed parts of genes.
*Extracted from featureCount.tsv.summary*

**Non-exonic_map**: The number or percentage of reads that aligned to non-exonic (intronic or intergenic) regions of the genome. While typically one would expect RNA-seq reads to map to exons, reads can also map to introns due to pre-mRNA or unspliced transcripts.
*Extracted from featureCount.tsv.summary*

**Gene_count**: The number of genes that had at least one read aligned to them. This can give an idea of the complexity and depth of the sequencing.
*Extracted from featureCount.tsv.summary*

### Output from samtools

**Uniq_map_ERCC (deprecated)**: If ERCC spike-ins (External RNA Controls Consortium) were used, this column would represent the number of reads that uniquely aligned to the ERCC controls. ERCC spike-ins are synthetic RNA molecules added to samples to provide a measure of technical variability.
*Output from samtools view Aligned.sortedByCoord.out.bam*

**ChrM_count (deprecated)**: The number of reads that aligned to the mitochondrial chromosome (ChrM). High levels of mitochondrial reads can indicate cellular stress or poor-quality cells.
*Output from samtools view Aligned.sortedByCoord.Processed.out.bam*

**Duplication_rate(incl_multimap)**: The percentage of reads that are considered duplicates, including those that map to multiple locations. Duplicates can arise due to PCR amplification during library preparation.
*Extracted from Aligned.sortedByCoord.Processed.out.bam.flagstat*

### Output from picard CollectRnaSeqMetrics
https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-

**MEDIAN_CV_COVERAGE**: The median coefficient of variation (CV) of coverage across the genome or transcriptome. CV is a measure of the relative variability and is calculated as the ratio of the standard deviation to the mean. A high CV indicates greater dispersion of coverage.
- The coefficient of variation (CV) is a measure of relative variability. In the context of sequencing coverage, it represents the dispersion of coverage values relative to the mean coverage. Specifically, \( CV = \frac{Standard \ Deviation}{Mean} \times 100 \).
- A smaller CV indicates that the coverage across different genes or regions is more consistent and uniform. This is generally desired because it means that there is even sequencing depth across the transcriptome, allowing for more reliable detection and quantification of transcripts.
- A higher CV indicates more variability in coverage, which can be due to various reasons such as library preparation biases, sequencing biases, or sample quality. High variability might mean that some regions are over-represented (possibly leading to wasted sequencing depth on those regions) while others are under-represented (which might hinder the detection of low-abundance transcripts).

**MEDIAN_5PRIME_BIAS and MEDIAN_3PRIME_BIAS**: These metrics measure the bias towards the 5' end or 3' end of transcripts. In some RNA-seq protocols, especially those that involve fragmentation, you might observe an over-representation of reads towards one end of the transcripts. This can be an indicator of degradation or specific library preparation biases.

**MEDIAN_5PRIME_TO_3PRIME_BIAS**: The comparison metric between the 5' and 3' biases mentioned above. A value close to 1 might indicate uniform coverage, whereas a deviation from 1 indicates a bias.

### Output from qualimap
Note: The reads are already adapter&quality trimmed. The attributes have ERCC reads included because ERCC sequences were added to the reference genome. 

**Uniq_map_rate_filtered**: the number of uniquely mapped read pairs after filtering out unmapped reads and secondary alignments
*Extracted from qualimap output: read pairs aligned*

**exonic_rate, intronic_rate, intergenic_rate**: 
*Extracted from qualimap output*
