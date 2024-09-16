#!/usr/bin/env nextflow

/*
 * Author: Lab for Functional Epigenetics, KU Leuven
 * Date of update: 2024/08
 */

nextflow.enable.dsl=2
params.workingDir=file("file://${baseDir}")
workingDir=params.workingDir

log.info """\
	Single cell/nucleus NMT-seq cDNA pipeline V0
	============================================
	"""
	.stripIndent()


process FASTQC {
	tag "Performing FastQC..."
	publishDir "${workingDir}/FASTQC", mode: 'copy'
	
	input:
	path(reads)
	
	output:
	path "*.html", emit: htmls
	path "*.zip", emit: zips

	script:
	"""
	fastqc -t ${task.cpus} $reads
	"""
}

process MULTIFASTQC {
	tag "Performing MultiFastQC..."
	publishDir "${workingDir}/FASTQC", mode: 'copy'
	
	input:
	path "*"

	output:
	path "*.html"

	script:
	"""
	multiqc . -n multiqc_report_raw_data.html
	"""
}

process TRIM {
	tag "Trimming RNA reads..."
	publishDir "${workingDir}/TRIM", mode: 'copy'
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*_1.fq.gz"), path("*_2.fq.gz"), emit: trimmed_reads

	script:
	"""
	trim_galore ${params.extra_trim_params_RNA} --length ${params.trim_length_RNA} -paired ${reads[0]} ${reads[1]}
	"""
}

process FASTQC_TRIM {
    tag "Performing FastQC on trimmed reads..."
    publishDir "${workingDir}/FASTQC", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)

    output:
    path "*.html", emit: htmls
    path "*.zip", emit: zips

    script:
    """
    fastqc -t ${task.cpus} $reads
    """
}

process MULTIFASTQC_TRIM {
        tag "Performing MultiFastQC..."
        publishDir "${workingDir}/FASTQC", mode: 'copy'

        input:
        path "*"

        output:
        path "*.html"

        script:
        """
        multiqc . -n multiqc_report_trimmed_data.html
        """
}

process STAR_ALIGN {
	tag "Aligning using STAR..."
	publishDir "${workingDir}/STAR", mode: 'copy'
	
	input:
	path star_index
	path sjdbgtf
	tuple val(sample), path(file1), path(file2)

	output:
	tuple val(sample), path("*/*.Aligned.sortedByCoord.out.bam"), path("*/*.Aligned.toTranscriptome.out.bam"), path('*/*Log.final.out'), emit: general_output
	path "*/*Aligned.sortedByCoord.out.bam", emit: featurecount_output
	path "*/*Aligned.toTranscriptome.out.bam", emit: transcriptome_bam_file
	path "*/*Unmapped*", emit: unmap

	script:
	"""
	mkdir ${sample}
	/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN ${task.cpus} --runMode alignReads --genomeDir $star_index --sjdbGTFfile $sjdbgtf ${params.star} --readFilesIn $file1 $file2 --outFileNamePrefix ${sample}/${sample}.
	"""
}

process STAR_MARKDUP {
	tag "Marking duplicates for STAR bam files..."
	publishDir "${workingDir}/STAR", mode: 'copy'
	
	input:
	tuple val(sample), path(bam_file), path(trans_file), path(log)
	
	output:
	tuple val(sample), path("*/*.Processed.out.bam"), emit: general_output
	path '*/*.Processed.out.bam', emit: featurecount_output
	
	script:
	"""
	mkdir -p ${sample}
	/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN ${task.cpus} --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdentical --inputBAMfile $bam_file --outFileNamePrefix ${sample}/${sample}.Aligned.sortedByCoord.
	"""
}

process STAR_FILTER {
	tag "Keeping only the uniquely mapped reads..."
	publishDir "${workingDir}/STAR", mode: 'copy'
	
	input:
	tuple val(sample), path(bam_file)

	output:
	tuple val(sample), path("*/filtered.${bam_file}"), path("*/*.bai"), emit: general_output
	path "*/filtered.${bam_file}", emit: bam
	
	script:
	"""
	mkdir -p ${sample}
	samtools view -bS ${params.samtools_filter_star} $bam_file > tmp
	mv tmp ${sample}/filtered.$bam_file
	samtools index ${sample}/filtered.$bam_file
	"""
}

process STAR_FEATURECOUNTS {
	tag "Feature counting for exonic reads..."
	publishDir "${workingDir}/FEATURECOUNTS" , mode: 'copy'
	
	input:
	path '*'
	path gtf

	output:
	path "*featureCount.tsv*"

	script:
	"""
	featureCounts -p --countReadPairs --primary -T ${task.cpus} -g gene_id -a $gtf -o featureCount.tsv *.bam
	cat featureCount.tsv | tail -n +2 > tmp.tsv
	mv tmp.tsv featureCount.tsv
	"""
}

process STAR_FEATURECOUNTS_GENE {
        tag "Quantification at gene level..."
        publishDir "${workingDir}/FEATURECOUNTS", mode: 'copy'

        input:
        path '*'
        path gtf

        output:
        path "*featureCount.gene.tsv*"

        script:
        """
        featureCounts -t gene -p --primary -T ${task.cpus} -g gene_id -a $gtf -o featureCount.gene.tsv *.bam
	"""
}

process KB_COUNT {
	tag "kb count..."
	publishDir "${workingDir}/kb_count", mode: 'copy'

	input:
	path '*'

	output:
	path 'result/*'

	script:
	"""
	sh create_kb_batch.sh
	kb count -i ${params.kb_index} -g ${params.kb_t2g} -c1 ${params.kb_cdna_t2c} -c2 ${params.kb_intron_t2c} -x SMARTSEQ2 --parity paired -o result -t ${task.cpus} --sum=nucleus -w NONE --batch-barcodes batch.txt --workflow=nac --h5ad
	"""
}

process QUALITY_METRICS {
	tag "Extracting quality metrics..."
	publishDir "${workingDir}/STAR", mode: 'copy'
	
	input:
	tuple val(sample), path(file1), path(file2), path(bam), path(trans_file), path(logf), path(processed_bam), path(filtered_bam), path(bai)
	path refflat_comprehensive
	path gtf
	path picard

	output:
	path '*/*.flagstat'
	path '*/*.insertSizeMetrics'
	path '*/*.RNA_Metrics'
	path '*/*.pdf'
	path '*/rnaseq_qc_results.txt'
	path '*/*Processed.out.bam.bai'
	
	script:
	"""
	mkdir -p $sample
	samtools index $processed_bam
	mv *bam.bai $sample
	samtools flagstat $processed_bam > ${sample}/${processed_bam}.flagstat
	qualimap rnaseq -pe --java-mem-size=8G -bam $filtered_bam -gtf $gtf -outformat PDF -outdir ${sample} -outfile ${filtered_bam}.qualimap
	java -Xmx4G -jar picard.jar CollectRnaSeqMetrics I=$filtered_bam O=${sample}/${filtered_bam}.comprehensive.RNA_Metrics REF_FLAT=$refflat_comprehensive STRAND=NONE MINIMUM_LENGTH=500 CHART_OUTPUT=${sample}/${filtered_bam}.comprehensive.RNA_Metrics.pdf
	java -Xmx4G -jar picard.jar CollectInsertSizeMetrics I=$filtered_bam O=${sample}/${filtered_bam}.insertSizeMetrics DEVIATIONS=50 INCLUDE_DUPLICATES=true HISTOGRAM_FILE=${sample}/${filtered_bam}.insertSizeMetrics.pdf
	"""
}

process RSEM {
	tag "RSEM transcript quantification..."
	publishDir "${workingDir}/RSEM", mode: 'copy'
	
	input:
	path(bamfile)
	path(rsem_ref)
	
	output:
	path("*rsem*")
	
	script:
	sample = bamfile.baseName - ".Aligned.toTranscriptome.out.bam"
	"""
	rsem-calculate-expression --single-cell-prior --estimate-rspd --strandedness non --alignments --paired-end --num-threads 9 --seed 42 $bamfile $rsem_ref/Gencode_filtercomp.gtf ${sample}.rsem
	"""
}

process BAMTOBED {
	tag "Converting bam to bed and quantifying repeats..."
	publishDir "${workingDir}/BED", mode: 'copy'

	input:
	path(bam_file)
	path(repeat_bed)

	output:
	path "*bed.gz"

	script:
	"""
	sample=\$(basename "$bam_file" | cut -d. -f2)
	bedtools bamtobed -i $bam_file > \${sample}.bed
	bedtools intersect -a \${sample}.bed -b $repeat_bed -wa | sort | uniq | gzip > \${sample}.reads_in_repeats.bed.gz
	gzip \${sample}.bed
	"""
}

process BLAST_UNMAP {
	container null
	tag "NT-blast on unmapped reads..."
	publishDir "${workingDir}/BLAST", mode: 'copy'
	
	input:
	path bamfile

	output:
	path "*unmapped.fastq"
	path "*unmapped.extracted.fa"
	path "*unmapped.blastn_out*"

	script:
	sample = bamfile.baseName - ".Aligned.sortedByCoord.Processed.out.bam"
	"""
	${params.samtools} fastq -f12 ${bamfile} > ${sample}.unmapped.fastq
	${params.seqtk} sample ${sample}.unmapped.fastq -s100 ${params.sample_reads_n} > ${sample}.unmapped.random.fa
	${params.seqtk} seq -A ${sample}.unmapped.random.fa | cut -d ' ' -f 1 > ${sample}.unmapped.extracted.fa && rm ${sample}.unmapped.random.fa
	${params.blastn} -query ${sample}.unmapped.extracted.fa -out ${sample}.unmapped.blastn_out.xml  -evalue 1e-5 -max_target_seqs 1 -outfmt 5 -db ${params.dbnt} -num_threads 32
	${params.python} ${params.parse_blast_xml} ${sample}.unmapped.blastn_out.xml > ${sample}.unmapped.blastn_out_simplified.tsv
	sh ${params.blastn_out_process} ${sample}.unmapped.blastn_out_simplified.tsv ${sample}.unmapped.blastn_out_count.tsv 
	"""
}

process BLAST_INTERGENIC{
	container null
	tag "NT-blast on intergenic reads..."
	publishDir "${workingDir}/BLAST", mode: 'copy'

	input:
	path bamfile

	output:
	path "*intergenic.bam"
	path "*intergenic.extracted.fa"
	path "*intergenic.blastn_out*"

	script:
	"""
	sample=\$(basename "$bam_file" | cut -d. -f2)
	${params.tools_path}/bedtools intersect -abam ${bamfile} -b ${params.intergenic_bedfile} > ${sample}.intergenic.bam
	${params.samtools} fastq ${sample}.intergenic.bam > ${sample}.intergenic.fastq
	${params.seqtk} sample ${sample}.intergenic.fastq -s100 ${params.sample_reads_n} > ${sample}.intergenic.random.fa
	${params.seqtk} seq -A ${sample}.intergenic.random.fa | cut -d ' ' -f 1 > ${sample}.intergenic.extracted.fa && rm ${sample}.intergenic.random.fa
	${params.blastn} -query ${sample}.intergenic.extracted.fa -out ${sample}.intergenic.blastn_out.xml  -evalue 1e-5 -max_target_seqs 1 -outfmt 5 -db ${params.dbnt} -num_threads 32
	${params.python} ${params.parse_blast_xml} ${sample}.intergenic.blastn_out.xml > ${sample}.intergenic.blastn_out_simplified.tsv
	sh ${params.blastn_out_process} ${sample}.intergenic.blastn_out_simplified.tsv ${sample}.intergenic.blastn_out_count.tsv
	"""
}

process SUMMARY_STAT {
	tag "Generating summary statistics file..."
	publishDir "${workingDir}/summary", mode: 'copy'

	input:
	path "*"

	output:
	path "*tsv"

	script:
        """
	python collect.RNAalign.nf.py --workdir ${params.workingDir} --outdir . --samtools ${params.samtools}
	"""
}

// This is the main pipeline
workflow pipeline {
	params.all_reads="${params.workingDir}/rawdata/*R*fastq.gz"
	all_reads = Channel.fromPath(params.all_reads)
	params.reads="${params.workingDir}/rawdata/*R{1,2}*fastq.gz"
	reads = Channel.fromFilePairs(params.reads)

	TRIM(reads)
	trim_rna_ch=TRIM.out.trimmed_reads
	
        if (params.run_fastqc) {
		FASTQC(all_reads)
		html_ch=FASTQC.out.htmls
		zip_ch=FASTQC.out.zips
		FASTQC_TRIM(trim_rna_ch)
		trimmed_html_ch=FASTQC_TRIM.out.htmls
		trimmed_zip_ch=FASTQC_TRIM.out.zips
		MULTIFASTQC(html_ch.mix(zip_ch).collect())
		MULTIFASTQC_TRIM(trimmed_html_ch.mix(trimmed_zip_ch).collect())
	}
	
	STAR_ALIGN(params.star_index, params.sjdbgtf, trim_rna_ch)
	star_align_ch=STAR_ALIGN.out.general_output
	transcriptome_bam_file=STAR_ALIGN.out.transcriptome_bam_file

	STAR_MARKDUP(star_align_ch)
	STAR_MARKDUP_ch=STAR_MARKDUP.out.general_output

	STAR_FILTER(STAR_MARKDUP_ch)
	star_filtered_ch=STAR_FILTER.out.general_output
	star_filtered_ch_featurecount=STAR_FILTER.out.bam

	trim_rna_ch.join(star_align_ch,by:0).set{first_rna_join_ch}
	first_rna_join_ch.join(STAR_MARKDUP_ch,by:0).set{second_rna_join_ch}
	second_rna_join_ch.join(star_filtered_ch,by:0).set{combined_rna_ch}

	BAMTOBED(star_filtered_ch_featurecount, params.repeat_bedfile)
	QUALITY_METRICS(combined_rna_ch, params.refflat_comprehensive, params.final_gtf, params.picard)
	metrics_ch=QUALITY_METRICS.out
	
	if (params.featurecount_gene) {
		quant_ch=STAR_FEATURECOUNTS_GENE(star_filtered_ch_featurecount.collect(), params.sjdbgtf)
	}
		
	if (params.rsem) {
		quant_ch=RSEM(transcriptome_bam_file, params.rsem_ref)
	}
}

# Standalone workflows
workflow raw_qc {
        params.all_reads="${params.workingDir}/rawdata/*R*fastq.gz"
        all_reads = Channel.fromPath(params.all_reads)

        FASTQC(all_reads)
        html_ch=FASTQC.out.htmls
        zip_ch=FASTQC.out.zips
        MULTIFASTQC(html_ch.mix(zip_ch).collect())
}

workflow summary {
        SUMMARY_STAT()
}

workflow featurecounts_gene {
        params.filtered_bam_file = "${params.workingDir}/*STAR/*/filtered*Aligned.sortedByCoord.Processed.out.bam"
        filtered_bam_file = Channel.fromPath(params.filtered_bam_file, checkIfExists: true)
	STAR_FEATURECOUNTS_GENE(filtered_bam_file.collect(), params.sjdbgtf)
}

workflow rsem {
	params.transcriptome_bam_file = "${params.workingDir}/*STAR/*/*Aligned.toTranscriptome.out.bam"
	transcriptome_bam_file = Channel.fromPath(params.transcriptome_bam_file, checkIfExists: true)
	RSEM(transcriptome_bam_file, params.rsem_ref, params.rsem_gtf)
}

# Entry point
workflow {
	pipeline()
}

workflow.onComplete {
        log.info (workflow.success? "Done." : "Something went wrong.")
}

