#!/usr/bin/env nextflow

/*
 * Lab: Lab for Functional Epigenetics, KU Leuven
 * Date: 2024/08
 */

nextflow.enable.dsl=2
workingDir=params.workingDir
params.all_reads="${params.workingDir}/${params.rawdataDir}/*R*fastq.gz"
all_reads = Channel.fromPath(params.all_reads)
params.reads="${params.workingDir}/${params.rawdataDir}/*R{1,2}*fastq.gz"
reads = Channel.fromFilePairs(params.reads)

log.info """\
	scNMT-seq & sciMET gDNA Processing Pipeline V0
	==============================================
	Working dir: ${params.workingDir}
	"""
	.stripIndent()

// Demultiplexing (for sciMET)
process demux_fq {
	tag "Demultiplexing fq files of $sample"
	publishDir "${workingDir}/trim_sc", mode: 'copy'

	input:
	tuple val(sample), file(reads)

	output:
	path("*R*fq.gz"), emit: fq_sc

	script:
	"""
	python /demux_sciMET_fq.py --fq1 ${reads[0]} --fq2 ${reads[1]}
	"""
}

process demux_bam {
	tag "Seperating bam files containing multiple cells into single-cell bam files"
	publishDir "${workingDir}/bam", mode: 'copy'

	input:
	path(bam)

	output:
	path("*sc.bam"), emit: sc_bams

	script:
	"""
	python /demux_sciMET_bam.py -i $bam
	"""
}

process trim_NMT {
	tag "Trimming raw reads"
	publishDir "${workingDir}/trim", mode: 'copy'

	input:
	tuple val(sample), path(reads)
	
	output:
	tuple val(sample), path("*_1.fq.gz"), path("*_2.fq.gz")
	
	script:
	"""
	trim_galore ${params.extra_trim_params_DNA} --length ${params.trim_length_DNA} --clip_R1 ${params.clip_DNA_R1} --clip_R2 ${params.clip_DNA_R2} --three_prime_clip_R1 ${params.three_prime_clip_DNA_R1} --three_prime_clip_R2 ${params.three_prime_clip_DNA_R2} --paired ${reads[0]} ${reads[1]}
	"""
}

process trim_sciMET {
	tag "Trimming raw reads"
	publishDir "${workingDir}/trim", mode: 'copy'

	input:
	tuple val(sample), path(reads)
	
	output:
	tuple(val(sample), path("${sample}_val_{1,2}.fq.gz")), emit: trimmed_reads
	tuple(val(sample), path("${sample}.fq.gz_trimming_report_R*.txt"))

	script:
	"""
	mv ${reads[0]} ${sample}_R1.fq.gz
	mv ${reads[1]} ${sample}_R2.fq.gz
	trim_galore --fastqc -a CTATCTCTTATA -a2 AGATCGGAAGAGC --three_prime_clip_R2 10 --basename ${sample} --paired ${sample}_R1.fq.gz ${sample}_R2.fq.gz
	rm ${sample}_R{1,2}.fq.gz
	mv ${sample}_R2.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R2.txt
	mv ${sample}_R1.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R1.txt
	"""
}

process fastqc {
	tag "Performing fastqc on the raw reads"
	publishDir "${workingDir}/fastqc/fastqc_raw", mode: 'copy'
	
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

process multiqc {
	tag "Performing multiqc"
	publishDir "${workingDir}/fastqc/multiqc_raw", mode: 'copy'

	input:
	path "*"

	output:
	path "*.html"

	script:
	"""
	multiqc . -n multiqc_report_raw_data.html
	"""
}

process fastqc_trim {
	tag "Performing fastqc on trimmed reads"
	publishDir "${workingDir}/fastqc/fastqc_trim", mode: 'copy'

	input:
	tuple val(sample), file(reads)

	output:
	path "*.html", emit: htmls
	path "*.zip", emit: zips

	script:
	"""
	fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
	"""
}

process multiqc_trim {
	tag "Performing multiqc on trimmed reads"
	publishDir "${workingDir}/fastqc/multiqc_trim", mode: 'copy'

	input:
	path "*"

	output:
	path "*.html"

	script:
	"""
	multiqc . -n multiqc_report_trimmed_data.html
	"""
}

// Bismark local alignment with Bowtie2
process align {
	tag "Aligning $sample"
	publishDir "${params.workingDir}/bam/${sample}", mode: 'copy'

	input:
	tuple val(sample), file(reads)
	path(genome_folder)

	output:
	tuple val(sample), path("*pe.bam"), emit: bam
	path("*PE_report.txt")

	script:
	"""
	bismark --local --non_directional --gzip --genome $genome_folder --output_dir . --basename $sample -1 ${reads[0]} -2 ${reads[1]}
	"""
}

// To examine mapping strategy and hybrid fragments in the libraries
process diagnosis {
	tag "Checking $sample"
		publishDir "${params.workingDir}/diagnosis/${sample}", mode: 'copy'

		input:
		tuple val(sample), file(reads)
		path(genome_folder)

		output:
		path("*remap_merged_sorted.bam")
	path("*chimeric_reads.sam")
	path("*_report.txt")
	path("*diagnosis.txt")

	script:
	"""
	# map in PE mode
	bismark --local --non_directional --gzip --genome $genome_folder --output_dir . --basename $sample --unmapped -1 ${reads[0]} -2 ${reads[1]}
	# remap R1 in SE mode
	bismark --local --non_directional --genome $genome_folder --output_dir . --basename ${sample}_unmapped_R1 ${sample}_unmapped_reads_1.fq.gz
	# remap R2 in SE mode
	bismark --local --non_directional --genome $genome_folder --output_dir . --basename ${sample}_unmapped_R2 ${sample}_unmapped_reads_2.fq.gz
	# combine remapped reads
	samtools merge ${sample}_remap_merged.bam ${sample}_unmapped_R1.bam ${sample}_unmapped_R2.bam
	samtools sort -n ${sample}_remap_merged.bam -o ${sample}_remap_merged_sorted.bam && rm ${sample}_remap_merged.bam
	# extract read names
	samtools view ${sample}_remap_merged_sorted.bam | awk '{ split(\$1, a, "_"); print a[1] }' | sort | uniq -c | awk '\$1 == 2 { print \$2 }' > remap.paired_reads.txt
	paired_reads_count=\$(wc -l < remap.paired_reads.txt)
	echo "Number of read pairs after SE alignment for the unmapped R1 and R2: \$paired_reads_count" > ${sample}.diagnosis.txt
	# extract reads
	samtools view ${sample}_remap_merged_sorted.bam | grep -F -f remap.paired_reads.txt > remap.paired_reads.sam
	# identify chimeric reads
	awk '
	{
		if (NR % 2 == 1) {
			read1 = \$0; chr1 = \$3; pos1 = \$4;
		} else {
			read2 = \$0; chr2 = \$3; pos2 = \$4;
			if (chr1 != chr2 || (chr1 == chr2 && sqrt((pos2 - pos1)^2) > 10000)) {
				print read1 "\\n" read2
			}
		}
	}' remap.paired_reads.sam > ${sample}_chimeric_reads.sam
	# count chimeric reads
	chimeric_reads_count=\$(( \$(wc -l < ${sample}_chimeric_reads.sam) / 2 ))
	echo "Number of chimeric read pairs: \$chimeric_reads_count" >> ${sample}.diagnosis.txt
	"""
	}

// To generate a saturation curve for evaluating sequencing depth
process bam_downsample {
	tag "Downsampling bam file"
	publishDir "${workingDir}/bam", mode: 'copy'
	
	input:
	path(bam_file)

	output:
	path("*_downsampled/*_downsampled.bam"), emit: downsampled_bams
	path "${sample}*downsampled"

	script:
	sample = bam_file.baseName - "_pe"

	"""
	parallel -k "
	mkdir -p ${sample}_0{}'_downsampled'
	samtools view -bs 100{} ${bam_file} > '${sample}_0{}'_downsampled/${sample}_0{}'_downsampled.bam'
	" ::: .1 .2 .4 .6 .8
	"""
}

process bam_dedup {
	tag "Deduplicating $sample"
	publishDir "${workingDir}/bam/${sample}", mode: 'copy'

	input:
	tuple val(sample), path(bam_file)

	output:
	tuple val(sample), path("*deduplicated.bam"), emit: dedup_bam
	path "*deduplicated.bam", emit: dedup_bam_file
	path "*.txt", emit: dedup_report

	script:
	"""
	deduplicate_bismark -p -o $sample --output_dir . $bam_file
	"""
}

process dedup_from_bam {
	tag "Deduplicating"
	publishDir "${workingDir}/bam", mode: 'copy'
		
	input:
		path(bam_file)

		output:
		tuple val(sample), path("*/*deduplicated.bam"), emit: dedup_bam
		path "*/*deduplicated.bam", emit: dedup_bam_file

		script:
		sample = bam_file.baseName
		"""
		deduplicate_bismark -p -o $sample --output_dir $sample $bam_file
		"""
}

process bam_filter {
	tag "Filtering reads for $sample"
	publishDir "${workingDir}/bam/${sample}", mode: 'copy'
	
	input:
	tuple val(sample), path(bam_file)
	
	output:
	tuple val(sample), path("*deduplicated.filtered.sorted.bam"), path("*deduplicated.filtered.sorted.bam.bai"), emit: final_bam_bs_input
	tuple path("*deduplicated.filtered.sorted.bam"), path("*deduplicated.filtered.sorted.bam.bai"), emit: final_bam
	path "*deduplicated.filtered.sorted.bam", emit: final_bam_file
	
	script:
	"""
	samtools view -bS -@ ${task.cpus} ${params.samtools_filter_bs} -o ${sample}.merge.deduplicated.filtered.bam $bam_file
	samtools sort -O bam -o ${sample}.merge.deduplicated.filtered.sorted.bam ${sample}.merge.deduplicated.filtered.bam
	samtools index ${sample}.merge.deduplicated.filtered.sorted.bam
	"""
}

// sciMET: demultiplex, dedup, filter and compute metrics
process bam_demux_process {
	tag "Processing single-cell bam files of $sample"
	publishDir "${workingDir}/bam_sc", mode: 'copy'

	input:
	tuple val(sample), path(bam)

	output:
	path "*/*sc.bam", emit: sc_bam, optional: true
	path "*/*sc.deduplicated.bam", emit: dedup_bam_file, optional: true
	tuple path("*/*sc.deduplicated.filtered.sorted.bam"), path("*/*sc.deduplicated.filtered.sorted.bam.bai"), emit: final_bam, optional: true
	path "*/*sc.deduplicated.filtered.sorted.bam", emit: final_bam_file, optional: true
	path "*/*deduplication_report.txt", optional: true
	path "*mapping_summary.tsv", emit: mapping_metrics

	script:
	"""
	python /demux_sciMET_bam.py -i $bam
	
	for file in *.sc.bam; do
		prefix=\$(echo "\$file" | awk -F'.sc.bam' '{print \$1}')
	
		# dedup
		deduplicate_bismark -p --bam \$file
	
		# filter (unsorted bam is needed for bismark_methylation_extractor)
		samtools view -bS -@ ${task.cpus} ${params.samtools_filter_bs} -o \${prefix}.sc.deduplicated.filtered.bam \${prefix}.sc.deduplicated.bam

		# sort
		samtools sort -O bam -o \${prefix}.sc.deduplicated.filtered.sorted.bam \${prefix}.sc.deduplicated.filtered.bam && rm \${prefix}.sc.deduplicated.filtered.bam
		# index
		samtools index \${prefix}.sc.deduplicated.filtered.sorted.bam

		# flagstat
		samtools flagstat \$file > \${prefix}.mapped.flagstat
		samtools flagstat \${prefix}.sc.deduplicated.bam > \${prefix}.dedup.flagstat
		samtools flagstat \${prefix}.sc.deduplicated.filtered.sorted.bam > \${prefix}.filter.flagstat
		python /summarize_bam_flagstat.py --flagstat1 \${prefix}.mapped.flagstat --flagstat2 \${prefix}.dedup.flagstat --flagstat3 \${prefix}.filter.flagstat --output_prefix \${prefix} && rm *flagstat
		
		# filter single-cell output files (remove redundancy) based on the number of uniquely mapped reads
		filtered_unique=\$(awk -F'\\t' '{print \$4}' "\${prefix}.mapping_summary.tsv")
		if [[ \$filtered_unique -lt ${params.uniq_cutoff} ]]; then
			rm -f \${prefix}*.bam* \${prefix}*deduplication_report.txt
		else
			mkdir -p \${prefix}
			mv \${prefix}*sc*bam* \${prefix}*deduplication_report.txt \${prefix}
		fi
	done
	"""
}

// Compute genome-wide coverage
process coverage_unfiltered {
	tag "Computing coverage"
	publishDir "${workingDir}/cov", mode: 'copy'
	
	input:
	path bam_file
	
	output:
	path "*_unfiltered.Coverage.tsv"
	path "*_unfiltered.Coverage.tsv.genome.wide", emit: unfiltered_cov
	
	script:
	"""
	sample=\$(basename "$bam_file" | awk -F '.deduplicated.bam' '{print \$1}' | awk -F '.merge|.sc' '{print \$1}')
	samtools sort -O bam -T \${sample}.sort -o \${sample}.deduplicated.sorted.bam $bam_file
	samtools coverage \${sample}.deduplicated.sorted.bam | head -26 > \${sample}_unfiltered.Coverage.tsv
	sh /getGenomeWideCoverage.sh \${sample}_unfiltered.Coverage.tsv \$sample && rm \${sample}.deduplicated.sorted.bam
	"""
}

process coverage_filtered {
	tag "Computing coverage"
	publishDir "${workingDir}/cov", mode: 'copy'
	
	input:
	path bam_file
	
	output:
	path "*_filtered.Coverage.tsv"
	path "*_filtered.Coverage.tsv.genome.wide", emit: filtered_cov
	
	script:
	"""
	sample=\$(basename "$bam_file" | sed 's/\\.merged\\.deduplicated\\.filtered\\.sorted\\.bam.*//')
	samtools coverage $bam_file | head -26 > \${sample}_filtered.Coverage.tsv
	sh /getGenomeWideCoverage.sh \${sample}_filtered.Coverage.tsv \$sample
	"""
}

process coverage_merge {
		tag "Merging coverage results"
		publishDir "${workingDir}/summary", mode: 'copy'

		input:
		path(unfiltered_cov)
		path(filtered_cov)

		output:
		path("combined.coverage.genome.wide.tsv")

		script:
		"""
		cat *_unfiltered.Coverage.tsv.genome.wide | awk 'BEGIN{FS=OFS="\\t"} {split(\$1,a,"\\.merge|\\.sc"); \$1=a[1]; print}' | sort -t \$'\\t' -k 1,1 > unfiltered.coverage.genome.wide.tsv
		cat *_filtered.Coverage.tsv.genome.wide | awk 'BEGIN{FS=OFS="\\t"} {split(\$1,a,"\\.merge|\\.sc"); \$1=a[1]; print}' | sort -t \$'\\t' -k 1,1 > filtered.coverage.genome.wide.tsv
		join -t \$'\\t' -1 1 -2 1 unfiltered.coverage.genome.wide.tsv filtered.coverage.genome.wide.tsv > combined.coverage.genome.wide.tsv
		"""
}

process coverage_merge_local {
	container null
	tag "Merging coverage results"
	publishDir "${workingDir}/summary", mode: 'copy'
	
	input:
	path(unfiltered_cov)
	path(filtered_cov)
	
	output:
	path 'combined.coverage.genome.wide.downsampled.tsv'

	script:
	"""
	cat ${workingDir}/cov/*_unfiltered.Coverage.tsv.genome.wide | awk 'BEGIN{FS=OFS="\\t"} {split(\$1,a,"\\.merge|\\.sc"); \$1=a[1]; print}' | sort -t \$'\\t' -k 1,1 > unfiltered.coverage.genome.wide.tsv
	cat ${workingDir}/cov/*_filtered.Coverage.tsv.genome.wide | awk 'BEGIN{FS=OFS="\\t"} {split(\$1,a,"\\.merge|\\.sc"); \$1=a[1]; print}' | sort -t \$'\\t' -k 1,1 > filtered.coverage.genome.wide.tsv
	join -t \$'\\t' -1 1 -2 1 unfiltered.coverage.genome.wide.tsv filtered.coverage.genome.wide.tsv > combined.coverage.genome.wide.downsampled.tsv 
	"""
}

process bam_to_bed {
	tag "Converting bam to bed"
	publishDir "${workingDir}/bedfile", mode: 'copy'

	input:
	path(bam)

	output:
	path "*bed.gz"

	script:
	"""
	sample=\$(echo "$bam" | awk -F '.bam' '{print \$1}' | awk -F '.merge|.sc' '{print \$1}')
	bedtools bamtobed -i $bam | gzip > \${sample}.bed.gz
	"""
}

process insert_size {
	tag "Collecting insert size"
	publishDir "${workingDir}/insert", mode: 'copy'

	input:
	path(bam_file)
	path(picard)

	output:
	path "*metrics.txt"
	path "*metrics.pdf"

	script:
	"""
	prefix=\$(echo "$bam_file" | awk -F'.bam' '{print \$1}')
	java -Xmx1G -jar $picard CollectInsertSizeMetrics DEVIATIONS=50 INCLUDE_DUPLICATES=false O=\${prefix}.insert_size_metrics.txt H=\${prefix}.insert_size_metrics.pdf I=$bam_file
	"""
}

// Extract methylation status by Bismark
process bs_extract_NMT {
	tag "Extracting methylation using bismark_methylation_extractor"
	publishDir "${workingDir}/bismark_extract", mode: 'copy'

	input:
	tuple path(bam), path(bai)
	path(genome_folder)

	output:
	path "*splitting_report.txt"
	path "*M-bias.txt"
	path "*bedGraph.gz"
	path "*bismark.cov.gz"
	path "*cytosine_context_summary.txt"

	script:
	"""
	samtools view -h -b $bam \$(echo chr{1..22} chrX chrY) > tmp.bam && rm $bam
	samtools sort -n -o $bam tmp.bam && rm tmp.bam
	prefix=\$(echo "$bam" | awk -F'.bam' '{print \$1}')
	bismark_methylation_extractor --CX -p --no_overlap --bedGraph --gzip --parallel ${task.cpus} $bam
	coverage2cytosine --nome-seq --gzip --genome_folder $genome_folder -o \${prefix} \${prefix}.bismark.cov.gz	
	"""
}

process bs_extract_sciMET {
	tag "Extracting mC by Bismark"
	publishDir "${workingDir}/bismark_extract"

	input:
	tuple path(bam)
	path(genome_folder)

	output:
	path "*splitting_report.txt"
	path "*M-bias.txt"
	path "*bedGraph.gz"
	path "*.cov.gz"
	//path "*cytosine_context_summary.txt" // coverage2cytosine output
	//path "*CX_report.txt.gz" // coverage2cytosine output

	script:
	"""
	samtools index $bam && samtools view -h -b $bam \$(echo chr{1..22} chrX chrY) > tmp.bam && rm $bam
	samtools sort -n -o $bam tmp.bam && rm tmp.bam
	prefix=\$(echo "$bam" | awk -F'.bam' '{print \$1}')
	bismark_methylation_extractor --CX --comprehensive -p --no_overlap --bedGraph --gzip --parallel ${task.cpus} --genome_folder $genome_folder $bam
	#coverage2cytosine --CX_context --gzip --genome_folder $genome_folder -o \${prefix} \${prefix}.bismark.cov.gz
	"""
}

// Extract ALLC
process extract_allc {
	container params.allc_container
	tag "Extracting overall mC info using methylpy"
	publishDir "${workingDir}/allc", mode: 'copy'
	
	input:
	path ref_fasta
	path ref_fai
	path(bam)
	
	output:
	path "*allc.tsv.gz", emit: allc

	script:
	"""
	samtools index $bam && samtools view -h -b $bam \$(echo chr{1..22} chrX chrY) > input.bam
	sample=\$(echo "$bam" | awk -F '.bam' '{print \$1}' | awk -F '.merge|.sc' '{print \$1}')
	methylpy call-methylation-state --num-upstream-bases 1 --input-file input.bam --paired-end True --sample \$sample --ref-fasta $ref_fasta --num-procs ${task.cpus}
	mv allc_\${sample}.tsv.gz \${sample}.allc.tsv.gz	
	"""
}

// Extract cytosine context-specific ALLC
process extract_allc_contexts_NMT {
	container params.allc_container
	tag "Extracting CpG and GpC from overall allc files"
	publishDir "${workingDir}/allc_context", mode: 'copy'
	
	input:
	path allc_file
	path chromsize_file
	
	output:
	path("*WCG*allc.tsv.gz"), emit: wcg
	path("*GCH*allc.tsv.gz"), emit: gch
	path("*HCH*allc.tsv.gz"), emit: hch
	tuple path("*WCG*.prof"),  path("*WCG*.tbi"), emit: input_mcds_wcg
	tuple path("*GCH*.prof"),  path("*GCH*.tbi"), emit: input_mcds_gch
	tuple path("*HCH*.prof"),  path("*HCH*.tbi"), emit: input_mcds_hch
	path "*Both.allc.tsv.gz.tbi", emit: allc_tbi

	script:
	"""
	sample=\$(basename "$allc_file" | awk -F'.allc.tsv.gz' '{print \$1}')

	allcools standardize-allc --allc_path $allc_file --chrom_size_path $chromsize_file
	
	# NMT cytosine contexts
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts WCGN
	allcools profile-allc --allc_path \${sample}.WCGN-Both.allc.tsv.gz --output_path \${sample}.WCGN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.WCGN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts GCHN
	allcools profile-allc --allc_path \${sample}.GCHN-Both.allc.tsv.gz --output_path \${sample}.GCHN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.GCHN-Both.allc.tsv.gz --chrom_size_path $chromsize_file

	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts HCHN
	allcools profile-allc --allc_path \${sample}.HCHN-Both.allc.tsv.gz --output_path \${sample}.HCHN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.HCHN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	"""
}

process extract_allc_contexts_sciMET {
	container params.allc_container
	tag "Extract mC contexts from allc files"
	publishDir "${workingDir}/allc_context", mode: 'copy'

	input:
	path allc_file
	path chromsize_file

	output:
	path("*allc.tsv.gz")
	path("*.prof")
	path("*.tbi")

	script:
	"""
	sample=\$(basename "$allc_file" | awk -F'.allc.tsv.gz' '{print \$1}')

	allcools standardize-allc --allc_path $allc_file --chrom_size_path $chromsize_file

	# canonical cytosine contexts
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts NCGN
	allcools profile-allc --allc_path \${sample}.NCGN-Both.allc.tsv.gz --output_path \${sample}.NCGN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.NCGN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts NCHN
	allcools profile-allc --allc_path \${sample}.NCHN-Both.allc.tsv.gz --output_path \${sample}.NCHN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.NCHN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	
	# NMT cytosine contexts
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts WCGN
	allcools profile-allc --allc_path \${sample}.WCGN-Both.allc.tsv.gz --output_path \${sample}.WCGN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.WCGN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts GCHN
	allcools profile-allc --allc_path \${sample}.GCHN-Both.allc.tsv.gz --output_path \${sample}.GCHN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.GCHN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	allcools extract-allc --allc_path $allc_file --output_prefix \${sample}. --chrom_size_path $chromsize_file --mc_contexts HCHN
	allcools profile-allc --allc_path \${sample}.HCHN-Both.allc.tsv.gz --output_path \${sample}.HCHN-Both.allc.tsv.gz.prof
	allcools standardize-allc --allc_path \${sample}.HCHN-Both.allc.tsv.gz --chrom_size_path $chromsize_file
	"""
}


// Main workflow for scNMT-seq data processing
workflow pipeline {
	trim_output_ch = trim_NMT(reads)
	trim_NMT_ch = trim_output_ch.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
	
	if (params.run_fastqc) {
		fastqc(all_reads)
		html_ch=fastqc.out.htmls
		zip_ch=fastqc.out.zips
		fastqc_trim(trim_NMT_ch)
		trimmed_html_ch=fastqc_trim.out.htmls
		trimmed_zip_ch=fastqc_trim.out.zips
		multiqc(html_ch.mix(zip_ch).collect())
		multiqc_trim(trimmed_html_ch.mix(trimmed_zip_ch).collect())
	}
	
	align(trim_NMT_ch, params.genome_folder)
	align_ch=align.out.bam
	
	bam_dedup(align_ch)
	dedup_bam_ch = bam_dedup.out.dedup_bam
	dedup_bam_file = bam_dedup.out.dedup_bam_file	
	
	bam_filter(dedup_bam_ch)
	final_bam = bam_filter.out.final_bam
	final_bam_file = bam_filter.out.final_bam_file
	
	if (params.bismark_extract) {
		bs_extract_NMT(final_bam, params.genome_folder)
	}

	if (params.generate_bedfile) {
		bam_to_bed(final_bam_file)
	}

	insert_size(final_bam_file, params.picard)

	coverage_unfiltered(dedup_bam_file)
	coverage_filtered(final_bam_file)
	unfiltered_cov = coverage_unfiltered.out.unfiltered_cov
	filtered_cov = coverage_filtered.out.filtered_cov
	coverage_merge(unfiltered_cov.collect(), filtered_cov.collect())

	if (params.allc_extract) {
		extract_allc(params.ref_fasta, params.ref_fai, final_bam)
		allc=extract_allc.out.allc
		extract_allc_contexts_NMT(allc, params.chromsizes)
	}
}

// Main workflow for sciMET data processing
workflow pipeline_sciMET {
	trim_sciMET(reads)
	trim_dna_ch = trim_sciMET.out.trimmed_reads
	//trim_dna_ch = Channel.fromFilePairs("${params.workingDir}/trim/*/*_val_{1,2}*fq.gz", checkIfExists: true)

		if (params.run_fastqc) {
				fastqc(all_reads)
				html_ch=fastqc.out.htmls
				zip_ch=fastqc.out.zips
				fastqc_trim(trim_NMT_ch)
				trimmed_html_ch=fastqc_trim.out.htmls
				trimmed_zip_ch=fastqc_trim.out.zips
				multiqc(html_ch.mix(zip_ch).collect())
				multiqc_trim(trimmed_html_ch.mix(trimmed_zip_ch).collect())
		}
	
	align(trim_dna_ch, params.genome_folder)
	bs_align_ch=align.out.bam

	bam_demux_process(bs_align_ch)
	dedup_bam_file = bam_demux_process.out.dedup_bam_file.flatten()
	final_bam_file = bam_demux_process.out.final_bam_file.flatten()
	final_bam = bam_demux_process.out.final_bam.map { tuple -> return tuple }

	if (params.bismark_extract) {	
		bs_extract_sciMET(final_bam_file, params.genome_folder)
	}

	if (params.generate_bedfile) {
		bam_to_bed(final_bam_file)
	}

		insert_size(final_bam_file, params.picard)

		coverage_unfiltered(dedup_bam_file)
		coverage_filtered(final_bam_file)
		unfiltered_cov = coverage_unfiltered.out.unfiltered_cov
		filtered_cov = coverage_filtered.out.filtered_cov
		coverage_merge(unfiltered_cov.collect(), filtered_cov.collect())

		if (params.allc_extract) {	
		extract_allc(params.ref_fasta, params.ref_fai, final_bam_file)
		allc=extract_allc.out.allc
		extract_allc_contexts_sciMET(allc, params.chromsizes)
	}
}

// Standalone workflows (feel free to customize)
workflow raw_qc {
	fastqc(all_reads)
	html_ch=fastqc.out.htmls
	zip_ch=fastqc.out.zips
	multiqc(html_ch.mix(zip_ch).collect())
}

workflow diagnosis_chimeric {
		trim_output_ch = trim_NMT(reads)
		trim_NMT_ch = trim_output_ch.map { sample, read1, read2 -> tuple(sample, [read1, read2]) }
	diagnosis(trim_NMT_ch, params.genome_folder)
}

workflow rerun_bam_process_sciMET { // in the case that a single-cell bam file is found to be truncated
	bs_align_ch = Channel.fromPath("${params.workingDir}/bam_rerun/*.bam", checkIfExists: true)
				 .map { it -> tuple(it.baseName, it) }
	bam_demux_process(bs_align_ch)
	dedup_bam_file = bam_demux_process.out.dedup_bam_file.flatten()
}

workflow downsample {
	bam_file_path = "${params.workingDir}/bam/*/*_pe.bam"
	bam_file = Channel.fromPath(bam_file_path, checkIfExists: true)	
	bam_downsample(bam_file)
	downsampled_bam_ch = bam_downsample.out.downsampled_bams.flatMap { it }.map { bam_file -> bam_file }

	dedup_from_bam(downsampled_bam_ch)
	dedup_bam_ch = dedup_from_bam.out.dedup_bam
	dedup_bam_file = dedup_from_bam.out.dedup_bam_file

	bam_filter(dedup_bam_ch)
	dedup_filtered_sorted_bam_ch = bam_filter.out.final_bam_file

	coverage_unfiltered(dedup_bam_file)
	coverage_filtered(dedup_filtered_sorted_bam_ch)
	unfiltered_cov = coverage_unfiltered.out.unfiltered_cov
	filtered_cov = coverage_filtered.out.filtered_cov
	coverage_merge_local(unfiltered_cov.collect(), filtered_cov.collect())
}

workflow coverage {
	dedup_bam_ch = Channel.fromPath("${params.workingDir}/bam/*/*deduplicated.bam", checkIfExists: true)
	dedup_filtered_sorted_bam_ch = Channel.fromPath("${params.workingDir}/*Bismark/*/*deduplicated.filtered.sorted.bam", checkIfExists: true)
	coverage_unfiltered(dedup_bam_ch)
	coverage_filtered(dedup_filtered_sorted_bam_ch)
	unfiltered_cov = coverage_unfiltered.out.unfiltered_cov
	filtered_cov = coverage_filtered.out.filtered_cov
	coverage_merge(unfiltered_cov.collect(), filtered_cov.collect())
}

workflow convert_bam_to_bed {
	params.dedup_filtered_bam_file = "${params.workingDir}/bam/*/*deduplicated.filtered.sorted.bam"
	dedup_filtered_bam_file = Channel.fromPath(params.dedup_filtered_bam_file, checkIfExists: true)
	bam_to_bed(dedup_filtered_bam_file)
}

workflow demux_fq_file {
	trim_dna_ch = Channel.fromFilePairs("${params.workingDir}/trim/*/*_val_{1,2}*fq.gz", checkIfExists: true)
	DEMUX_FQ(trim_dna_ch)
}

workflow demux_bam_file {
	bam = Channel.fromPath("${params.workingDir}/bamDeDup/*/*dedup.nsrt.bam", checkIfExists: true)
	DEMUX_SC(bam)
}

workflow analysis_NMT {
	final_bam = Channel.fromFilePairs("${params.workingDir}/bam/*/*deduplicated.filtered.sorted.bam", 
										size: 2, 
										checkIfExists: true) { file -> 
											file.name.replace('.bam', '')
										}

	extract_allc(params.ref_fasta, extract_allc(params.ref_fasta, final_bam))
	allc=extract_allc.out.allc
	extract_allc_contexts_NMT(allc, params.chromsizes)
}

// Entry point
workflow {
	pipeline()
}

workflow.onComplete {
	println ( workflow.success ? """
	Pipeline execution summary
	---------------------------
	Completed at: ${workflow.complete}
	Duration	: ${workflow.duration}
	Success	 : ${workflow.success}
	Directory   : ${workflow.workingDir}
	Exit status : ${workflow.exitStatus}
	""" : """
	Failed: ${workflow.errorReport}
	Exit status : ${workflow.exitStatus}
	"""
	)
}
