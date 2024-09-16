#!/usr/bin/env python
import os
import re
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor

def process_sample(sampleid, workdir, qc_dir,map_dir,bed_dir,count_dir, samtools):

    rawreadpairs = totalfrags = totalreads = avglen = uniqmap = uniqmaprate = reads_in_filtered = readpairs_in_filtered = avgmaplen = duprate = ercccount = genecount = mediancv = fiveprime = threeprime = five2three = annosplice = gtsplice = gcsplice = atsplice = ncansplice = mismatch = multirate = manyrate = mismrate = shortrate = otherrate = chimericrate = exonrate = intronrate = intergenicrate = overlappingexonrate = exoncount = introncount = intergeniccount = overlappingexoncount = spliced_reads = reads_in_repeats = -1
    
    procline0 = '******\nProcessing sample: %s\n******' % (sampleid)
    print(procline0)
    
    zip_file = os.path.join(qc_dir, sampleid+'_R1_001_fastqc.zip')
    procline00 = 'Processing fastqc output: %s' % (zip_file)
    print(procline00)
    if os.path.isfile(zip_file):
        command = f"unzip -p {zip_file} */fastqc_data.txt | grep 'Total Sequences' | cut -f 2"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            rawreadpairs = int(stdout.decode().strip())
        else:
            print(f"Error processing {zip_file}: {stderr.decode()}")
    else:
        print("File does not exist")

    sample_sub_workdir = os.path.join(map_dir, sampleid) # One folder for each sample in the STAR alignment path
    
    filtered_bam = '%s/filtered.%s.Aligned.sortedByCoord.Processed.out.bam' % (sample_sub_workdir, sampleid)
    procline1 = 'Processing filtered bam file: %s' % (filtered_bam)
    print(procline1)
    if os.path.isfile(samtools):
        command_splice = f"{samtools} view {filtered_bam} | awk '$6 ~ /N/ {{print $1}}' | sort | uniq | wc -l"
        process = subprocess.run(command_splice, shell=True, stdout=subprocess.PIPE)
        spliced_reads = process.stdout.decode().strip()
    else:
        print("Path for samtools does not exist - proceed without counting spliced reads")
    
    rep_bed = os.path.join(bed_dir, sampleid+".reads_in_repeats.bed.gz")
    procline1 = 'Processing filtered bed file: %s' % (rep_bed)
    print(procline1)
    if os.path.isfile(rep_bed):
        command_rep = f"zcat {rep_bed} | wc -l"
        process = subprocess.run(command_rep, shell=True, stdout=subprocess.PIPE)
        reads_in_repeats = process.stdout.decode().strip()
    else:
        print("File does not exist")
    
    samplelogfh0 = '%s/%s.Log.final.out' % (sample_sub_workdir,sampleid)
    procline1 = 'Processing STAR summary log file: %s' % (samplelogfh0)
    print(procline1)
    if os.path.isfile(samplelogfh0):
        with open(samplelogfh0, 'r') as fh0:
            for l0 in fh0:
                if 'Number of input reads' in l0:
                    totalfrags = l0.strip().split("|")[1]
                    totalfrags = re.sub(r"\s+", "", totalfrags, flags=re.UNICODE)
                elif 'Average input read length' in l0:
                    avglen = l0.strip().split("|")[1]
                    avglen = re.sub(r"\s+", "", avglen, flags=re.UNICODE)
                elif 'Uniquely mapped reads number' in l0:
                    uniqmap = l0.strip().split("|")[1]
                    uniqmap = re.sub(r"\s+", "", uniqmap, flags=re.UNICODE)
                elif 'Uniquely mapped reads %' in l0:
                    uniqmaprate = l0.strip().split("|")[1]
                    uniqmaprate = re.sub(r"\s+","", uniqmaprate, flags=re.UNICODE)
                elif 'Average mapped length' in l0:
                    avgmaplen = l0.strip().split("|")[1]
                    avgmaplen = re.sub(r"\s+", "", avgmaplen, flags=re.UNICODE)
                elif 'sjdb' in l0:
                    annosplice = l0.strip().split("|")[1]
                    annosplice = re.sub(r"\s+", "", annosplice, flags=re.UNICODE)
                elif 'GT/AG' in l0:
                    gtsplice = l0.strip().split("|")[1]
                    gtsplice = re.sub(r"\s+", "", gtsplice, flags=re.UNICODE)
                elif 'GC/AG' in l0:
                    gcsplice = l0.strip().split("|")[1]
                    gcsplice = re.sub(r"\s+", "", gcsplice, flags=re.UNICODE)
                elif 'AT/AC' in l0:
                    atsplice = l0.strip().split("|")[1]
                    atsplice = re.sub(r"\s+", "", atsplice, flags=re.UNICODE)
                elif 'Non-canonical' in l0:
                    ncansplice = l0.strip().split("|")[1]
                    ncansplice = re.sub(r"\s+", "", ncansplice, flags=re.UNICODE)
                elif 'Mismatch rate per base' in l0:
                    mismatch = l0.strip().split("|")[1]
                    mismatch = re.sub(r"\s+", "", mismatch, flags=re.UNICODE)
                elif '% of reads mapped to multiple loci' in l0:
                    multirate = l0.strip().split("|")[1]
                    multirate = re.sub(r"\s+", "", multirate, flags=re.UNICODE)
                elif '% of reads mapped to too many loci' in l0:
                    manyrate = l0.strip().split("|")[1]
                    manyrate = re.sub(r"\s+", "", manyrate, flags=re.UNICODE)
                elif '% of reads unmapped: too many mismatches' in l0:
                    mismrate = l0.strip().split("|")[1]
                    mismrate = re.sub(r"\s+","", mismrate, flags=re.UNICODE)
                elif '% of reads unmapped: too short' in l0:
                    shortrate = l0.strip().split("|")[1]
                    shortrate = re.sub(r"\s+", "", shortrate, flags=re.UNICODE)
                elif '% of reads unmapped: other' in l0:
                    otherrate = l0.strip().split("|")[1]
                    otherrate = re.sub(r"\s+", "", otherrate, flags=re.UNICODE)
                elif '% of chimeric reads' in l0:
                    chimericrate = l0.strip().split("|")[1]
                    chimericrate = re.sub(r"\s+","",chimericrate, flags=re.UNICODE)
    else:
        print("File does not exist")

    samplelogfh1 = '%s/%s.Aligned.sortedByCoord.Processed.out.bam.flagstat' % (sample_sub_workdir,sampleid)
    procline1 = 'Processing bam flagstat file: %s' % (samplelogfh1)
    print(procline1)
    if os.path.isfile(samplelogfh1):
        with open(samplelogfh1, 'r') as fh1:
            for l1 in fh1:
                if 'in total' in l1:
                    totalreads = int(l1.strip().split(' ')[0])
                elif 'duplicates' in l1:
                    dupreads = int(l1.strip().split(' ')[0])

            duprate = dupreads/totalreads
            duprate = '{:.2%}'.format(duprate)
    else:
        print("File does not exist")

    samplelogfh3 = '%s/filtered.%s.Aligned.sortedByCoord.Processed.out.bam.comprehensive.RNA_Metrics' % (sample_sub_workdir, sampleid)
    procline3 = 'Processing picard RNA Metrics: %s\n' % (samplelogfh3)
    print(procline3)
    if os.path.isfile(samplelogfh3):
        with open(samplelogfh3, 'r') as fh3:
            for l3 in fh3:
                if l3.startswith('PF_BASES'):
                    indexline = l3.strip().split('\t')
                    mediancv_pos = int(indexline.index('MEDIAN_CV_COVERAGE'))
                    fiveprime_pos = int(indexline.index('MEDIAN_5PRIME_BIAS'))
                    threeprime_pos = int(indexline.index('MEDIAN_3PRIME_BIAS'))
                    five2three_pos = int(indexline.index('MEDIAN_5PRIME_TO_3PRIME_BIAS'))
                    valueline = next(fh3, '').split('\t')
                    mediancv = valueline[mediancv_pos]
                    fiveprime = valueline[fiveprime_pos]
                    threeprime = valueline[threeprime_pos]
                    five2three = valueline[five2three_pos]
    else:
        print("File not exist")

    samplelogfh4 = '%s/featureCount.tsv' % (count_dir)
    procline4 = 'Processing featureCounts: %s\n' % (samplelogfh4)
    print(procline4)
    genecount = ercccount = 0
    if os.path.isfile(samplelogfh4):
        with open(samplelogfh4, 'r') as fh4:
            for l4 in fh4:
                if l4.startswith('Geneid'):
                    bams = l4.strip().split('\t')[6:]
                    samples = {b.split('.Aligned')[0].split('filtered.')[1] : (i+6) for i, b in enumerate(bams)}
                    idx = samples[sampleid]
                elif l4.startswith("ENSG"):
                    transcripthit = int(l4.strip().split('\t')[idx])
                    if (transcripthit > 0):
                        genecount += 1
                elif any(l4.startswith(prefix) for prefix in ["ERCC", "NIST"]):
                    ercc_count = int(l4.strip().split('\t')[idx])
                    ercccount += ercc_count
    else:
        print("File does not exist")

    samplelogfh5 = '%s/rnaseq_qc_results.txt' % (sample_sub_workdir)
    procline5 = 'Processing RNA qualimap: %s\n' % (samplelogfh5)
    print(procline5)
    if os.path.isfile(samplelogfh5):
        with open(samplelogfh5, 'r') as fh5:
            for l5 in fh5:
                cl5 = l5.lstrip(' ')
                if cl5.startswith('total alignments'):
                    reads_in_filtered0 = cl5.strip().split('=')[1].strip()
                    reads_in_filtered = reads_in_filtered0.replace(",", "")
                elif cl5.startswith('read pairs aligned'):
                    readpairs_in_filtered0 = cl5.strip().split('=')[1].strip()
                    readpairs_in_filtered = readpairs_in_filtered0.replace(",", "")
                elif cl5.startswith('exonic'):
                    exoncount = cl5.strip().split('(')[0].strip().split('=')[1].strip().replace(",", "")
                    exonrate = cl5.strip().split('(')[1].split('%')[0]
                elif cl5.startswith('intronic'):
                    introncount = cl5.strip().split('(')[0].strip().split('=')[1].strip().replace(",", "")
                    intronrate = cl5.strip().split('(')[1].split('%')[0]
                elif cl5.startswith('intergenic'):
                    intergeniccount = cl5.strip().split('(')[0].strip().split('=')[1].strip().replace(",", "")
                    intergenicrate = cl5.strip().split('(')[1].split('%')[0]
                elif cl5.startswith('overlapping exon'):
                    overlappingexoncount = cl5.strip().split('(')[0].strip().split('=')[1].strip().replace(",", "")
                    overlappingexonrate = cl5.strip().split('(')[1].split('%')[0]

    metrics_ori = (sampleid, rawreadpairs, totalfrags, totalreads, avglen, uniqmap, uniqmaprate, reads_in_filtered, readpairs_in_filtered, avgmaplen, duprate, ercccount, genecount, mediancv, fiveprime, threeprime, five2three,annosplice, gtsplice, gcsplice, atsplice, ncansplice, mismatch, multirate, manyrate, mismrate, shortrate, otherrate, chimericrate, exonrate, intronrate, intergenicrate, overlappingexonrate, exoncount, introncount, intergeniccount, overlappingexoncount, spliced_reads, reads_in_repeats)
    metrics = tuple(
            element.replace('%', '') if isinstance(element, str) else element
            for element in metrics_ori
            )
    pfline = ('%s\t' * (len(metrics)-1) + '%s\n') % metrics
    return pfline

def run_stats(workdir, outdir, qcdir, mapdir, beddir, countsdir, samtools):
    workdir = os.path.abspath(workdir)
    print('Working dir: '+workdir)
    sampleids_list = os.listdir(os.path.join(workdir, mapdir))
    sampleids = [s for s in sampleids_list if 'downsampled' not in s]
    #outdir = os.path.join(workdir, outdir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    qc_dir = os.path.join(workdir, qcdir)
    map_dir = os.path.join(workdir, mapdir)
    bed_dir = os.path.join(workdir, beddir)
    count_dir = os.path.join(workdir, countsdir)
    outfh = os.path.join(outdir, 'RNAalign.batch.stats.tsv')

    with open(outfh, 'w') as fo:
        fo.write('Sample\tRaw_readpairs\tTrimmed_readpairs\tTrimmed_reads\tAvg_length\tUniq_map_readpairs\tUniq_map_rate\tUniq_map_reads_filtered\tUniq_map_readpairs_filtered\tAvg_map_length\tDuplication_rate\tUniq_map_ERCC\tGene_count_exonic\tMEDIAN_CV_COVERAGE\tMEDIAN_5PRIME_BIAS\tMEDIAN_3PRIME_BIAS\tMEDIAN_5PRIME_TO_3PRIME_BIAS\tAnnotate_splices\tGT/AG_splices\tGC/AG_splices\tAT/AC_splices\tNon-canonical_splices\tMismatch_rate\tMulti_map_rate\tToomanyloci_rate\tUnmap_mismatch_rate\tUnmap_tooshort_rate\tUnmap_other_rate\tChimeric_rate\tExonic_rate\tIntronic_rate\tIntergenic_rate\tOverlapping_exon_rate\tExonic_count\tIntronic_count\tIntergenic_count\tOverlapping_exon_count\tSpliced_reads\tReads_in_repeats\n')
        
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_sample, sampleid, workdir, qc_dir, map_dir, bed_dir, count_dir, samtools) for sampleid in sampleids]
            for future in futures:
                try:
                    pfline = future.result()
                    fo.write(pfline)
                except Exception as e:
                    print(f"An error occurred: {e}")

    print(f"\nSummary file generated: {outfh}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, default='.', help='Path to the root working directory')
    parser.add_argument('--qcdir', type=str, default='FASTQC', help='Name of the fastqc directory')
    parser.add_argument('--mapdir', type=str, default='STAR', help='Name of the mapping directory')
    parser.add_argument('--beddir', type=str, default='BED', help='Name of the bed files directory')
    parser.add_argument('--countsdir', type=str, default='FEATURECOUNTS', help='Name of the featureCounts directory')
    parser.add_argument('--outdir', type=str, default='summary', help='Name of the output directory')
    parser.add_argument('--samtools', type=str, default='samtools', help='Path to samtools executable')
    args = parser.parse_args()
    run_stats(args.workdir, args.outdir, args.qcdir, args.mapdir, args.beddir, args.countsdir, args.samtools)
