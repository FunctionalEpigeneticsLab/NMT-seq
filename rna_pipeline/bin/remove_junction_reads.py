import pysam
from argparse import ArgumentParser
import multiprocessing

def parse_gtf(gtf_file):
    gene_boundaries = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, feature_type, start, end, attributes = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]
            gene_id = [attr for attr in attributes.split(';') if 'gene_id' in attr][0].split('"')[1]


            if feature_type in ['exon', 'intron']:
                key = (chrom, gene_id)
                if key not in gene_boundaries:
                    gene_boundaries[key] = set()
                gene_boundaries[key].add(start)
                gene_boundaries[key].add(end)
    return gene_boundaries

def filter_reads_for_gene(args):
    bam_file, output_file, gene_boundaries, gene_key = args
    bam = pysam.AlignmentFile(bam_file, "rb")
    out_bam = pysam.AlignmentFile(output_file, "wb", template=bam)
    total_reads = junction_reads = output_reads = 0

    boundaries = gene_boundaries[gene_key]
    chrom, gene_id = gene_key

    for read in bam.fetch(chrom):
        total_reads += 1
        if any(boundary - 1 <= read.reference_start < boundary + 1 or boundary - 1 < read.reference_end <= boundary + 1 for boundary in boundaries):
            junction_reads += 1
        else:
            out_bam.write(read)
            output_reads += 1

    summary = {"gene_id": gene_id, "total_reads": total_reads, "junction_reads": junction_reads, "output_reads": output_reads, "temp_output_file": output_file }
    bam.close()
    out_bam.close()
    return summary

def main():
    parser = ArgumentParser(description="Filter reads not spanning exon-intron junctions and generate summary.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file with exon and intron annotations")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    args = parser.parse_args()

    gene_boundaries = parse_gtf(args.gtf)

    # Correcting the preparation of arguments for multiprocessing
    # Use gene_key to extract chromosome and gene_id for filename
    pool_args = [(args.bam, f"{args.output}_{gene_key[1]}.bam", gene_boundaries, gene_key) for gene_key in gene_boundaries.keys()]

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(filter_reads_for_gene, pool_args)

    # merging BAM files
    output_bam_filenames = [result["temp_output_file"] for result in results]
    merged_output_filename = args.output  # final, merged BAM file name
    
    # Merge BAM files and clean up
    #pysam.merge("-f", merged_output_filename, *output_bam_filenames)
    #for filename in output_bam_filenames:
    #    os.remove(filename)
    # Index the merged BAM file
    #pysam.index(merged_output_filename)

    # Aggregate results and write summary
    total_reads = sum(result["total_reads"] for result in results)
    junction_reads = sum(result["junction_reads"] for result in results)
    output_reads = sum(result["output_reads"] for result in results)
    
    summary_file = f"{merged_output_filename}.junctions.summary"
    with open(summary_file, 'w') as f:
        f.write(f"Total reads processed\t{total_reads}\n")
        f.write(f"Reads overlapping with junctions\t{junction_reads}\n")
        f.write(f"Reads in final output (non-junctional)\t{output_reads}\n")

if __name__ == "__main__":
    main()
