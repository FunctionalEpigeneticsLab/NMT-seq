import pysam
from argparse import ArgumentParser

def parse_gtf(gtf_file):
    """
    Parses a GTF file and returns a dictionary with gene-based exon and intron boundaries.
    The dictionary keys are tuples of (chromosome, gene_id), and the values are sets of boundary positions.
    """
    gene_boundaries = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split('\t')
            chrom, feature_type, start, end, attributes = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]
            gene_id = [attr for attr in attributes.split(';') if 'gene_id' in attr][0].split('"')[1]

            if feature_type in ['exon', 'intron']:  # Use your GTF's specific feature types
                key = (chrom, gene_id)
                if key not in gene_boundaries:
                    gene_boundaries[key] = set()
                # Add start and end positions to the set of boundaries for this gene
                gene_boundaries[key].add(start)
                gene_boundaries[key].add(end)
    return gene_boundaries


def filter_reads(bam_file, gene_boundaries, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    #out_bam = pysam.AlignmentFile(output_file, "wb", template=bam)
    junction_bam = pysam.AlignmentFile(output_file, "wb", template=bam)

    total_reads = junction_reads = output_reads = 0

    for read in bam.fetch():
        total_reads += 1
        chrom = read.reference_name
        is_junction_read = False
        for (gene_chrom, gene_id), boundaries in gene_boundaries.items():
            if chrom == gene_chrom and any(boundary - 1 <= read.reference_start < boundary + 1 or boundary - 1 < read.reference_end <= boundary + 1 for boundary in boundaries):
                junction_reads += 1
                junction_bam.write(read)
                is_junction_read = True
                break
        if not is_junction_read:
            #out_bam.write(read)
            output_reads += 1

    bam.close()
    #out_bam.close()
    junction_bam.close()

    # Write summary
    summary_file = f"{output_file}.junctions.summary"
    with open(summary_file, 'w') as f:
        f.write(f"Total processed reads\n{total_reads}\n")
        f.write(f"Reads overlapping with junctions\n{junction_reads}\n")
        f.write(f"Reads in output (non-junctional)\n{output_reads}\n")

def main():
    parser = ArgumentParser(description="Filter reads not spanning exon-intron junctions and generate summary.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file with exon and intron annotations")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    args = parser.parse_args()

    gene_boundaries = parse_gtf(args.gtf)
    
    junction_reads_file = "junction_reads.bam"  # Define the filename for junctional reads
    filter_reads(args.bam, gene_boundaries, args.output)

if __name__ == "__main__":
    main()
