import pysam
import os
import argparse

def demultiplex_bam(input_bam_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    bam = pysam.AlignmentFile(input_bam_path, "rb")
    output_bams = {}
    
    try:
        for read in bam:
            barcode = read.query_name.split(':')[0]
            barcode = barcode.replace('+', '_')  # replace '+' with '_' in the filename
            
            if barcode not in output_bams:
                output_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(input_bam_path))[0].split('_pe')[0]}.{barcode}.sc.bam")
                output_bams[barcode] = pysam.AlignmentFile(output_path, "wb", template=bam)

            output_bams[barcode].write(read)
    
    finally:
        for bam_file in output_bams.values():
            bam_file.close()
        bam.close()

def main():
    parser = argparse.ArgumentParser(description='Demultiplex a BAM file (from sciMET) based on barcodes')
    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output', default='.', help='Output directory for demultiplexed BAM files')
    args = parser.parse_args()
    demultiplex_bam(args.input, args.output)

if __name__ == '__main__':
    main()

