import gzip
import os
import argparse
from collections import defaultdict

def demultiplex_paired_fastq(fq1_path, fq2_path, prefix, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    file_handles = defaultdict(lambda: {'fq1': None, 'fq2': None})
    
    def get_file_handles(barcode):
        if file_handles[barcode]['fq1'] is None:
            fq1_file = os.path.join(output_dir, f"{prefix}.{barcode}_R1.fq.gz")
            fq2_file = os.path.join(output_dir, f"{prefix}.{barcode}_R2.fq.gz")
            file_handles[barcode]['fq1'] = gzip.open(fq1_file, 'wt')
            file_handles[barcode]['fq2'] = gzip.open(fq2_file, 'wt')
        return file_handles[barcode]
    
    def close_all_handles():
        for handles in file_handles.values():
            handles['fq1'].close()
            handles['fq2'].close()
    
    with gzip.open(fq1_path, 'rt') as fq1, gzip.open(fq2_path, 'rt') as fq2:
        while True:
            try:
                r1_lines = [fq1.readline().strip() for _ in range(4)]
                r2_lines = [fq2.readline().strip() for _ in range(4)]
                
                if len(r1_lines) < 4:
                    break  
                
                read_name = r1_lines[0]
                if read_name=='':
                    break

                parts = read_name.split(':')
                
                if len(parts) < 2:
                    print(f"Warning: read name does not contain a barcode: {read_name}")
                    continue
                
                barcode = parts[0].split('@')[1].replace('+','_')
                
                handles = get_file_handles(barcode)
                
                handles['fq1'].write('\n'.join(r1_lines) + '\n')
                handles['fq2'].write('\n'.join(r2_lines) + '\n')
            
            except StopIteration:
                break

    close_all_handles()

def main():
    parser = argparse.ArgumentParser(description="Demultiplex PE fastq files that contain multiple cells into single-cell fastq")
    parser.add_argument('--fq1', required=True, help='Path to the input fastq R1')
    parser.add_argument('--fq2', required=True, help='Path to the input fastq R2')
    parser.add_argument('--outdir', default='.', help='Directory to output demultiplexed fastq files')
    args = parser.parse_args()
   
    pattern="_val_"
    prefix = args.fq1.split(pattern)[0]
    demultiplex_paired_fastq(args.fq1, args.fq2, prefix, args.outdir)

if __name__ == "__main__":
    main()

