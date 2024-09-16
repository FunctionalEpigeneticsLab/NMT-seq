import argparse
import os

def parse_flagstat(flagstat_file):
    with open(flagstat_file) as f:
        for line in f:
            if "mapped (" in line and "primary" in line:
                count = int(line.split(" + ")[0])
                return count
    return 0

def main():
    parser = argparse.ArgumentParser(description="Extract the number of uniquely mapped reads from flagstat files.")
    parser.add_argument("--flagstat1")
    parser.add_argument("--flagstat2")
    parser.add_argument("--flagstat3")
    parser.add_argument("--output_prefix", help="Output prefix for the summary file")

    args = parser.parse_args()
    
    initial_flagstat=args.flagstat1
    deduplicated_flagstat=args.flagstat2
    filtered_flagstat=args.flagstat3
    cellname=args.output_prefix

    initial_unique = parse_flagstat(initial_flagstat)
    deduplicated_unique = parse_flagstat(deduplicated_flagstat)
    filtered_unique = parse_flagstat(filtered_flagstat)

    output_file = f"{args.output_prefix}.mapping_summary.tsv"

    with open(output_file, "w") as out_f:
        #out_f.write(f"cell\tuniq_map\tuniq_dedup\tuniq_filtered\n")
        out_f.write(f"{cellname}\t{initial_unique}\t{deduplicated_unique}\t{filtered_unique}\n")

if __name__ == "__main__":
    main()

