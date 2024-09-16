import gffutils
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Compute gene lengths from a GTF file and output to a TSV file.')
parser.add_argument('-i', '--input', required=True, help='Input GTF file path.')
parser.add_argument('-o', '--output', required=True, help='Output TSV file path.')

args = parser.parse_args()

gtf_file = args.input
output_file = args.output

db_file = gtf_file + '.db'
db = gffutils.create_db(gtf_file, dbfn=db_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True,disable_infer_genes=True,disable_infer_transcripts=True)

# Open the database
#db = gffutils.FeatureDB(db_file)

genes_info = []
# Iterate over gene entries in the GTF database
for gene in db.features_of_type('gene'):
    gene_id = gene.attributes['gene_id'][0]
    gene_name = gene.attributes.get('gene_name', [''])[0] 
    gene_length = gene.end - gene.start + 1
    genes_info.append([gene_id, gene_name, gene_length])


genes_df = pd.DataFrame(genes_info, columns=['GeneID', 'GeneName', 'GeneLength'])

genes_df.to_csv(output_file, sep='\t', index=False)
print(f'Gene lengths have been saved to {output_file}')
