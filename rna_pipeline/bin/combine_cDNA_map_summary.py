import os
import glob
import pandas as pd
import argparse

wd="/staging/leuven/stg_00064/Xinran/project/snNMT_PD"
os.chdir(wd)

def to_numeric(summary):
    for col in summary.columns:
        if summary[col].dtype == 'object':
            first_non_null = summary[col].dropna().iloc[0]
            if '%' in first_non_null:
                summary[col] = summary[col].str.rstrip('%').astype(float)
            else:
                try:
                    summary[col] = summary[col].str.replace(',', '').astype(float)
                except ValueError:
                    pass
    return summary

def get_combined_data(tsv_files):
    combined_data = pd.DataFrame()
    for tsv_file in tsv_files:
        file_path = tsv_file
        print('Reading '+file_path)
        current_data = pd.read_csv(file_path, sep='\t', header=0) #skiprows=[0]
        current_data['Sample'] = current_data['Sample'].apply(lambda x: x.split('_S')[0])
        combined_data = pd.concat([combined_data, current_data], axis=0)
    combined_data = to_numeric(combined_data)
    return combined_data

def get_combined_data_all(tsv_files):
    combined_data = pd.DataFrame()
    for tsv_file in tsv_files:
        file_path = os.path.join(folder_path, tsv_file)
        print('Reading '+file_path)
        current_data = pd.read_csv(file_path, sep='\t', header=0)
        current_data['Sample'] = current_data['Sample'].apply(lambda x: x.split('_S')[0])
        combined_data = pd.concat([combined_data, current_data], axis=0)
    return combined_data

def write_summary(folder_pattern, formatted_date):
    tsv_files = []
    for folder_path in glob.glob(folder_pattern):
        tsv_files.extend([os.path.join(folder_path, f) for f in os.listdir(folder_path) if (f.endswith('.tsv') and f.startswith('RNAalign.batch.stats'))])
    print('Number of mapping summary files found:',len(tsv_files))

    #for f in tsv_files:
    #    print(f)
    combined_data = get_combined_data(tsv_files)
    print('Number of libraries in combined summary:', combined_data.shape[0])

    map_summary_file = os.path.join("RNA_"+formatted_date+"/combined_RNAalign.batch.stats."+formatted_date+".tsv")
    combined_data.to_csv(os.path.join(wd, map_summary_file), sep='\t', index=False)
    print('Written to:', os.path.join(wd, map_summary_file))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Example program with a -d argument')
    parser.add_argument('-d', '--date', type=str, default="", help='Input date (prefix) for analysis (default: "")')
    args = parser.parse_args()
    folder_pattern = "./*/summary"
    write_summary(folder_pattern, args.date)

