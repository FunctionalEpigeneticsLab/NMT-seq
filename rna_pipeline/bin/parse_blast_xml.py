import xml.etree.ElementTree as ET
import sys

def process_file(file_path):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Dictionary to hold query definitions and their hits
    queries = {}

    # Iterate through each Iteration element in the XML
    for iteration in root.findall('.//Iteration'):
        query_def = iteration.find('Iteration_query-def').text if iteration.find('Iteration_query-def') is not None else ""
        hits = [hit.find('Hit_def').text.replace('&gt;', ' ') for hit in iteration.findall('.//Hit') if hit.find('Hit_def') is not None]
        
        # Check if hits list is empty and assign "No hits found" if true
        if not hits:
            hits = ["No hits found"]
        
        queries[query_def] = hits

    # Print the results
    for query_def, hits in queries.items():
        print(f"{query_def}\t{'\t'.join(hits)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    process_file(file_path)

