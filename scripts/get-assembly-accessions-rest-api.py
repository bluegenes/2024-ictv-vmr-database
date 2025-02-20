import os
import sys
import argparse
import csv
import time
import requests

# Set API key for NCBI access
API_KEY = os.environ.get("NCBI_API_KEY")
if not API_KEY:
    print("Please set the NCBI_API_KEY environment variable.")
    sys.exit(1)
# HEADERS = {"Accept": "application/json"} # explicitly specify we want JSON
HEADERS = {}


def load_processed_rows(output_file):
    """
    Load processed rows from the output file.
    """
    processed = {}
    if os.path.exists(output_file):
        with open(output_file, 'r') as outF:
            reader = csv.DictReader(outF, delimiter='\t')
            for row in reader:
                if row["Virus GENBANK accession"] and row["GenBank Assembly ID"]:
                    acc = row["Virus GENBANK accession"]
                    processed[acc] = row
    return processed


def retrieve_assembly_accession(identifier):
    """
    Retrieve the NCBI Assembly accession using the REST API.
    """
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/sequence_accession/{identifier}/sequence_assemblies"
    params = {"api_key": API_KEY}
    try:
        response = requests.get(url, params=params, headers=HEADERS)
        
        if response.status_code == 200:
            data = response.json()
            if "accessions" in data:
                acc = data["accessions"]
                if acc:
                    print(f"Found {identifier}: {acc[0]}")
                    return acc[0]
        else:
            print(f"Error retrieving {identifier}: HTTP {response.status_code} - {response.text}")
    except requests.RequestException as e:
        print(f"Request error for {identifier}: {e}")
        return None
    time.sleep(0.1)  # Limit requests to 10 per second


# def extract_accs(id_col):
#     """Split multiple NCBI accessions; extract segment accessions."""
#     split_ids = []
#     if id_col:
#         ids = id_col.split(";")
#         for id in ids:
#             if ':' in id:
#                 acc_only = id.split(':')[1].strip()
#                 split_ids.append(acc_only)
#             else:
#                 split_ids.append(id.strip())
#     return split_ids

def extract_accs(id_col):
    "split muliple NCBI accs; split names from segment accessions"
    split_ids = []
    if id_col:
        ids = id_col.split(";")
        for id in ids:
            if '(' in id:
                acc_only = id.split('(')[0]
                acc_only = acc_only.strip()
            if ':' in id:
                acc_only = id.split(':')[1]
                acc_only = acc_only.strip()
                split_ids.append(acc_only)
            else:
                id = id.strip()
                split_ids.append(id)
    return split_ids



def retrieve_acc(acc_col):
    """Retrieve the GenBank Assembly accession for a given list of accessions."""
    all_accs = set()
    acc = ""
    failure_reasons = []
    
    # ignore ranged accessions here
    if '(' in acc_col or ')' in acc_col:
        failure_reasons.append("parentheses")
    else:
        ids = extract_accs(acc_col)
        if not ids:
            failure_reasons.append('no_accession')
    
        for id in ids:
            try:
                if '.' not in id:
                    for i in range(1, 21):
                        acc = retrieve_assembly_accession(f"{id}.{i}")
                        if acc:
                            break
                else:
                    acc = retrieve_assembly_accession(id)
            except requests.RequestException:
                failure_reasons.append("retrieval")
                continue
            if not acc:
                failure_reasons.append('no_assembly')
            else:
                all_accs.add(acc)
                if len(all_accs) > 1:
                    failure_reasons.append("multiple_acc")
                    acc = ""
    return acc, failure_reasons


def find_assembly_accessions(row, n_found):
    """Find the GenBank Assembly accession corresponding to the given GenBank/RefSeq IDs."""
    genbank_assembly_id, genbank_failures = retrieve_acc(row["Virus GENBANK accession"])
    if genbank_assembly_id:
        n_found += 1
        row["GenBank Assembly ID"] = genbank_assembly_id
    row["GenBank Failures"] = ";".join(genbank_failures)
    return row, n_found


def main(args):
    if args.input_vmr.endswith('xlsx'):
        import pandas as pd
        vmr = pd.read_excel(args.input_vmr, sheet_name=args.sheet_name)
        vmr_tsv = args.input_vmr.replace('.xlsx', '.tsv')
        vmr.to_csv(vmr_tsv, index=False, sep='\t')
        if args.only_convert:
            sys.exit(0)
        args.input_vmr = vmr_tsv

    if args.output_directsketch:
        basename = args.output_vmr.split('.tsv')[0]
        # open a file for gbsketch (directsketch) output
        ds = open(f"{basename}.gbsketch.csv", 'w')
        ds.write("accession,name\n")
        # open a file for curated directsketch output
        curated_ds = open(f"{basename}.urlsketch.csv", 'w')
        curated_ds.write("accession,name,moltype,md5sum,download_filename,url,range\n")
    
    processed = load_processed_rows(args.existing_vmr) if args.existing_vmr else {}
    n_found = 0
    with open(args.input_vmr, 'r') as inF, open(args.output_vmr, 'w', newline='') as out_acc:
        reader = csv.DictReader(inF, delimiter='\t')
        fieldnames = reader.fieldnames + ["GenBank Assembly ID", "GenBank Failures"]
        writer = csv.DictWriter(out_acc, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for n, row in enumerate(reader):
            gb_col = row["Virus GENBANK accession"]
            row["GenBank Assembly ID"] = ""
            row["GenBank Failures"] = ""
            if gb_col in processed:
                writer.writerow(processed[gb_col])
                continue 
            
            if n % 500 == 0:
                print(f"Processed {n} accessions...")
            
            row, n_found = find_assembly_accessions(row, n_found)
            # seqacc = row["Virus GENBANK accession"]
            # # genbank_assembly_id, genbank_failures = retrieve_acc(seqacc)
            # print(f"Found {genbank_assembly_id} with failures: {genbank_failures}")
            # row["GenBank Assembly ID"] = genbank_assembly_id
            # row["GenBank Failures"] = ";".join(genbank_failures)
            if args.output_directsketch:
                virus_name = f"{row['Virus name(s)']} {row['Virus isolate designation']}".strip().replace(',', ';')
                if row['GenBank Assembly ID']:
                    name = f"{row['GenBank Assembly ID']} {row['Virus name(s)']} {row['Virus isolate designation']}"
                    ds.write(f"{row['GenBank Assembly ID']},{name},\n")
                else:
                    if gb_col == "":
                        # skip anything where we don't have any sequence
                        continue
                    vmr_acc = f"{basename}_{row['Species Sort']}_{row['Isolate Sort']}"
                    curated_name = f"{vmr_acc} {virus_name}".strip()
                    dl_filename = f"{vmr_acc}.fna.gz"
                    # here, handle multiple acc input styles:
                    gb_acc = extract_accs(gb_col)
                
                    dl_links = []
                    for gba in gb_acc:
                        if gba:
                            dl_links.append(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gba}&rettype=fasta&retmode=text")
                    dl_info = ";".join(dl_links)
                    
                    range = ""
                    # if there's a range provided, add it.
                    # NOTE: accs with range should not have multiple accs to join
                    if "(" in gb_col and ")" in gb_col:
                        range = gb_col[gb_col.find("(")+1:gb_col.find(")")]
                        range = range.replace('.', '-')
                    if dl_info is not None:
                        curated_ds.write(f"{vmr_acc},{curated_name},DNA,,{dl_filename},{dl_info},{range}\n")
                        
                    if args.output_directsketch:
                        pass
                        # prep and write urlsketch file
            writer.writerow(row)
    print(f"Found {n_found} assembly accessions.")


def cmdline(sys_args):
    """Command line entry point."""
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input-vmr", default="inputs/VMR_MSL39.v4_20241106.xlsx")
    p.add_argument("-o", "--output-vmr", default='inputs/VMR_MSL39.v4_20241106.acc.tsv')
    p.add_argument("-p", "--existing-vmr", default='inputs/VMR_MSL39.v4_20241106.acc.bak.tsv')
    p.add_argument("-s", "--sheet-name", default='VMR MSL39')
    p.add_argument("-c", "--only-convert", action='store_true', default=False, help="Only convert excel to tsv, then exit.")
    p.add_argument("--output-directsketch", action="store_true", help="Output directsketch gbsketch file and urlsketch file(s).")
    args = p.parse_args()
    return main(args)


if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)