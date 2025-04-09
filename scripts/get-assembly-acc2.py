import os
import sys
import argparse
import csv
import time
import re
import requests
import concurrent.futures
import threading
import signal
from tenacity import retry, stop_after_attempt, wait_exponential
import multiprocessing
# avoid seeing openpyxl warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")


API_KEY = os.environ.get("NCBI_API_KEY")
HEADERS = {}

NUM_CORES = multiprocessing.cpu_count()
THREADS = min(10, NUM_CORES * 2)  # Use 2x cores but cap at 10

# Semaphore to limit requests to 10 per second
semaphore = threading.Semaphore(10)
lock = threading.Lock()  # Prevents race conditions when printing

# global var to track if script should exit
should_exit = False

# handle ctrl-c
def signal_handler(sig, frame):
    global should_exit
    print("\nReceived Ctrl+C. Shutting down, this may take a few minutes...")
    should_exit = True

# register singal handler
signal.signal(signal.SIGINT, signal_handler)


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

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=2, min=4, max=10))
def retrieve_assembly_accession(identifier):
    """
    Retrieve the NCBI Assembly accession using the REST API.
    """
    global should_exit
    if should_exit:
        return None  # Stop processing if exiting

    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/sequence_accession/{identifier}/sequence_assemblies"
    params = {"api_key": API_KEY}
    
    with semaphore:  # Ensures no more than 10 requests per second
        try:
            response = requests.get(url, params=params, headers=HEADERS, timeout=10)
            response.raise_for_status() # raise an exception for HTTP errors (e.g. 404, 500)
            
            if response.status_code == 200:
                data = response.json()
                if "accessions" in data:
                    acc = data["accessions"]
                    if acc:
                        return acc[0]
            else:
                with lock:
                    print(f"Error retrieving {identifier}: HTTP {response.status_code} - {response.text}")
        except requests.RequestException as e:
            with lock:
                print(f"Request error for {identifier}: {e}")
        time.sleep(0.1)  # Enforce rate limit
    return None


def retrieve_acc(acc_col):
    """Retrieve the GenBank Assembly accession for a given list of accessions."""
    all_accs = set()
    failure_reasons = []

    if '(' in acc_col or ')' in acc_col:
        failure_reasons.append("parentheses")
    else:
        ids = extract_accs(acc_col)
        if not ids:
            failure_reasons.append('no_accession')

        for id in ids:
            for i in range(1, 21):
                identifier = f"{id}.{i}" if '.' not in id else id
                if should_exit:
                    return "", ["interrupted"]
                result = retrieve_assembly_accession(identifier)
                if result:
                    all_accs.add(result)

    if not all_accs:
        failure_reasons.append('no_assembly')
        return "", failure_reasons

    # Separate GCA and GCF accessions by base accession and version
    acc_versions = {"GCA": {}, "GCF": {}}
    for acc in all_accs:
        match = re.match(r"^(GCF|GCA)_(\d+)\.(\d+)$", acc)
        if match:
            prefix, base, version = match.groups()
            base_acc = f"{prefix}_{base}"
            version = int(version)
            if base_acc not in acc_versions[prefix] or version > acc_versions[prefix][base_acc]:
                acc_versions[prefix][base_acc] = version

    # Prefer highest GCF version, if available
    if acc_versions["GCF"]:
        best_gcf = max(
            acc_versions["GCF"].items(), key=lambda kv: kv[1]
        )  # (base_acc, version)
        return f"{best_gcf[0]}.{best_gcf[1]}", []

    # Else fall back to highest GCA version
    if acc_versions["GCA"]:
        best_gca = max(
            acc_versions["GCA"].items(), key=lambda kv: kv[1]
        )  # (base_acc, version)
        return f"{best_gca[0]}.{best_gca[1]}", []

    failure_reasons.append("no_valid_accession")
    return "", failure_reasons

def find_assembly_accessions(row, n_found):
    """Find the GenBank Assembly accession corresponding to the given GenBank/RefSeq IDs."""
    genbank_assembly_id, genbank_failures = retrieve_acc(row["Virus GENBANK accession"])
    if genbank_assembly_id:
        n_found += 1
        row["GenBank Assembly ID"] = genbank_assembly_id
    row["GenBank Failures"] = ";".join(genbank_failures)
    return row, n_found


def write_directsketch_output(row, basename, ds=None, curated_ds=None):
    virus_name = f"{row.get('Virus name(s)', '')} {row.get('Virus isolate designation', '')}".strip().replace(',', ';')
    gb_col = row.get("Virus GENBANK accession", "")
    gb_assembly = row.get("GenBank Assembly ID")

    if ds and gb_assembly:
        name = f"{gb_assembly} {virus_name}"
        ds.write(f"{row['GenBank Assembly ID']},{name}\n")
    elif curated_ds and gb_col:
        vmr_acc = f"{basename}_{row.get('Species Sort', '')}_{row.get('Isolate Sort', '')}"
        curated_name = f"{vmr_acc} {virus_name}".strip()
        dl_filename = f"{vmr_acc}.fna.gz"

        gb_acc = extract_accs(gb_col)
        dl_links = [
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gba}&rettype=fasta&retmode=text"
            for gba in gb_acc if gba
        ]

        range_info = ""
        if "(" in gb_col and ")" in gb_col:
            range_info = gb_col[gb_col.find("(")+1:gb_col.find(")")]
            range_info = range_info.replace('.', '-')

        if dl_links:
            curated_ds.write(f"{vmr_acc},{curated_name},DNA,,{dl_filename},{';'.join(dl_links)},{range_info}\n")


def main(args):
    global should_exit
    # Need to API key for 10 NCBI requests/s
    if not API_KEY:
        print("Please set the NCBI_API_KEY environment variable.")
        sys.exit(1)
    
    processed = load_processed_rows(args.existing_vmr) if args.existing_vmr else {}
    num_processed = len(processed)
    if num_processed:
        print(f"loaded {num_processed} pre-processed rows.")

    if args.input_vmr.endswith('xlsx'):
        import pandas as pd
        vmr = pd.read_excel(args.input_vmr, sheet_name=args.sheet_name)
        vmr_tsv = args.input_vmr.replace('.xlsx', '.tsv')
        vmr.to_csv(vmr_tsv, index=False, sep='\t')
        print(f"Converted '{args.input_vmr}' to '{vmr_tsv}'")
        if args.only_convert:
            sys.exit(0)
        args.input_vmr = vmr_tsv

    ds = None
    curated_ds = None
    if args.output_directsketch:
        basename = args.output_vmr.split('.tsv')[0]
        ds = open(f"{basename}.gbsketch.csv", 'w')
        ds.write("accession,name\n")
        curated_ds = open(f"{basename}.urlsketch.csv", 'w')
        curated_ds.write("accession,name,moltype,md5sum,download_filename,url,range\n")

    with open(args.input_vmr, 'r') as inF, open(args.output_vmr, 'w', newline='') as out_acc:
        reader = csv.DictReader(inF, delimiter='\t')
        fieldnames = reader.fieldnames + ["GenBank Assembly ID", "GenBank Failures"]
        writer = csv.DictWriter(out_acc, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        n_found = 0
        n_to_search = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as executor:
            future_to_row = {}
            try:
                for row in reader:
                    if should_exit:
                        break # stop processing if exiting
                    gb_col = row["Virus GENBANK accession"]
                    # init these rows for unprocessed accs
                    row["GenBank Assembly ID"] = ""
                    row["GenBank Failures"] = ""

                    if gb_col in processed.keys():
                        processed_row = processed[gb_col]
                        writer.writerow(processed_row)
                        if args.output_directsketch:
                            write_directsketch_output(processed_row, basename, ds, curated_ds)
                        continue
                    n_to_search += 1

                    # Submit row for processing in parallel
                    future = executor.submit(find_assembly_accessions, row, n_found)
                    future_to_row[future] = row

                print(f"Loaded all rows, found {n_to_search} rows to search.")
                search_count = 0

                # Collect results as they finish
                for future in concurrent.futures.as_completed(future_to_row):
                    if should_exit:
                        break  # Stop processing if exiting

                    row, n_found = future.result()
                    
                    writer.writerow(row)
                    search_count += 1
                    if search_count % 100 == 0:
                        print(f"Searched {search_count}/{n_to_search} accessions...")
                        sys.stdout.flush() # ensure progress updates are visible
                        out_acc.flush() # flush every 100 rows
                        if args.output_directsketch:
                            ds.flush()
                            curated_ds.flush()

                    # Handle directsketch output
                    if args.output_directsketch:
                        write_directsketch_output(row, basename, ds, curated_ds)

            except KeyboardInterrupt:
                print("\nInterrupted. Shutting down workers...")
                executor.shutdown(wait=False, cancel_futures=True)
                sys.exit(1)

    if args.output_directsketch:
        ds.close()
        curated_ds.close()

    if not should_exit:
        print(f"Processing complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-vmr", default="inputs/VMR_MSL40.v1.20250307.xlsx")
    parser.add_argument("-o", "--output-vmr", default='inputs/VMR_MSL40.v1.20250307.acc.tsv')
    parser.add_argument("-p", "--existing-vmr")
    parser.add_argument("-s", "--sheet-name", default='VMR MSL40')
    parser.add_argument("-c", "--only-convert", action='store_true', default=False, help="Only convert excel to tsv, then exit.")
    parser.add_argument("--output-directsketch", action="store_true", help="Output directsketch gbsketch file and urlsketch file(s).")
    args = parser.parse_args()

    main(args)

