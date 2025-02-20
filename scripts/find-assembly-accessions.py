import os
import sys
import argparse
import csv
from Bio import Entrez

# Set email address for NCBI access
Entrez.email = "ntpierce@ucdavis.edu"

def load_processed_rows(output_file):
    """
    Load processed rows from the output file.

    Args:
        output_file (str): Path to the output TSV file.

    Returns:
        set: A set of processed 'Virus GENBANK accession' values.
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
    Retrieve the NCBI Assemblies dataset identifier for the given RefSeq or GenBank identifier.

    Args:
        identifier (str): The RefSeq or GenBank identifier.

    Returns:
        str: The NCBI Assemblies dataset identifier, or None if no Assembly is found for the given identifier.
    """
    handle = Entrez.esearch(db="assembly", term=identifier)
    record = Entrez.read(handle)
    handle.close()
    if record["Count"] == "0":
        return None
    else:
        assembly_id = record["IdList"][0]
        handle = Entrez.efetch(db="assembly", id=assembly_id, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
        dataset_id = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
        return dataset_id


def extract_accs(id_col):
    "split muliple NCBI accs; split names from segment accessions"
    split_ids = []
    if id_col:
        ids = id_col.split(";")
        for id in ids:
            if ':' in id:
                acc_only = id.split(':')[1]
                acc_only = acc_only.strip()
                split_ids.append(acc_only)
            else:
                id = id.strip()
                split_ids.append(id)
    return split_ids

def retrieve_acc(acc_col):
    """
    Given a list of NCBI accessions, retrieves the corresponding GenBank/RefSeq Assemblies Dataset Accession.
    If no definitive accession can be retrieved, returns an empty string and the failure reason(s) encountered
    during accession retrieval.

    Parameters:
    -----------
    acc_col : str
        A column containing one or more NCBI accession IDs (;-separated) to be searched for corresponding GenBank/RefSeq assembly IDs.

    Returns:
    --------
    acc : str
        The Assemblies dataset ID corresponding to the first accession ID in the input list that can be successfully
        retrieved and does not contain parentheses, or the assembly ID corresponding to multiple IDs if they all
        map to the same assembly ID. If none of the input IDs can be retrieved, returns an empty string.

    failure_reasons: list
         Failures are provided to aid in manual resolution of these issues.
         Options:
            no_accession - no NCBI accession was provided.
            no_assembly  - no corresponding Assemblies dataset exists.
            parentheses  - Parentheses in an id, e.g "LK928904 (2253..10260)", indicate that the viral genome
                           is this sequence range within the host genome identified by the accession. We do not
                           want to mischaracterize the entire host genome as this virus; return a failure.
            retrieval    - some error was encountered while attempting to find the accession via entrez.
                           This may occur due to removal or modification of the accession.
            multiple_acc - multiple assembly accessions were returned via the list of nucleotide accessions.
                           The multiple accessions in any entry should correspond to segments of a segmented virus,
                           all of which should be collected into a single 'Assembly Dataset' Accession we can retrieve.
                           If we encounter multiple accessions, it may indicate an accession removal or modification.
    """
    all_accs = set()
    acc = ""
    failure_reasons = []
    ids = extract_accs(acc_col)
    if not ids:
        failure_reasons=['no_accession']
    for id in ids:
        if "(" in id or ')' in id:
            failure_reasons.append("parentheses")
            continue
        try:
            acc= retrieve_assembly_accession(id)
        except RuntimeError:
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


def find_assembly_accessions(row):
    """
    Find the GenBank Assemblies dataset accession  corresponding to the GenBank/RefSeq IDs in the given row.

    Parameters:
    -----------
    row : dict
        A dictionary containing the key "Virus GENBANK accession".

    Returns:
    --------
    dict
        The input dictionary with two additional keys:
            - "GenBank Assembly ID"
            - "GenBank Failures"
    """
    # Retrieve the Assembly dataset accession, if available
    genbank_assembly_id, genbank_failures = retrieve_acc(row["Virus GENBANK accession"])
    # Add the assembly IDs to the row
    row["GenBank Assembly ID"] = genbank_assembly_id
    row["GenBank Failures"] = ";".join(genbank_failures)
    return row


def main(args):
    # if we have excel format, convert using pandas
    input_vmr = args.input_vmr
    if input_vmr.endswith('xlsx'):
        import pandas as pd
        # ICTV VMR Reference file
        vmr = pd.read_excel(input_vmr, sheet_name=args.sheet_name)
        vmr_tsv = args.input_vmr.split('xlsx')[0] + 'tsv'
        vmr.to_csv(vmr_tsv, index=False, sep='\t') # there are commas in some columns, so use tab
        if args.only_convert:
            sys.exit(0)
        input_vmr = vmr_tsv

    if args.output_directsketch:
        ds = open(args.output_directsketch, 'w')
        ds.write("accession,name,ftp_path\n")
    
    # Load processed rows if the output file exists
    if args.existing_vmr:
        processed = load_processed_rows(args.existing_vmr)

    # read in the file with tsv, loop through to link to assembly dataset information
    with open(input_vmr, 'r') as inF, open(args.output_vmr, 'w', newline='') as out_acc:
        reader = csv.DictReader(inF, delimiter='\t')
        fieldnames = reader.fieldnames + ["GenBank Assembly ID", "GenBank Failures"]
        writer = csv.DictWriter(out_acc, fieldnames=fieldnames, delimiter='\t') # there are commas in some columns, so use tab
        writer.writeheader()
        for n, row in enumerate(reader):
            gb_acc = row["Virus GENBANK accession"]
            if gb_acc in processed.keys():
                processed_row = processed[gb_acc]
                writer.writerow(processed_row)
                continue  # Skip re-processing already processed rows
            if n % 500 == 0:
                print(f"Processed {n} accessions...")
            row = find_assembly_accessions(row)
            if row['GenBank Assembly ID'] and args.output_directsketch:
                name = row['GenBank Assembly ID'] + ' ' +  row['Virus name(s)'] + ' ' + row['Virus isolate designation'] 
                ds.write(f"{row['GenBank Assembly ID']},{name},\n")
            writer.writerow(row)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input-vmr", default= "inputs/VMR_MSL39.v4_20241106.xlsx")
    p.add_argument("-o", "--output-vmr", default='inputs/VMR_MSL39.v4_20241106.acc.tsv')
    p.add_argument("-p", "--existing-vmr", default='inputs/VMR_MSL39.v4_20241106.acc.bak.tsv')
    p.add_argument("-s", "--sheet-name", default='VMR MSL39')
    p.add_argument("-c", "--only-convert", action='store_true', default=False, help="Only convert excel to tsv, then exit.")
    p.add_argument("--output-directsketch", help="optionally, output file for sourmash_plugin_directsketch use")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
