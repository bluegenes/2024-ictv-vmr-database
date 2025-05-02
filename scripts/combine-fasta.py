import os
import gzip
import csv
import screed
import argparse

def main(args):
    # Collect fasta files from urlsketch_csv
    fastafiles = []
    with open(args.urlsketch_csv) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            fastafiles.append(os.path.join(args.curated_fasta_dir, row['download_filename']))

    # Collect fasta files from gbsketch_csv
    with open(args.gbsketch_csv) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ident = row['name'].split(' ')[0]
            fastafiles.append(os.path.join(args.assembly_fasta_dir, ident + "_genomic.fna.gz"))

    # use seen_accs to fix a bug where we were including both GCA and GCF genomes for some accesions
    seen_accs = set()
    # Write combined fasta and lengths
    with gzip.open(args.combined, 'wt') as combined_file, open(args.lengths, 'w') as lens_file:
        lens_file.write("filename,sequence_name,length\n")
        for ff in fastafiles:
            acc = ff.split('_genomic.fna.gz')[0]
            if acc.startswith("GC"):
                acc = acc.split("_", 1)[1]
            if acc in seen_accs:
                continue
            seen_accs.add(acc)
            with screed.open(ff) as seqfile:
                seen_names = set()
                for record in seqfile:
                    bp = len(record.sequence)
                    if record.name not in seen_names: # handle bug where we appended the same record to the file
                        lens_file.write(f"{os.path.basename(ff)},{record.name},{bp}\n")
                        combined_file.write(f">{record.name}\n{record.sequence}\n")
                        seen_names.add(record.name)

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Combine FASTA files and save sequence length information.")
    p.add_argument("--urlsketch-csv", required=True, help="CSV file with URLs and filenames for curated FASTA files")
    p.add_argument("--gbsketch-csv", required=True, help="CSV file with GenBank sketch data")
    p.add_argument("--combined", required=True, help="Output combined FASTA file (gzipped)")
    p.add_argument("--lengths", required=True, help="Output CSV file with sequence length information")
    p.add_argument("--curated-fasta-dir", required=True, help="Directory containing curated FASTA files")
    p.add_argument("--assembly-fasta-dir", required=True, help="Directory containing assembly FASTA files")

    args = p.parse_args()

    main(args)