import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Define the argument parser
parser = argparse.ArgumentParser(description="Extract 16S rRNA sequences from GenBank files.")
parser.add_argument("input_directory", help="Input directory containing subdirectories with GenBank files")

# Parse the command-line arguments
args = parser.parse_args()

# Create an output directory to store 16S rRNA FASTQ files
output_directory = os.path.join(args.input_directory, "16s_XTrakTA_output")
os.makedirs(output_directory, exist_ok=True)


def extract_16s_sequences(gbk_file, subdirectory):
    gbk_path = os.path.join(args.input_directory, subdirectory, gbk_file)
    print(f"Processing GenBank file: {gbk_path}")
    gbank = SeqIO.parse(open(gbk_path, "r"), "genbank")

    rRNAs = []
    for genome in gbank:
        for feature in genome.features:
            if feature.type == "rRNA" and "product" in feature.qualifiers:
                product = feature.qualifiers["product"][0]
                if "16S ribosomal RNA" in product:
                    # Check if 'db_xref' qualifier exists, set to 'Unknown' if missing
                    ID = feature.qualifiers.get('db_xref', ["Unknown"])[0]
                    desc = feature.qualifiers.get('locus_tag', ["Unknown"])[0]
                    seq = feature.extract(genome.seq)

                    # Generate mock quality scores (e.g., all quality scores set to a constant value)
                    #mock_quality_scores = [30] * len(seq)
                    record = SeqRecord(seq, id=ID, description=desc)
                    #record.letter_annotations["phred_quality"] = mock_quality_scores

                    rRNAs.append(record)

    if rRNAs:
        # Save each 16S rRNA sequence as an individual FASTQ file
        for i, rRNA_record in enumerate(rRNAs, start=1):
            filename = f"{subdirectory}.16s_{i}.fasta"  # Change the extension to .fastq
            output_path = os.path.join(output_directory, filename)
            SeqIO.write([rRNA_record], output_path, "fasta")  # Write in FASTQ format


# Iterate through subdirectories in the input directory
for subdirectory in os.listdir(args.input_directory):
    subdirectory_path = os.path.join(args.input_directory, subdirectory)
    if os.path.isdir(subdirectory_path):
        print(f"Processing subdirectory: {subdirectory}")
        for gbk_file in os.listdir(subdirectory_path):
            if gbk_file.endswith(".gbk"):
                extract_16s_sequences(gbk_file, subdirectory)

print("16S rRNA sequences extraction and saving in FASTQ format completed.")
