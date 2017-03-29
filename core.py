import argparse, collections
from levenshtein import levenshtein
from Bio import SeqIO


# Parsing fastq
def parse_fastq(input_file, quality):
    for record in SeqIO.parse(input_file, "fastq"):
        if min(record.letter_annotations["phred_quality"]) >= quality:
            yield record.seq

# Creating a reference dictionary
def parse_reference(ref_file):
    record_dict = collections.defaultdict(str)
    for rec in SeqIO.parse(ref_file, "fasta"):
        record_dict[rec.id] = rec.seq
    return record_dict



def build_alignment(input_file, ref_file_norm, ref_file_mut, quality):
    # here we will have a dict with gene names and how many reads aligns to norm and mutant form
    gene_list = collections.defaultdict(list)
    ref_norm_dict = parse_reference(ref_file_norm)
    ref_mut_dict = parse_reference(ref_file_mut)
    for seq in parse_fastq(input_file, quality):
        min_dist = 1000
        ref_min = ""
        seq_small = seq
        #print('SEQ %s' % (seq))
        for ref_norm in ref_norm_dict:
            dist = levenshtein(ref_norm_dict[ref_norm], seq_small)
            #print(dist, ref_norm)
            if dist < min_dist:
                min_dist = dist
                ref_min = ref_norm
        mut_dist = levenshtein(ref_mut_dict[ref_min], seq_small)
        gene_list[ref_min] = [min_dist, mut_dist]
    print(gene_list)




def main():
    parser = argparse.ArgumentParser(description='Clinically relevant somatic mutations caller tool')
    parser.add_argument('-i', '--input', help='Input fastq file', metavar='File', type=argparse.FileType(),
                                    required=True)
    parser.add_argument('-n', '--nref', help='Normal amplicon reference in fasta-format', metavar='File',
                                    type=argparse.FileType(), required=True)
    parser.add_argument('-m', '--mref', help='Mutant amplicon reference in fasta-format', metavar='File',
                                    type=argparse.FileType(), required=True)
    parser.add_argument('-q', '--quality', help='Quality filtering (default: 30)', metavar='Int', type=int,
                                    default=30)
    #parser.add_argument('-o', '--output', help='Output vcf', metavar='File', type=argparse.FileType('w'),
     #                               required=True)
    args = parser.parse_args()

    #print(parse_reference(args.nref))
    print(build_alignment(args.input, args.nref, args.mref, args.quality))


if __name__ == '__main__':
    main()
