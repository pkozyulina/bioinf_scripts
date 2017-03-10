from Bio import Seq
import collections, argparse

parser = argparse.ArgumentParser(description='Generate graph from sequence reads and create a .dot file')
parser.add_argument('-k', '--kmer',  default=55,  help='Set the k-mer size', type=int, metavar='INT')
parser.add_argument('-o', '--output', help='Output file name', required=True, metavar='FILE')
parser.add_argument('-i', '--input', help='Input file name', required=True, metavar='FILE')

args = parser.parse_args()


# Creat a complementary chain
def twin(read):
    return Seq.reverse_complement(read)

# Check input file format
def filetype(file):
    return str(file.readline).strip().startswith(">")

# Splits sequence to kmers
def split_read_to_kmer(seq,k):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

# Generates next or previous step of a graph possibilities
def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]

# Parses all reads in input file file and creates a dict with kmers
def fileparse(file, k=55, limit=1):
    d = collections.defaultdict(int)
    if filetype(file):
        for line in file:
            read = line.strip()
            for km in split_read_to_kmer(read, k):
                d[km] += 1
            read = twin(read)
            for km in split_read_to_kmer(read, k):
                d[km] += 1
    else:
        for line in file:
            _ = line
            read = next(file).strip()
            next(file)
            next(file)
            for km in split_read_to_kmer(read, k):
                d[km] += 1
            read = twin(read)
            for km in split_read_to_kmer(read, k):
                d[km] += 1

    return d

# Merges list to string contig
def contig_to_string(contig):
    return contig[0] + ''.join(x[-1] for x in contig[1:])

# Get linear contig out of graph
def get_contig(graph, start_km):
    contig_fw = [start_km]
    while True:
        current_km = contig_fw[-1]
        if sum(x in graph for x in fw(current_km)) != 1:
            break # stops if a joint consists of multiple possible paths or no paths
        cand = [x for x in fw(current_km) if x in graph][0]
        if cand == start_km or cand == twin(start_km):
            break  # break out of cycles or mobius contigs
        #if cand == twin(current_km):
         #   break  # break out of hairpins
        if sum(x in graph for x in bw(cand)) != 1:
            break
        contig_fw.append(cand)
    return contig_fw

# Creates a complete contig by combining fw and bw parts
def combine_fw_and_bw(graph, start_km):
    contig_fw = get_contig(graph, start_km)
    contig_bw = get_contig(graph, twin(start_km))
    contig = [twin(x) for x in contig_bw[-1:0:-1]] + contig_fw
    return contig_to_string(contig) #, contig

# Finds all possible contigs
def find_all_contigs(graph):
    contigs = set()
    checked = set()
    for km in graph:
        if km not in checked:
            checked.add(km)
            current_con = combine_fw_and_bw(graph, km)
            contigs.add(current_con)
            contigs.add(twin(current_con))
        else:
            pass
    return contigs

# Creates a .dot file with starts and ends of each contig
def draw_dot(contigs, output_file, k):
    with open(output_file, "w") as outp:
        outp.write("digraph {\n")
        for contig in contigs:
            start = contig[0:(k-1)]
            end = contig[(len(contig) - (k-1)):len(contig)]
            L = len(contig[k:(len(contig) - (k-1))])
            outp.write('%s -> %s[label="C L %i"]\n' % (start, end, L))
        outp.write("}")


# Read input file and create a dict
def main(input_file, output_file, k):
    with open(input_file, "r") as fq:
        kmers = fileparse(fq, k)
    contigs = find_all_contigs(kmers)
    draw_dot(contigs, output_file, k)



######## SCRIPT BODY #################
if __name__ == "__main__":
    main(args.input, args.output, args.kmer)

