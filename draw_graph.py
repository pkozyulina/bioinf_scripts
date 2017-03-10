from Bio import Seq
import collections, argparse

class Graph:
    def __init__(self):
        self.graph = collections.defaultdict(int)
        self.contigs = collections.defaultdict(int)
        self.k = 55

    # Creat a complementary chain
    def twin(self, read):
        return Seq.reverse_complement(read)

    # Check input file format
    def is_fasta(self, file):
        return file.readline().strip().startswith(">")

    # Splits sequence to kmers
    def split_read(self, seq):
        for i in range(len(seq) - self.k + 1):
            yield seq[i:i + self.k]

    # Generates fasta reads
    def split_read_fasta(self, file):
        read = ""
        cnt = 0
        for line in file:
            current_line = line.strip()
            if not current_line.startswith(">"):
                read += current_line
                continue
            current_line = ""
            if read != "":
                yield read
                read = ""
        yield read

    # Generates fastq reads
    def split_read_fastq(self, file):
        for line in file:
            read = line.strip()
            next(file)
            next(file)
            next(file)
            yield read

    # Generates next or previous step of a graph possibilities
    def fw(self, km):
        for x in 'ACGT':
            yield km[1:]+x

    def bw(self, km):
        for x in 'ACGT':
            yield x + km[:-1]

    # Parses all reads in input file file and creates a dict with kmers
    def fileparse(self, file, k=55):
        self.k = k
        reads = self.split_read_fasta(file) if self.is_fasta(file) else self.split_read_fastq(file)
        for read in reads:
            for km in self.split_read(read):
                self.graph[km] += 1
            read = self.twin(read)
            for km in self.split_read(read):
                self.graph[km] += 1
        print('Done parsing file')

    # Merges list to string contig
    def contig_to_string(self, contig):
        return contig[0] + ''.join(x[-1] for x in contig[1:])

    # Get linear contig out of graph
    def get_contig(self, start_km):
        contig_fw = [start_km]
        while True:
            current_km = contig_fw[-1]
            if sum(x in self.graph for x in self.fw(current_km)) != 1:
                break # stops if a joint consists of multiple possible paths or no paths
            cand = [x for x in self.fw(current_km) if x in self.graph][0]
            if cand == start_km or cand == self.twin(start_km):
                break  # break out of cycles or mobius contigs
            #if cand == twin(current_km):
             #   break  # break out of hairpins
            if sum(x in self.graph for x in self.bw(cand)) != 1:
                break
            contig_fw.append(cand)
        return contig_fw

    # Creates a complete contig by combining fw and bw parts
    def combine_fw_and_bw(self, start_km):
        contig_fw = self.get_contig(start_km)
        contig_bw = self.get_contig(self.twin(start_km))
        contig = [self.twin(x) for x in contig_bw[-1:0:-1]] + contig_fw
        coverage = 0
        for con in contig:
            if con in self.graph:
                coverage += self.graph[con]
        return self.contig_to_string(contig), coverage, contig

    # Finds all possible contigs
    def find_all_contigs(self):
        checked = set()
        for km in self.graph:
            if km not in checked:
                current_con, coverage, done_contig_list = self.combine_fw_and_bw(km)
                for con in done_contig_list:
                    checked.add(con)
                self.contigs[current_con] = coverage
                self.contigs[self.twin(current_con)] = coverage
        print('Done getting all contigs')

    # Creates a .dot file with starts and ends of each contig
    def draw_dot(self, output_file):
        with open(output_file, "w") as outp:
            outp.write("digraph {\n")
            for contig in self.contigs:
                start = contig[0:(self.k-1)]
                end = contig[(len(contig) - (self.k-1)):len(contig)]
                length_c = len(contig[self.k:len(contig)])
                if length_c != 0:
                    coverage = self.contigs[contig]/length_c
                else:
                    coverage = self.contigs[contig]
                outp.write('%s -> %s[label="C %i L %i"]\n' % (start, end, coverage, length_c))
            outp.write("}")
        print('Done drawing graph!!!')



# Read input file and create a dict
def main(input_file, output_file, k):
    contigs = Graph()
    with open(input_file, "r") as fq:
        contigs.fileparse(fq, k)
    contigs.find_all_contigs()
    contigs.draw_dot(output_file)



######## SCRIPT BODY #################

parser = argparse.ArgumentParser(description='Generate graph from sequence reads and create a .dot file')
parser.add_argument('-k', '--kmer',  default=55,  help='Set the k-mer size', type=int, metavar='INT')
parser.add_argument('-o', '--output', help='Output file name', required=True, metavar='FILE')
parser.add_argument('-i', '--input', help='Input file name', required=True, metavar='FILE')

args = parser.parse_args()

if __name__ == "__main__":
    main(args.input, args.output, args.kmer)


