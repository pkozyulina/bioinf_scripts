#!/usr/bin/env python3

import argparse
import collections
from tqdm import tqdm

class Edge:
    show_sequences = True

    def __init__(self, v1, v2): # Edge instance should contain starting and finishing vertices, coverage and edge sequence
        if str(v1)[1:] == str(v2)[:-1]:
            self.edge = v1 + v2[-1]
        else:
            self.edge = str(v1[:-Graph.k]) + str(v2)
        self.vertex_left = v1
        self.vertex_right = v2
        self.coverage = 1

    def print_vertecies(self, k):
        if Edge.show_sequences:
            return '%s -> %s [label = "C %i\nL %i\n %s"];\n' % (self.vertex_left, self.vertex_right, self.coverage/(len(self.edge) - k), (len(self.edge) - k), self.edge)
        return '%s -> %s [label = "C %i\nL %i"];\n' % (self.vertex_left, self.vertex_right, self.coverage /(len(self.edge) - k), (len(self.edge) - k))


    def inc_coverage(self, cov=1):
        self.coverage += cov


    def merge(self, following_edge):

        if self.vertex_right == following_edge.vertex_left:
            new_edge = Edge(self.edge, following_edge)
            new_edge.vertex_left = self.vertex_left
            new_edge.vertex_right = following_edge.vertex_right
            new_edge.inc_coverage(self.coverage + following_edge.coverage - 2)

        return new_edge


    def get_vertex(self, direction):
        if direction == 0:
            return str(self.vertex_left)
        else:
            return str(self.vertex_right)


    def __eq__(self, other):
        return self.edge == other.edge

    def __getitem__(self, i):
        return self.edge[i]

    def __hash__(self):
        hashh = 200
        for el in self.edge:
            hashh += hashh*ord(el)
        return hashh

    def __len__(self):
        return len(self.edge)

    def __str__(self):
        return str(self.edge)

    def __repr__(self):
        return str(self.edge)


class Vertex:

    def __init__(self, seq):
        self.vertex = seq
        self.edges = [{}, {}]


    def add_edge(self, other, direction):

        if direction == 0:
            edge = str(other + self[-1])
        else:
            edge = str(self + other[-1])

        if edge in self.edges[direction]:
            self.edges[direction][edge].inc_coverage()
            return

        if direction == 0:
            self.edges[direction][edge] = other.edges[1][edge]
        else:
            self.edges[direction][edge] = Edge(self, other)


    def extend_edge(self, new_edge, old_edge, direction):
            self.edges[direction].pop(str(old_edge))
            self.edges[direction][str(new_edge)] = new_edge


    def compress(self):
        if len(self.edges[0]) != 1 or len(self.edges[1]) != 1:
            return None

        in_edge = list(self.edges[0].values())[0]
        out_edge = list(self.edges[1].values())[0]
        new_edge = in_edge.merge(out_edge)

        return new_edge


    def is_edge(self, direction):
        return self.edges[direction] != {}

    def __add__(self, other):
        return self.vertex + other

    def __eq__(self, other):
        return self.vertex == other

    def __getitem__(self, i):
        return self.vertex[i]

    def __hash__(self):
        hashh = 500
        for el in self.vertex:
            hashh += hashh*ord(el)
        return hashh

    def __str__(self):
        return str(self.vertex)

    def __repr__(self):
        return str(self.vertex)

    def __len__(self):
        return len(self.vertex)


class Graph:
    k = None
    show_vertex = True

    def __init__(self):
        self.graph = collections.defaultdict(Vertex)

    def is_in_graph(self, seq):
        if seq not in self.graph:
            self.graph[seq] = Vertex(seq)

    def add_edge(self, seq1, seq2):
        self.is_in_graph(seq1)
        self.is_in_graph(seq2)

        self.graph[seq1].add_edge(self.graph[seq2], 1)
        self.graph[seq2].add_edge(self.graph[seq1], 0)


    def split_read(self, seq):
        for i in range(len(seq) - Graph.k):
            yield seq[i:i + Graph.k], seq[i+1:i + 1 + Graph.k]


    def add_seq(self, read): # Adds edges between all k-mers in the sequence
        for kmer1, kmer2 in self.split_read(read):
            self.add_edge(kmer1, kmer2)


    def compress(self):
        to_delete = []  # List of redundant vertices

        for kmer, vertex in self.graph.items():
            new_edge = vertex.compress()
            if new_edge:
                to_delete.append(kmer)
                self.graph[str(new_edge.vertex_left)].extend_edge(new_edge, list(vertex.edges[0].values())[0], 1)
                self.graph[str(new_edge.vertex_right)].extend_edge(new_edge, list(vertex.edges[1].values())[0], 0)

        for kmer in to_delete:
            del(self.graph[kmer])


    def save_dot(self, output_file):
        printed_vertecies = set()
        with open(output_file.name, "w") as outp:
            outp.write("digraph {\n")
            if not Graph.show_vertex:
                outp.write('node[label=""];')
            for seq, vertex in self.graph.items():
                if seq not in printed_vertecies:
                    for in_edge in vertex.edges[0]:
                        outp.write(vertex.edges[0][in_edge].print_vertecies(Graph.k))
                    printed_vertecies.add(seq)
            outp.write("}")


    def __len__(self):
        return len(self.graph)

    def __str__(self):
        return str(self.graph)



complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def reverse_complement(seq):
    return ''.join(complement[nt] for nt in seq[::-1])


def read_fastq(f):
    for line in f:
        name = line.strip()
        seq = next(f).strip()
        next(f)
        next(f)
        yield name, seq


def read_fasta(f):
    name = None
    seq = None
    for line in f:
        if line.startswith('>'):
            if name:
                yield name, seq
            name = line.lstrip('>').strip()
            seq = ''
        else:
            seq += line.strip()
    yield name, seq


def read(f):
    if f.name.endswith('a'):
        return read_fasta(f)
    else:
        return read_fastq(f)


def main():
    parser = argparse.ArgumentParser(description='De Bruijn graph')
    parser.add_argument('-i', '--input', help='Input fastq', metavar='File',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-k', help='k-mer size (default: 55)', metavar='Int',
                        type=int, default=55)
    parser.add_argument('-o', '--output', help='Output dot', metavar='File',
                        type=argparse.FileType('w'), required=True)
    parser.add_argument('-c', '--compress', help='Shrink graph', action='store_true')
    parser.add_argument('--vertex', help='Show vertex sequences', action='store_true')
    parser.add_argument('--edge', help='Show edge sequences', action='store_true')
    args = parser.parse_args()

    Graph.k = args.k

    Graph.show_vertex = args.vertex
    Edge.show_sequences = args.edge

    graph = Graph()
    for name, seq in tqdm(read(args.input)):
        graph.add_seq(seq)
        graph.add_seq(reverse_complement(seq))

    if args.compress:
        graph.compress()

    graph.save_dot(args.output)


if __name__ == '__main__':
    main()
