#!/usr/bin/env python3

import argparse
import collections
from tqdm import tqdm

class Edge:
    show_sequences = False

    def __init__(self, v1, v2): # Edge instance should contain starting and finishing vertices, coverage and edge sequence
        if str(v1)[1:] == str(v2)[:-1]:
            self.edge = v1 + v2[-1]
        else:
            self.edge = str(v1[:-Graph.k]) + str(v2)
        self.vertex_left = v1
        self.vertex_right = v2
        self.coverage = 1


    def inc_coverage(self, cov=1):
        self.coverage += cov

    def merge(self, following_edge):
        if self.vertex_right == following_edge.vertex_left:
            new_edge = Edge(self.edge, following_edge)
            new_edge.vertex_left = self.vertex_left
            new_edge.vertex_right = following_edge.vertex_right
            new_edge.inc_coverage(self.coverage + following_edge.coverage)

        return new_edge

    def get_vertices(self, direction):
        if direction == 0:
            return self.vertex_left
        else:
            return self.vertex_right

    def get_coverage(self):
        return self.coverage

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


            #  def __repr__(self):
   #     return str(self.edge)


class Vertex:
    show_sequences = False

    def __init__(self, seq):
        self.vertex = seq
        self.edges = [{}, {}]


    def add_edge(self, other, direction):
        if direction == 0:
            edge = Edge(other, self)
        else:
            edge = Edge(self, other)

        if str(edge) not in self.edges[direction]:
            self.edges[direction][str(edge)] = edge
        else:
            self.edges[direction][str(edge)].inc_coverage()


    def extend_edge(self, new_edge, old_edge, direction):
            cov = self.edges[direction].pop(str(old_edge)).coverage
            new_edge.inc_coverage(cov)
            self.edges[direction][str(new_edge)] = new_edge


    def compress(self):
        if len(self.edges[0]) != 1 or len(self.edges[1]) != 1:
            return None

        in_edge = list(self.edges[0].values())[0]
        out_edge = list(self.edges[1].values())[0]
        new_edge = in_edge.merge(out_edge)

        return new_edge


    def get_edges(self, direction):
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

    def __init__(self):
        # Contains all vertices
        self.graph = collections.defaultdict(Vertex)

    def add_edge(self, seq1, seq2):
        # Increases coverage if the edge already exists

        if seq1 in self.graph and seq2 in self.graph:
            self.graph[seq1].add_edge(self.graph[seq2], 1)
            self.graph[seq2].add_edge(self.graph[seq1], 0)
        elif seq1 in self.graph:
            ver2 = Vertex(seq2)
            self.graph[seq1].add_edge(ver2, 1)
            ver2.add_edge(self.graph[seq1], 0)
            self.graph[seq2] = ver2
        elif seq2 in self.graph:
            ver1 = Vertex(seq1)
            self.graph[seq2].add_edge(ver1, 0)
            ver1.add_edge(self.graph[seq2], 1)
            self.graph[seq1] = ver1
        else:
            ver1, ver2 = Vertex(seq1), Vertex(seq2)
            ver1.add_edge(ver2, 1)
            ver2.add_edge(ver1, 0)
            self.graph[seq1] = ver1
            self.graph[seq2] = ver2

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
            for kmer, vertex in self.graph.items():
                if vertex.get_edges(0):
                    in_edge = list(vertex.edges[0].values())[0]
                    if str(in_edge.vertex_left) not in printed_vertecies and kmer not in printed_vertecies:
                        outp.write("%s -> %s [label = C%i];\n" % (in_edge.vertex_left, vertex, in_edge.coverage/(len(in_edge) - Graph.k)))
                        printed_vertecies.add(kmer)

                if vertex.get_edges(1):
                    out_edge = list(vertex.edges[1].values())[0]
                    if str(out_edge.vertex_right) not in printed_vertecies and kmer not in printed_vertecies:
                        outp.write("%s -> %s [label = C%i];\n" % (vertex, out_edge.vertex_right, out_edge.coverage/(len(out_edge) - Graph.k)))
                        printed_vertecies.add(kmer)
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
    Vertex.show_sequences = args.vertex
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
