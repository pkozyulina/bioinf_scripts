#!/usr/bin/env python3

import argparse, collections
from tqdm import tqdm


class Edge:
    show_sequences = False

    def __init__(self, v1, v2): # Edge instance should contain starting and finishing vertices, coverage and edge sequence
        self.edge = v1 + v2[-1]
        self.vertex_left = v1
        self.vertex_right = v2
        self.coverage = 1

    def inc_coverage(self, cov=1):
        self.coverage += cov

    def merge(self, following_edge):
        new_edge = Edge(self.edge, following_edge)
        new_edge.vertex_left = self.vertex_left
        new_edge.vertex_right = following_edge.vertex_right
        new_edge.inc_coverage(self.coverage + following_edge.coverage)
        return new_edge

    def get_vertices(self, direction):
        if direction == "in":
            return self.vertex_left
        else:
            return self.vertex_right

    def get_coverage(self):
        return self.coverage

    def __eq__(self, other):
        return self.edge == other

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
    show_sequences = False

    def __init__(self, seq):
        self.vertex = seq
        self.edges = {'in': {}, 'out': {}}


    def add_edge(self, other, direction):
        if direction == 'in':
            edge = Edge(other, self)
        else:
            edge = Edge(self, other)

        if str(edge) not in self.edges[direction]:
            self.edges[direction][str(edge)] = edge
        else:
            self.edges[direction][str(edge)].inc_coverage()


    def extend_edge(self, edge, k, direction):
        if direction == 'in':
            key = str(edge)[len(str(edge))-k-1:]
        if direction == 'out':
            key = str(edge)[:k+1]
        #print(key)
        if key in self.edges[direction]:
            cov = self.edges[direction][key].coverage # IN
            edge.inc_coverage(cov)
            self.edges[direction].pop(key)
            self.edges[direction][str(edge)] = edge
#            print('Added new Edge %s with coverage %i to the Vertex %s' % (edge, edge.coverage, self.vertex))


    def compress(self):
        #print(self.vertex, self.edges)
        if len(self.edges['in']) != 1 or len(self.edges['out']) != 1:
            return False

        for in_edge in self.edges['in']:
            for out_edge in self.edges['out']:
                new_edge = self.edges['in'][in_edge].merge(self.edges['out'][out_edge])

        return new_edge


    def get_edges(self, direction):
        return self.edges[direction]

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


class Graph:
    k = None

    def __init__(self):
        # Contains all vertices
        self.graph = collections.defaultdict(Vertex)

    def add_edge(self, seq1, seq2):
        # Increases coverage if the edge already exists
        ver1, ver2 = Vertex(seq1), Vertex(seq2)
        if seq1 in self.graph and seq2 in self.graph:
            self.graph[seq1].add_edge(self.graph[seq2], 'out')
            self.graph[seq2].add_edge(self.graph[seq1], 'in')
        elif seq1 in self.graph:
            self.graph[seq1].add_edge(ver2, 'out')
            ver2.add_edge(self.graph[seq1], 'in')
            self.graph[seq2] = ver2
        elif seq2 in self.graph:
            self.graph[seq2].add_edge(ver1, 'in')
            ver1.add_edge(self.graph[seq2], 'out')
            self.graph[seq1] = ver1
        else:
            ver1.add_edge(ver2, 'out')
            ver2.add_edge(ver1, 'in')
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
                self.graph[str(new_edge.vertex_left)].extend_edge(new_edge, Graph.k, "out")
                self.graph[str(new_edge.vertex_right)].extend_edge(new_edge, Graph.k, "in")

        for kmer in to_delete:
            del(self.graph[kmer])

    def save_dot(self, output_file):
        dict = collections.defaultdict(list)
        with open(output_file.name, "w") as outp:
            outp.write("digraph {\n")
            for key, vertex in self.graph.items():
                if vertex.get_edges('in'):
                    for edge in vertex.edges['in']:
                        outp.write("%s -> %s [label = %s];\n" % (vertex.edges['in'][edge].get_vertices("in"), vertex, edge))
                    dict[vertex.edges['in'][edge].get_vertices("in")].append(vertex)
                if vertex.get_edges('out'):
                    for edge in vertex.edges['out']:
                        outp.write("%s -> %s [label = %s];\n" % (vertex, vertex.edges['out'][edge].get_vertices("out"), edge))
                    dict[vertex].append(vertex.edges['out'][edge].get_vertices("out"))
            outp.write("}")
        for key, value in dict.items():
            print(key, value)

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

