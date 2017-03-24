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
        self.edge += following_edge[-1]
        self.vertex_right = following_edge.vertex_right
        self.coverage += following_edge.coverage
        return

    def get_vertices(self):
        return self.vertex_left, self.vertex_right

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
            edge = Edge(other, self.vertex)
        if direction == 'out':
            edge = Edge(self.vertex,other)

        if str(edge) not in self.edges[direction]:
            self.edges[direction][str(edge)] = edge
        else:
            self.edges[direction][str(edge)].inc_coverage()

    #def extend_edge(self, edge, k, direction):
        #self.edges[direction][st[:k]]

    def extend_edge(self, edge, k, direction):
        print(self.vertex, edge, direction)
        if direction == 'in':
            key = str(edge)[len(str(edge))-k-1:]
        if direction == 'out':
            key = str(edge)[:k+1]
        print( self.edges[direction])
        if key in self.edges[direction]:
            cov = self.edges[direction][key].coverage # IN
            edge.inc_coverage(cov)
            self.edges[direction].pop(key)
            self.edges[direction][str(edge)] = edge
            print('result', self.edges[direction], self.edges[direction][str(edge)])

    def compress(self):
        #print(self.vertex, self.edges)
        if len(self.edges['in']) != 1 or len(self.edges['out']) != 1:
            return False
        for in_edge in self.edges['in']:
            for out_edge in self.edges['out']:
                self.edges['in'][in_edge].merge(self.edges['out'][out_edge])

        return self.edges['in'][in_edge]
        # Returns False, if cannot be compressed
        # Otherwise compresses this vertex and returns true

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
        if seq1 and seq2 in self.graph:
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

    def split_read(self, read):
        for i in range(len(read) - Graph.k):
            yield read[i:i + Graph.k], read[i + 1:i + 1 + Graph.k]

    def add_seq(self, read):
        for kmer1, kmer2 in self.split_read(read):
            self.add_edge(kmer1, kmer2)

    def compress(self):

        to_delete = []  # List of redundant vertices

        for kmer, vertex in self.graph.items():
            new_edge = vertex.compress()
            if new_edge:
                to_delete.append(kmer)
                self.graph[str(new_edge.vertex_right)].extend_edge(new_edge, Graph.k, 'in')
                self.graph[str(new_edge.vertex_left)].extend_edge(new_edge, Graph.k,'out')

        for kmer in to_delete:
            del(self.graph[kmer])

    # Delete redundant vertex

    def save_dot(self): #, output_file):
        #with open(output_file, "w") as outp:
            #outp.write("digraph {\n")
        for key, vertex in self.graph.items():
            if vertex.get_edges('in'):
                for edge in vertex.edges['in']:
                    print(vertex.edges['in'][edge].get_vertices())
            if vertex.get_edges('out'):
                for edge in vertex.edges['out']:
                    print(vertex.edges['out'][edge].get_vertices())
            #print('%s -> [label="EDGE = %s %s "]\n' % (value,  value.edges['in'], value.edges['out']))
                #outp.write('%s -> [label="EDGE = %s %s "]\n' % (value, value.edges['in'], value.edges['out']))

    def __str__(self):
        return str(self.graph)


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def reverse_complement(seq):
    return ''.join(complement[nt] for nt in seq[::-1])


def read_fastq(f):
    '''
    for line in f:
        name = line.strip()
        seq = next(f).strip()
        next(f)
        next(f)
        yield name, seq
    '''
    for line in f:
        name = line.strip()
        seq = line.strip()
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

    print(graph)
    #for kmer, vertex in graph.graph.items():
      #  print('VERTEX = %s, EDGES = %s %s' % (vertex, vertex.get_edges('in'), vertex.get_edges('out')))



   # for kmer, vertex in graph.graph.items():

        #print('VERTEX = %s, EDGES = %s %s' % (vertex, vertex.get_edges('in'), vertex.get_edges('out')))
    print(graph)

    if args.compress:
        graph.compress()
    #graph.save(args.output)
    graph.save_dot() #args.output)
    print(graph)


if __name__ == '__main__':
    main()
