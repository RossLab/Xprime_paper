#!/usr/bin/env python3

from sys import stderr
import sys

merged_dump_file = sys.argv[1]
kmer_prefix = 'data/kmers_k27'

A_file = kmer_prefix + '_A.fasta'
X_file = kmer_prefix + '_X.fasta'
Xp_file = kmer_prefix + '_Xp.fasta'

A_kmers = 0
X_kmers = 0
Xp_kmers = 0

with open(Xp_file, 'w') as Xp, open(X_file, 'w') as X, open(A_file, 'w') as A, open(merged_dump_file, 'r') as dump:
        for line in dump:
                kmer = line.rstrip().split('\t')
                male = int(kmer[1])
                female = int(kmer[2])
                if 110 < male and male < 190 and 10 < female and female < 40:
                        A_kmers += 1
                        kmer_fasta_record = '>A_k' + str(A_kmers) + '\n' + kmer[0] + '\n'
                        A.write(kmer_fasta_record)
                if 50 < male and male < 100 and 5 < female and female < 35:
                        X_kmers += 1
                        kmer_fasta_record = '>X_k' + str(X_kmers) + '\n' + kmer[0] + '\n'
                        X.write(kmer_fasta_record)
                if male < 5 and 5 < female:
                        Xp_kmers += 1
                        kmer_fasta_record = '>Xp_k' + str(Xp_kmers) + '\n' + kmer[0] + '\n'
                        Xp.write(kmer_fasta_record)
