#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   quant.py: genotypes from extracted chromosome 6 reads.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   This file is part of arcasHLA.
#
#   arcasHLA is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   arcasHLA is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with arcasHLA.  If not, see <https://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
import pandas as pd
import json
import pickle
import gzip
import numpy as np
from Bio import pairwise2

from argparse import RawTextHelpFormatter

import os
import sys
import re
import json
import pickle
import argparse
import logging as log
from scipy import stats

from reference import check_ref
from arcas_utilities import *

#-------------------------------------------------------------------------------

__version__     = '0.2'
__date__        = '2019-04-02'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir    = os.path.dirname(os.path.realpath(__file__)) + '/../'
parameters = rootDir + 'dat/info/parameters.p'
hla_fasta  = rootDir + 'dat/info/hla.fasta'
mapping_p  = rootDir + 'dat/info/mapping.p'

#-------------------------------------------------------------------------------
#   Quantification
#-------------------------------------------------------------------------------

def get_eq_class(read):
    def get_eq(read):
        i = 1
        kmer_prev = read[:k]

        while ('N' in kmer_prev or not db[kmer_prev]) and i < len(read) - k:
            kmer_prev = read[i:i+k]
            i += 1

        read_eqs = []

        while i < len(read) - k:
            kmer_curr = read[i:i+k]
            if 'N' in kmer_curr:
                i += 1
                break

            if kmer_curr not in db[kmer_prev]:
                break
            elif len(db[kmer_prev]) == 1 and i < len(read) - k - 2:
                read_eqs.append(eq[kmer_curr])
                while len(db[kmer_prev]) == 1 and i < len(read) - k - 2:
                    kmer_prev = list(db[kmer_prev])[0]
                    i += 1
            else:
                read_eqs.append(eq[kmer_curr])
                kmer_prev = kmer_curr
                i += 1

        if len(read_eqs) == 0: return set()
        return set.intersection(*read_eqs)
    
    return get_eq(read) | get_eq(reverse_complement(read))

def reverse_complement(string):
    complement = ''.maketrans('AGCT','TCGA')
    return string[::-1].translate(complement)

def build_de_bruijn(genotype, all_hlas):
    db = defaultdict(set)
    eq = defaultdict(set)

    for allele,allele_id in genotype.items(): 
        if allele not in all_hlas:
            continue
        seq = str(all_hlas[allele].seq)

        for i in range(len(seq) - k):
            kmer_curr = seq[i:i+k]
            kmer_right = None

            if i < len(seq):
                kmer_right = seq[i+1:i+k+1]

            db[kmer_curr].add(kmer_right)


            eq[kmer_curr].add(allele.split('*')[0])
            
    for allele in all_hlas:
        if allele.split('*')[0] in genotype_genes: continue

        seq = str(all_hlas[allele].seq)

        for i in range(len(seq) - k):
            kmer_curr = seq[i:i+k]
            kmer_right = None

            if i < len(seq):
                kmer_right = seq[i+1:i+k+1]

            db[kmer_curr].add(kmer_right)
            eq[kmer_curr].add(allele.split('*')[0])
    return db, eq

def hash_reads(readIDs, reads1, reads2, db, eq):
    eq_count = defaultdict(int)
    eq_reads = defaultdict(set)
    eq_read_by_gene = defaultdict(set)

    for readID in readIDs:
        read1 = reads1[readID]
        read2 = reads2[readID]

        intersection = get_eq_class(read1) & get_eq_class(read2)
        eq_count[','.join(sorted(intersection))] += 1
        if intersection:
            eq_reads[','.join(sorted(intersection))].add(readID)
            for gene in intersection:
                eq_read_by_gene[gene].add(readID)
                
    return eq_count, eq_reads, eq_read_by_gene

def find_mismatches(a1, a2):
    def get_kmers(idx, seq):
        kmers = set()
        start = max(idx - k, 0)
        stop = min(idx + k, len(seq) - k)
        for i in range(start,stop):
            kmers.add(seq[i:i+k])
        return kmers - shared_kmers
    
    a1 = str(all_hlas[a1].seq)
    a2 = str(all_hlas[a2].seq)

    alignments = pairwise2.align.localms(a1, a2, 2, 0, -10, -1)
    align1, align2, _, _, _ = alignments[0]

    a1_kmers = set()
    for i in range(len(a1)):
        a1_kmers.add(a1[i:i+k])
    a2_kmers = set()
    for i in range(len(a2)):
        a2_kmers.add(a2[i:i+k])

    shared_kmers = a1_kmers & a2_kmers

    kmers_mismatch = defaultdict(set)
    kmers_allele = defaultdict(set)
    mismatch_locations = defaultdict(list)

    n = 0
    gap = 0

    idx1 = 0
    idx2 = 0

    for idx in range(len(align1)):
        if align1[idx] == align2[idx]:
            if gap and idx1_kmers and idx2_kmers:
                mismatch_locations[n].append(idx1)
                mismatch_locations[n].append(idx2)
                
                
                for kmer in kmers1:
                    kmers_mismatch[kmer].add(n)
                    kmers_allele[kmer] |= {1}
                for kmer in kmers2:
                    kmers_mismatch[kmer].add(n)
                    kmers_allele[kmer] |= {2}

                n+=1

            gap = 0
            idx1 += 1
            idx2 += 1
            continue

        if gap == 0:
            kmers1 = set()
            kmers2 = set()

        idx1_kmers = get_kmers(idx1, a1)
        idx2_kmers = get_kmers(idx2, a2)

        if idx1_kmers and idx2_kmers:
            kmers1 |= idx1_kmers
            kmers2 |= idx2_kmers

        if align1[idx] == '-':
            idx2 += 1
        elif align2[idx] == '-':
            idx1 += 1
        else:
            idx1 += 1
            idx2 += 1

        gap += 1
    
    return kmers_mismatch, kmers_allele, n, mismatch_locations

def get_mismatch_coverage(gene, count_both_reads):
    def get_mismatch(read):
        mismatches = set()
        alleles = set()

        for i in range(len(read) - k):
            kmer = read[i:i+k]
            kmer_rc = reverse_complement(kmer)

            mismatches |=  kmers_mismatch[kmer] | kmers_mismatch[kmer_rc]
            alleles |= kmers_allele[kmer] | kmers_allele[kmer_rc]
        return mismatches, alleles
    
    mismatch_coverage = {i:defaultdict(int) for i in range(n)}
    mismatch_genes = {i:[[],[]] for i in range(n)}
    

    for readID in eq_read_by_gene[gene]:
        mismatches1, alleles1 = get_mismatch(reads1[readID])
        mismatches2, alleles2 = get_mismatch(reads2[readID])
        mismatches = mismatches1 | mismatches2
        alleles = alleles1 & alleles2
        
        genes1 = get_eq_class(reads1[readID])
        genes2 = get_eq_class(reads2[readID])

        genes =  genes1 & genes2

        contribution = eq_count[gene]/sum([eq_count[gene_] for gene_ in genes])

        if len(alleles) != 1: continue
            

        allele = alleles.pop()
        
        #for mismatch in mismatches1:
        #    mismatch_genes[mismatch]['allele' + str(allele)].extend(list(genes1))
        #for mismatch in mismatches2:
        #    mismatch_genes[mismatch]['allele' + str(allele)].extend(list(genes2))
        
        for mismatch in mismatches:
            mismatch_genes[mismatch][allele - 1].extend(list(genes))
        
        for mismatch in mismatches:
            if count_both_reads and mismatch in mismatches1 & mismatches2:
                mismatch_coverage[mismatch][allele] += 2*contribution
            else:
                mismatch_coverage[mismatch][allele] += contribution
            
    return mismatch_coverage, mismatch_genes

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    
    all_hlas = SeqIO.to_dict(SeqIO.parse(hla_fasta, 'fasta'))
    with open(mapping_p,'rb') as file:
        allele_mapping_single,allele_mapping_2field,allele_mapping_ggroup = pickle.load(file)
    
    parser = argparse.ArgumentParser(prog='arcasHLA quant',
                                 usage='%(prog)s [options] FASTQs',
                                 add_help=False,
                                 formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('files', 
                        help='list of extracted fastq files', 
                        nargs='*')
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('--subject',
                        help = 'subject name',
                        type = str,
                        default = None)
    
    parser.add_argument('--sample',
                        help = 'sample name',
                        type = str)
    
    parser.add_argument('--genotypes',
                        help = 'tsv file containing HLA genotypes',
                        type = str)
    
    parser.add_argument('-k',
                        help = 'k-mer length (default: 31)',
                        type = int,
                        default = 31)
    
    parser.add_argument('--mapping',
                        help = 'single, 2-field, g-group',
                        type = str,
                        default = 'g-group')

    
    parser.add_argument('--count_both_reads',
                        help = 'if both mates overlap the same mismatch, both reads contribute to read count.',
                        action = 'count')
    
    '''
    
    parser.add_argument('--purity',
                        help = 'purity',
                        type = float,
                        default = np.nan)
    
    parser.add_argument('--ploidy',
                        help = 'ploidy',
                        type = float,
                        default = np.nan)
                        
    '''

    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='out directory\n\n',
                        default='./', 
                        metavar='')
    
    parser.add_argument('--temp', 
                        type=str,
                        help='temp directory\n\n',
                        default='/tmp/', 
                        metavar='')

    args = parser.parse_args()
    k = args.k
    
    # Load genotype
    genotypes = pd.read_csv(args.genotypes, sep = '\t').set_index('subject')
    input_genotype = genotypes.loc[args.subject].to_dict()
    
    # Select mapping based on argument input
    if args.mapping.lower() == 'single':
        allele_mapping = allele_mapping_single
    elif args.mapping.lower() == '2-field':
        allele_mapping = allele_mapping_2field
    else:
        allele_mapping = allele_mapping_ggroup

    # Get all alleles based on mapping
    genotype = dict()
    for allele_id,allele in input_genotype.items():
        if type(allele) != str:
            continue
        for allele in allele_mapping[allele]:
            genotype[allele] = allele_id

    # Which genes to quantify
    genotype_genes = {allele.split('*')[0] for allele in genotype}
    
    # Load reads
    reads1 = dict()
    with gzip.open(args.files[0], 'rt') as file:
        for record in SeqIO.parse(file, 'fastq'):
            reads1[record.id.split('/')[0]] = str(record.seq)

    reads2 = dict()
    with gzip.open(args.files[1], 'rt') as file:
        for record in SeqIO.parse(file, 'fastq'):
            reads2[record.id.split('/')[0]] = str(record.seq)
            
    readIDs = list(reads1.keys())
        
    db, eq = build_de_bruijn(genotype, all_hlas)
                
    eq_count, eq_reads, eq_read_by_gene = hash_reads(readIDs, reads1, reads2, db, eq)
                
        
    results = defaultdict(dict)
    
    genes = [gene for gene in eq_count.keys() if len(gene.split(',')) == 1]
    
    results['gene_count'] = {gene:eq_count[gene] for gene in genes}
    total_count = sum(results['gene_count'].values())
    results['gene_abundance'] = {gene:eq_count[gene]/total_count for gene in genes}
    results['genotyped_genes'] = sorted(genotype_genes)
    
    for gene in sorted(genotype_genes):
        a1 = allele_mapping_single[input_genotype[gene + '1']][0]
        a2 = allele_mapping_single[input_genotype[gene + '2']][0]
        
        results[gene]['allele1'] = a1
        results[gene]['allele2'] = a1
            
        kmers_mismatch, kmers_allele, n, mismatch_locations = find_mismatches(a1, a2)

        results[gene]['mismatch_locations'] = mismatch_locations
        results[gene]['n_mismatches_sites'] = n
                 
        mismatch_coverage, mismatch_genes = get_mismatch_coverage(gene, args.count_both_reads)
        
        #results[gene]['mismatch_genes'] = {i:{k:Counter(l) for k,l in j.items()} for i, j in mismatch_genes.items()}
        results[gene]['mismatch_genes'] = {i:[Counter(k) for k in j] for i, j in mismatch_genes.items()}
        
        cov1 = [x[1] for x in mismatch_coverage.values()]
        cov2 = [x[2] for x in mismatch_coverage.values()]
        
        results[gene]['allele_count'] = [cov1, cov2]
        
        baf1 = [i/(i + j) for i,j in zip(cov1, cov2)]
        baf2 = [j/(i + j) for i,j in zip(cov1, cov2)]
        
        results[gene]['allele_freq'] = [baf1, baf2]
        
    json.dump(results, open(args.outdir + '/' + args.sample + '.json','w'))
#-------------------------------------------------------------------------------