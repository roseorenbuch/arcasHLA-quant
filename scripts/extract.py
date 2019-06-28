#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   extract.py: extracts chromosome 6 reads from a BAM file for HLA genotyping.
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

import os
import sys
import re
import pickle
import argparse
import logging as log

from datetime import date
from os.path import isfile
from argparse import RawTextHelpFormatter
from arcas_utilities import *

__version__     = '0.2.0'
__date__        = '2019-06-26'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir       = os.path.dirname(os.path.realpath(__file__)) + '/../'
hla_regions_p = rootDir + 'dat/info/hla_regions.p'
alts_p        = rootDir + 'dat/info/decoys_alts.p'

#-------------------------------------------------------------------------------
#   Extract Reads
#-------------------------------------------------------------------------------

def index_bam(bam):
    '''Attempts to index BAM if .bai file is not found.'''
    if not isfile(''.join([bam, '.bai'])):
        run_command(['samtools', 'index', bam], 
                    '[extract] indexing bam: ')
                    
    if not isfile(''.join([bam, '.bai'])):
        sys.exit('[extract] Error: unable to index bam file.')
        
    
def extract_hlas(bam, sample, chromosome, header, paired, unmapped, alts, hla_regions, outdir, temp, threads):
    '''Extracts reads from HLA regions on chromosome 6 and alts/decoys if applicable'''
    
    log.info(f'[extract] Extracting read IDs from HLA regions.')

    hla_reads = ''.join([temp, sample, '.reads.txt'])
    command = ['touch', hla_reads]
    run_command(command)
    
    hla_filtered = ''.join([temp, sample, '.hla.bam'])
    
    # Get read IDs from each HLA region
    for region in hla_regions:
        command = ['samtools', 'view', '-@'+threads, bam, chromosome+region, '|', 'cut', '-f1', '>>', hla_reads]
        run_command(command)
        
    for alt in alts:
        if alt not in header:
            continue
        message = '[extract] Extracting read IDs from {alt}: '
        command = ['samtools', 'view', '-@'+threads, bam, region, '|', 'cut', '-f1', '>>', hla_reads]
        run_command(command, message)
        
    if unmapped:
        message = '[extract] Extracting unmapped read IDs: '
        command = ['samtools', 'view', '-@'+threads]
        
        if paired: command.append('-f 12')
        else: command.append('-f 4')
        
        command.extend([bam, '|', 'cut', '-f1', '>>', hla_reads])
        run_command(command, message)
        
    command = ['cat', hla_reads, '|', 'sort', '-u', '>>', hla_reads]
    run_command(command)
    
    message = '[extract] Extracting reads by read ID: '
    command = ['picard', 'FilterSamReads', 'RLF=' + hla_reads, 'FILTER=includeReadList', 'I=' + bam, 'O=' + hla_filtered, 'SORT_ORDER=queryname', 'TMP_DIR=' + temp]
    run_command(command, message)
    
    bam_to_fastq(hla_filtered, sample, paired, outdir, threads)

def extract_chr6(bam, sample, chromosome, header, paired, unmapped, alts, outdir, temp, threads):
    '''Extracts reads from chromosome 6 and alts/decoys if applicable.'''

    hla_filtered = ''.join([temp, sample, '.hla.sam'])
    hla_filtered_bam = ''.join([temp, sample, '.hla.bam'])
    
    # Extract BAM header
    command = ['samtools', 'view', '-H', '-@'+threads, bam, '-o', hla_filtered]
    run_command(command)
    
    # Extracted reads mapped to chromosome 6
    message = '[extract] Extracting chromosome 6: '
    command = ['samtools', 'view', '-@'+threads]
    if paired: command.append('-f 2')
    else: command.append('-F 4')
    command.extend([bam, chromosome, '>>', hla_filtered])
    run_command(command, message)
    
    # Extract unmapped reads
    if unmapped:
        message = '[extract] Extracting unmapped reads: '
        command = ['samtools', 'view', '-@'+threads]
        
        if paired: command.append('-f 12')
        else: command.append('-f 4')
        
        command.extend([bam, '>>', hla_filtered])
        run_command(command, message)
    
    # Check for alts in header and extract reads if present
    for alt in alts:
        if alt not in header:
            continue
            
        message = '[extract] Extracting reads from {alt}: '
        command = ['samtools', 'view', '-@'+threads]

        if paired: command.append('-f 2')
        else: command.append('-F 4')

        command.extend([bam, alt+':', '>>', hla_filtered])
        run_command(command, message)

    # Convert SAM to BAM
    message = '[extract] Converting SAM to BAM: '
    command = ['samtools', 'view', '-Sb', '-@'+threads,
                hla_filtered, '>', hla_filtered_bam]    
    run_command(command, message)

    # Sort BAM
    hla_sorted = ''.join([temp, sample, '.hla.sorted.bam'])
    file_list.append(hla_sorted)
    file_list.append(hla_sorted + '.bai')
    message = '[extract] Sorting bam: '
    command = ['samtools', 'sort', '-n', '-@'+threads, 
                hla_filtered_bam, '-o', hla_sorted]
    run_command(command, message)
    
    bam_to_fastq(hla_filtered, sample, paired, outdir, threads)

    
    
def bam_to_fastq(bam, sample, paired, outdir, threads):
    ''' Convert BAM to FASTQ and compress. '''
    message = '[extract] Converting bam to fastq: '
    command = ['bedtools', 'bamtofastq', '-i', bam]
    if paired:
        fq1 = ''.join([outdir, sample, '.extracted.1.fq'])
        fq2 = ''.join([outdir, sample, '.extracted.2.fq'])
        command.extend(['-fq', fq1, '-fq2', fq2])
        run_command(command, message)
        
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq1])
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq2])
        
    else:
        fq = ''.join([outdir, sample, '.extracted.fq'])
        command.extend(['-fq', fq])
        run_command(command, message)
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq])
    
def extract_reads(bam, paired, chr6, unmapped, alts, hla_regions, outdir, temp, threads):
    '''Extracts reads from chromosome 6 and alts/decoys if applicable.'''
    
    log.info(f'[extract] Extracting reads from {bam}')
    
    sample = os.path.splitext(os.path.basename(bam))[0]
    
    # Index bam
    index_bam(bam)
        
    # Get bam header to check for chromosome nomenclature
    output = run_command(['samtools', 'view', '-@'+threads, '-H', bam])
    header = output.stdout.decode('utf-8')
    
    if 'SN:chr' in header: chromosome = 'chr6'
    else: chromosome = '6'
        
    if chr6:
        extract_chr6(bam, sample, chromosome, header, paired, unmapped, alts, outdir, temp, threads)
    else:
        extract_hlas(bam, sample, chromosome, header, paired, unmapped, alts, hla_regions, outdir, temp, threads)


#-------------------------------------------------------------------------------
#   Main
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA extract',
                                     usage='%(prog)s [options] BAM file',
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    
    
    parser.add_argument('bam', 
                        type=str, 
                        help='/path/to/sample.bam')
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
                        
    parser.add_argument('--log', 
                        type=str,
                        help='log file for run summary\n  '+
                        'default: sample.extract.log\n\n',
                        default=None, 
                        metavar='')
     
    parser.add_argument('-r', '--ref',
                        type = str,
                        help='GRCh37 or GRCh38\n\n',
                        metavar='')
    
    parser.add_argument('--paired', 
                        action = 'count',
                        help='paired-end reads\n\n',
                        default=False)
    
    parser.add_argument('--chr6', 
                        action = 'count',
                        help='extract all reads mapped to chromosome 6\n\n',
                        default=False)    
    
    parser.add_argument('--unmapped', 
                        action = 'count',
                        help='include unmapped reads\n\n',
                        default=False)
                        
    parser.add_argument('-o', '--outdir', 
                        type=str,
                        help='out directory\n\n',
                        default='./', metavar='')
                        
    parser.add_argument('--temp', 
                        type=str,
                        help='temp directory\n\n',
                        default='/tmp/', metavar='')
                        
    parser.add_argument('--keep_files',
                        action = 'count',
                        help='keep intermediate files\n\n',
                        default=False)
    
    parser.add_argument('-t',
                        '--threads', 
                        type = str,
                        default='1')
    
    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    outdir = check_path(args.outdir)
    temp = create_temp(args.temp)
    
    sample = os.path.basename(args.bam).split('.')[0]
    
    if args.log:
        log_file = args.log
    else:
        log_file = ''.join([outdir,sample,'.extract.log'])
        
    with open(log_file, 'w'):
        pass
    
    if args.verbose:
        handlers = [log.FileHandler(log_file), log.StreamHandler()]
        
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
    else:
        handlers = [log.FileHandler(log_file)]
            
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
        
    log.info('')
    hline()
    log.info(f'[log] Date: %s', str(date.today()))
    log.info(f'[log] Sample: %s', sample)
    log.info(f'[log] Input file: %s', args.bam)
    log.info('[log] Read type: {}-end'
             .format( 'paired' if args.paired else 'single'))
    hline()
    
    with open(alts_p, 'rb') as file:
        alts = pickle.load(file)
        
    with open(hla_regions_p, 'rb') as file:
        grch37, grch38 = pickle.load(file)
    if args.ref.lower() == 'grch37': 
        hla_regions = grch37
    else: 
        hla_regions = grch38

    extract_reads(args.bam,
                  args.paired,
                  args.chr6,
                  args.unmapped,
                  alts,
                  hla_regions,
                  outdir,
                  temp,
                  args.threads)
    
    remove_files(temp, args.keep_files)
    
    hline()
    log.info('')
#-------------------------------------------------------------------------------