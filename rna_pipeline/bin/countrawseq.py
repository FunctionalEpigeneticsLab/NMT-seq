#!/usr/bin/env python

from __future__ import division
import sys
import os
import gzip
import subprocess

def check_dirslash(d):
    if d.strip()[-1] == '/':
        return(d)
    else:
        d = d + '/'
        return(d)
    print(d)

def collect_raw_sequence_count(samplelist, fastqdir, outfh):
    '''
    samplelist - one sample id per line
    fastqdir - directory of fastq files
    outfh - output file for sequence count
    '''
    fastqdir = check_dirslash(fastqdir)

    with open(outfh, 'w') as fo:
        with open(samplelist, 'r') as fh:
            for line in fh:
                sampleid = line.strip()
                samR1count = 0
                samR2count = 0
                for fastqfile in os.listdir(fastqdir):
                    if fastqfile.startswith(sampleid):
                        if fastqfile.endswith('R1.fastq.gz'):
                            R1countcmd = f'zcat {fastqdir}{fastqfile} | wc -l'
                            cursamR1count = subprocess.run(R1countcmd, shell=True,stdout=subprocess.PIPE)
                            cursamR1count = cursamR1count.stdout.decode()
                            samR1count = samR1count + int(cursamR1count)/4
                        elif fastqfile.endswith('R2.fastq.gz'):
                            R2countcmd = f'zcat {fastqdir}{fastqfile} | wc -l'
                            cursamR2count = subprocess.run(R2countcmd, shell=True,stdout=subprocess.PIPE)
                            cursamR2count = cursamR2count.stdout.decode()
                            samR2count = samR2count + int(cursamR2count)/4

                pfline = f'{sampleid}\t{samR1count}\t{samR2count}\n'
                fo.write(pfline)

if __name__ == '__main__':
    collect_raw_sequence_count(sys.argv[1], sys.argv[2], sys.argv[3])
