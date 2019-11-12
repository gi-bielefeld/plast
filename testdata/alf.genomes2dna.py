#!/usr/bin/env python

from sys import stdout,stderr,exit,argv
from optparse import OptionParser
from os.path import basename
import logging

from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import re

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

FASTA_HEADER_PAT = re.compile('^(G\d+_SE\d+), sequence type: (.*), locus: (-?\d+)$')

DEFAULT_MIN_LENGTH = 100
DEFAULT_CHROMOSOME_NO = 6

if __name__ == '__main__':

    usage = 'usage: %prog [options] <{ALF SIMULATED GENOME}_dna.fa>'
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)

    
    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    
    cf = logging.FileHandler('%s.log' %(basename(argv[0]).rsplit('.', 1)[0]), mode='w')
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #

    
    genome = list()

    for rec in SeqIO.parse(open(args[0]), 'fasta'):
        m = FASTA_HEADER_PAT.match(rec.description)
        if not m:
            print >> stderr, 'Unable to parse FASTA sequence header "%s". Exiting' %(rec.description)
            exit(1)
        locus = int(m.group(3))
        if locus < 0:
            locus = -locus
            genome.append((locus, str(rec.seq.reverse_complement())))
        else:
            genome.append((locus, str(rec.seq)))

    genome.sort()
    _, dna = zip(*genome)

    gName = basename(args[0]).rsplit('_dna.fa', 1)[0]
    gRec = SeqRecord(Seq(''.join(dna), generic_dna), 
                id='%s|chromosome|1' %gName, description='')
    SeqIO.write(gRec, stdout, 'fasta')

