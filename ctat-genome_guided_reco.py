#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os, re
import argparse
import subprocess

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "PyLib"]))
from Pipeliner import Pipeliner, Command

import logging
FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
global logger
logger = logging.getLogger()
logging.basicConfig(filename='FusionInspector.log', format=FORMAT, filemode='w', level=logging.DEBUG)
# add a new Handler to print all INFO and above messages to stdout
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def main():
    
    genome_lib_dir = args_parsed.genome_lib_dir
    if genome_lib_dir is None:
        raise RuntimeError("Error, must specify --genome_lib_dir or set env var CTAT_GENOME_LIB")
    
    genome_lib_dir = os.path.abspath(genome_lib_dir)
    
    args_parsed.gtf_filename = os.path.sep.join([genome_lib_dir, "ref_annot.gtf"])
    args_parsed.genome_fasta_filename = os.path.sep.join([genome_lib_dir, "ref_genome.fa"])
    args_parsed.cdna_fasta_filename = os.path.sep.join([genome_lib_dir, "ref_cdna.fasta"])




if __name__ == '__main__':

    main()








