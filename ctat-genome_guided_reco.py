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

    ## cmd-line argument processing
    arg_parser = argparse.ArgumentParser(
        
        description = "runs STAR & Stringtie",
        formatter_class=argparse.RawTextHelpFormatter
        )

    arg_parser._action_groups.pop()
    required = arg_parser.add_argument_group('required arguments')
    optional = arg_parser.add_argument_group('optional arguments')
    
    required.add_argument("--left_fq", dest="left_fq", type=str, required=True, default=None,
                          help="left (or single) fastq file")
    
    required.add_argument("--right_fq", dest="right_fq", type=str, required=False, default="",
                          help="right fastq file (optional)")
    
    optional.add_argument("--genome_lib_dir", dest="genome_lib_dir", type=str, default=os.environ.get('CTAT_GENOME_LIB'),
                          help="genome lib directory - see http://FusionFilter.github.io for details.  Uses env var CTAT_GENOME_LIB as default")
    
    optional.add_argument("--out_dir", dest="out_dir", type=str, required=False, default="star_stringtie", help="output directory")
    
    optional.add_argument("--out_prefix", dest="out_prefix", type=str, default='sample', help="output filename prefix (ie. sample name")

    optional.add_argument("--CPU", dest="CPU", required=False, type=int, default=4,
                          help="number of threads for multithreaded processes")

    optional.add_argument("--stringtie_incl_ref_annot", action='store_true', default=False, help="include reference annotation as guide for stringtie")
    
    args_parsed = arg_parser.parse_args()

    ## done arg parsing

    
    genome_lib_dir = args_parsed.genome_lib_dir
    if genome_lib_dir is None:
        raise RuntimeError("Error, must specify --genome_lib_dir or set env var CTAT_GENOME_LIB")
    
    genome_lib_dir = os.path.abspath(genome_lib_dir)
    
            
    args_parsed.left_fq = os.path.abspath(args_parsed.left_fq)
    check_files_exist( [ args_parsed.left_fq ])
                                                   
    if args_parsed.right_fq:
        args_parsed.right_fq = os.path.abspath(args_parsed.right_fq)
        check_files_exist( [ args_parsed.right_fq ])


    if not os.path.exists(args_parsed.out_dir):
        os.mkdir(args_parsed.out_dir)

    os.chdir(args_parsed.out_dir)  ## Now in the workdir for all ops below

    checkpoints_dir = "__stringtie_star_chckpts"
    ## Construct pipeline
    pipeliner = Pipeliner(checkpoints_dir)
    
    ## run STAR
            
    star_cmd = str("STAR --genomeDir {}/ref_genome.fa.star.idx".format(genome_lib_dir) +
                   " --outReadsUnmapped None " +
                   " --runThreadN {} ".format(args_parsed.CPU) +
                   " --outSAMstrandField intronMotif " +
                   " --limitBAMsortRAM {}".format(40 * 1024**3) +  # use 40 G
                   " --outSAMtype BAM SortedByCoordinate " +
                   " --genomeLoad NoSharedMemory " +
                   " --twopassMode Basic " +
                   " --readFilesIn {} {}".format(args_parsed.left_fq, args_parsed.right_fq) )
    
    if re.search(".gz$", args_parsed.left_fq) is not None: 
        star_cmd = star_cmd + " --readFilesCommand 'gunzip -c' "
        
    
    pipeliner.add_commands([Command(star_cmd, "star_align.ok")])

    # rename according to the sample name

    target_bam_filename = "{}.Aligned.sortedByCoord.out.bam".format(args_parsed.out_prefix)
    cmd = "mv Aligned.sortedByCoord.out.bam {}".format(target_bam_filename)

    pipeliner.add_commands([Command(cmd, "renamed_bam.ok")])

    
    ## run Stringtie

    stringtie_gtf = "{}.stringtie.gtf".format(args_parsed.out_prefix)
    stringtie_cmd = "stringtie {} -o {}".format(target_bam_filename, stringtie_gtf)
    stringtie_chkpt = "stringtie.ok"
    if args_parsed.stringtie_incl_ref_annot:
        ref_annot_gtf_filename = os.path.sep.join([genome_lib_dir, "ref_annot.gtf"])
        stringtie_gtf = "{}.stringtie.incRefAnnot.gtf".format(args_parsed.out_prefix)
        stringtie_cmd = "stringtie {} -G {} -o {}".format(target_bam_filename, ref_annot_gtf_filename, stringtie_gtf)
        stringtie_chkpt = "stringtie.incRefAnnot.ok"

    pipeliner.add_commands([Command(stringtie_cmd, stringtie_chkpt)])
    
    
    
    pipeliner.run()
    
    exit(0)




def check_files_exist(files_lst):

    missing = False
    for file in files_lst:
        if not os.path.exists(file):
            print("Error, cannot locate file: {}\n".format(file), file=sys.stderr)
            missing = True

    if missing:
        raise RuntimeError("Error, missing files as indicated")




if __name__ == '__main__':

    main()








