# Code to read in two different proteomics results files from two different runs of DIA based MS proteomics data
# and selectively map only peptides for proteins where the peptides detected are identical at the peptide sequence
# level to enable an apples to apples comparison.

#!/usr/bin/python3
import argparse
import os
import re

# First, load in the CPTAC TMT10 data. TMT stands for Tandom Mass Tags, which are used for pooling samples and
# deconvolution of the peptides detected to a sample. CPTAC MS Proteomics was performed on each patient's tumor and
# normal sample, each pooled sample contains 10 samples, 4 tumor and 4 normal and 1 internal reference control sample.
# Each sample has 25 fractions worth of MS detected peptide data generated.

def get_parser():
    parser = argparse.ArgumentParser(description="Get to mapping and merging proteomics data!")
    parser.add_argument('-d1',
                        action='store',
                        required=True,
                        help='folder path to first dataset')
    parser.add_argument('-d1_type',
                        action='store',
                        required=True,
                        choices=['TMT10','DIA'],
                        help='choose the type of methodology used to generate peptide level intensity data')
    parser.add_argument('-d2',
                        action='store',
                        required=True,
                        help='folder path to second dataset')
    parser.add_argument('-d2_type',
                        action='store',
                        required=True,
                        choices=['TMT10','DIA'],
                        help='choose the type of methodology used to generate peptide level intensity data')
    return parser

def main(args=None):

    parser = get_parser()
    args = parser.parse_args(args)

    d1_dict = {}
    d2_dict = {}

    if args.d1 and args.d1_type:
        dir = os.fsencode(args.d1)
        for file in os.listdir(dir):
            filename = os.fsdecode(file)
            if filename.endswith(".psm") or filename.endswith(".csv"):
                print('Reading in data file ' + filename)
                df = []
                fn = args.d1 + "/" + filename
                fi = open(fn, 'r')
                df = fi.readlines()
                fi.close()
                # initially, we want to make a unique list of peptide sequences detected in the dataset.
                if args.d1_type == 'TMT10':
                    fill_detected_peptides_TMT10(df, d1_dict)
                elif args.d1_type == 'DIA':
                    fill_detected_peptides_DIA(df, d1_dict)
    if args.d2 and args.d2_type:
        print('Working on second dataset.')
        dir = os.fsencode(args.d2)
        for file in os.listdir(dir):
            filename = os.fsdecode(file)
            if filename.endswith(".psm") or filename.endswith(".csv"):
                print('Reading in data file ' + filename)
                df = []
                fn = args.d2 + "/" + filename
                fi = open(fn, 'r')
                df = fi.readlines()
                fi.close()
                # initially, we want to make a unique list of peptide sequences detected in the dataset.
                if args.d2_type == 'TMT10':
                    fill_detected_peptides_TMT10(df, d2_dict)
                elif args.d2_type == 'DIA':
                    fill_detected_peptides_DIA(df, d2_dict)

    d1_len = len(d1_dict.keys())
    d2_len = len(d2_dict.keys())
    print('Dataset 1: ' + str(d1_len) + ', Dataset 2: ' + str(d2_len))
    intersection = d1_dict.keys() & d2_dict.keys()
    num_of_matches = len(intersection)
    print('There are ' + str(num_of_matches) + ' matches between the two datasets at the peptide level.')
    test = next(iter(intersection))
    print(test + " TMT10: " + d1_dict[test] + " DIA: " + d2_dict[test])
    


def fill_detected_peptides_TMT10(rdata, dict):
    # TMT10 output by CPTAC in tab delimited .psm file format
    # filename and sample name required to deconvolute to patient id and clinical status (tumor or normal) 
    # via lookup with sample mapping file

    # TMT10 peptide sequence strings contain molecular weights for various static modifications due to the TMT reagents 
    # (+229.163 Da) added to lysines and N termini, carbamidomethyl (+57.021 Da) on cysteines and dynamic modifications
    # for oxidation of methionine residues (+15.9949 Da). We can use regex to remove non-amino acid characters from
    # the peptide sequence with confidence

    h = rdata[0].split('\t')
    s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 = h[23], h[24], h[25], h[26], h[27], h[28], h[29], h[30], h[31], h[32]

    for line in rdata[1:]:
        arr = line.split("\t")
        # Get rid of the modified molecular weight info in the string
        peptide_sequence = keepAminoAcidsOnly(arr[11]) 
        # Get rid of the pre and post amino acid annotation that indicates where the enzymatic digest cut
        protein = arr[13]
        if re.search('\(', protein):
            protein = re.sub(r'\(.+\)', '', protein)
        # print('Protein: ' + protein + '; Peptide Seq: ' + peptide_sequence)
        # Store this in a dictionary where the peptide_sequence detected should be unique
        dict[peptide_sequence] = protein
    return

def fill_detected_peptides_DIA(rdata, dict):
    # DIA output by Biognosis in comma delimited file
    # DIA peptide sequence strings contain [] for specific positions that include multiple amino acids at that position
    # in the peptide. Because of this, multiple peptide sequences should be created to represent the full complexity as
    # possibilities to match against when building the dictionary as there is ambiguity at that position.
    print('Working on a DIA formatted file.')
    for line in rdata[1:]:
        arr = line.split(',')
        peptide_sequence = re.sub('_', '', arr[4])
        peptide_sequence = re.sub('\\[CAM\\]', '', peptide_sequence)
        protein = arr[0]
        dict[peptide_sequence] = protein
        # print('Protein: ' + protein + '; Peptide seq: ' + peptide_sequence)
    return

def keepAminoAcidsOnly(s):
    t = ''
    for i in s:
        if(i.isalpha()):
            t+=i
    return t

if __name__ == '__main__':
    main()



