# Code to read in two different proteomics results files from two different runs of DIA based MS proteomics data
# and selectively map only peptides for proteins where the peptides detected are identical at the peptide sequence
# level to enable an apples to apples comparison.

#!/usr/bin/python3
import argparse
import os
import re
import pandas as pd
import time
from sortedcontainers import SortedList
from statistics import mean
from scipy.stats import zscore

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
    parser.add_argument('-d1_sampleMap',
                        action='store',
                        required=True,
                        help='sample mapping file that dereferences sample id to patient id')
    parser.add_argument('-d1_phenoData',
                        action='store',
                        help='sample annotation file to add to R ExpressionSet Object')
    parser.add_argument('-d2',
                        action='store',
                        required=True,
                        help='folder path to second dataset')
    parser.add_argument('-d2_type',
                        action='store',
                        required=True,
                        choices=['TMT10','DIA'],
                        help='choose the type of methodology used to generate peptide level intensity data')
    parser.add_argument('-d2_sampleMap',
                        action='store',
                        required=True,
                        help='sample mapping file that dereferences sample id to patient id')
    parser.add_argument('-d2_phenoData',
                        action='store',
                        help='sample annotation file to add to R ExpressionSet Object')
    return parser

def main(args=None):

    parser = get_parser()
    args = parser.parse_args(args)

    d1_dict = {}
    d2_dict = {}
    # Get the unique peptide sequences detected between the two datasets, 
    # then compare to get the intersection for our comparitive analysis
    get_unique_peptide_sequences(args.d1, args.d1_type, d1_dict)
    get_unique_peptide_sequences(args.d2, args.d2_type, d2_dict)

    d1_len = len(d1_dict.keys())
    d2_len = len(d2_dict.keys())
    print('Dataset 1: ' + str(d1_len) + ', Dataset 2: ' + str(d2_len))
    intersection = d1_dict.keys() & d2_dict.keys()
    intersection_list = list(intersection)
    num_of_matches = len(intersection)
    print('There are ' + str(num_of_matches) + ' matches between the two datasets at the peptide level.')

    # Now that we have significant overlap in peptide sequences, let's create
    # dataframes that store the raw peptide intensities and normalize each
    # dataset using z-score
    print('Creating raw intensity matrix for dataset 1')
    d1_raw_df = store_raw_intensities(args.d1, args.d1_type, args.d1_sampleMap)
    print('Writing raw intensity matrix for dataset 1 as csv file')
    d1_raw_df.to_csv('./full_raw_CPTAC_TMT10.csv')
    print('Creating raw intensity matrix for dataset 2')
    d2_raw_df = store_raw_intensities(args.d2, args.d2_type, args.d2_sampleMap)
    print('Writing raw intensity matrix for dataset 1 as csv file')
    d2_raw_df.to_csv('./full_raw_Earli_DIA.csv')
    print('Creating zscore matrix for dataset 1')
    d1_zscore_df = compute_zscore(d1_raw_df)
    print('Writing zscore matrix for dataset 1 to csv')
    d1_zscore_df.to_csv('./full_zscore_CPTAC_TMT10.csv')
    print('Creating zscore matrix for dataset 2')
    d2_zscore_df = compute_zscore(d2_raw_df)
    print('Writing zscore matrix for dataset 1 to csv')
    d2_zscore_df.to_csv('./full_zscore_Earli_DIA.csv')

    # Subset the raw intensities by the intersected peptide sequences for each dataframe
    print('Creating intersection raw intensity matrix for dataset 1')
    subset_d1_raw_df = d1_raw_df.loc[intersection_list]
    print('Writing intersection raw intensity matrix for dataset 1 as csv')
    subset_d1_raw_df.to_csv('./intersection_raw_CPTAC_TMT10.csv')
    print('Creating intersection raw intensity matrix for dataset 2')
    subset_d2_raw_df = d2_raw_df.loc[intersection_list]
    print('Writing intersection raw intensity matrix for dataset 2 as csv')
    subset_d2_raw_df.to_csv('./intersection_raw_Earli_DIA.csv')

    # Subset the z-scores by the intersected peptide_sequences for each dataframe
    print('Creating intersection zscore intensity matrix for dataset 1')
    subset_d1_zscore_df = d1_zscore_df.loc[intersection_list]
    print('Writing intersection zscore intensity matrix for dataset 1 as csv')
    subset_d1_zscore_df.to_csv('./intersection_zscore_CPTAC_TMT10.csv')
    print('Creating intersection zscore intensity matrix for dataset 2')
    subset_d2_zscore_df = d2_zscore_df.loc[intersection_list]
    print('Writing intersection zscore intensity matrix for dataset 2 as csv')
    subset_d2_zscore_df.to_csv('./intersection_zscore_Earli_DIA.csv')

    # Make full featureData.csv file
    featureData = open('./full_featureData_CPTAC_TMT10.csv', 'w')
    fD_header = "Peptide_Sequence,NCBI_Protein_ID,Source\n"
    featureData.write(fD_header)
    for k,v in d1_dict:
        line = k + "," + v + ",CPTAC\n"
        featureData.write(line)
    featureData.close()

    featureData = open('./full_featureData_Earli_DIA.csv', 'w')
    fD_header = "Peptide_Sequence,UniprotID,Gene_Symbol,Source\n"
    featureData.write(fD_header)
    for k,v in d2_dict:
        pid, gs = re.split(';', v)
        line = k + "," + pid + "," + gs + ",Earli\n"
        featureData.write(line)
    featureData.close()

    # Make subset featureData.txt file
    featureDataIntersection = open('./intersection_featureData.csv', 'w')
    fDI_header = "Peptide_Sequence,UniprotID,Gene_Symbol\n"
    for i in intersection_list:
        v = d2_dict[i]
        pid, gs = re.split(';', v)
        line = i + "," + pid + "," + gs + "\n"
        featureDataIntersection.write(line)
    featureDataIntersection.close()


def compute_zscore(pd_matrix_df):
    pd_zscore_df = pd_matrix_df.apply(zscore)
    return pd_zscore_df

def store_raw_intensities(ds_folder, ds_type, ds_sampleMap):
    sampleMap_dict = {}
    sampleMap_dict = create_sampleMap_dict(ds_type, ds_sampleMap, sampleMap_dict)
    raw_df = {}
    dir = os.fsencode(ds_folder)
    for file in os.listdir(dir):
        filename = os.fsdecode(file)
        if filename.startswith("phenoData"):
            next
        elif filename.endswith(".psm") or filename.endswith(".csv"):
            df = []
            fn = ds_folder + "/" + filename
            fi = open(fn, 'r')
            df = fi.readlines()
            fi.close()
            if ds_type == 'TMT10':
                raw_df = store_raw_with_sampleMap_TMT(df, sampleMap_dict, raw_df)
            elif ds_type == 'DIA':
                raw_df = store_raw_with_sampleMap_DIA(df, sampleMap_dict, raw_df)
    
    avg_raw_df = avg_raw_per_peptide_by_sample(raw_df)
    return avg_raw_df

def avg_raw_per_peptide_by_sample(raw_df):
    # This function is necesary to average the intensity values per peptide per sample before some
    # form of z-score based normalization can happen between different proteomics datasets. This is
    # especially true in the case of TMT style proteomics where multiple fractions are screened per
    # pool of samples and the same peptide may be identified across multiple fraction data files 
    # for the sample

    idx = []
    sorted_idx = SortedList()
    matrix = {}
    for k, v in raw_df.items():
        k_parts = re.split('_', k)
        peptide_sequence = k_parts[0]
        sample_name = k_parts[1]
        if len(v) == 1:
            avg_intensity = v[0]
        else:
            avg_intensity = mean(v)

        # Need to speed up search by running binary on a sorted list of peptide sequences.
        # Using SortedList method in SortedContainers library for this! .add function allows
        # for inplace insertion of new values into the list while .__contains__ on the SortedList
        # returns True if found. If not found, then add the peptide sequence to the list.

        # print('Running binary search')
        # start = time.time()

        if sorted_idx.__contains__(peptide_sequence) == True:
            next
        else:
            idx.append(peptide_sequence)
            sorted_idx.add(peptide_sequence)
            idx_len = len(idx)
            if idx_len % 10000 == 0:
                print('Peptide Sequences: ' + str(idx_len))
        
        # end = time.time()
        # dur = round(end - start, ndigits = 10)
        # print('Binary search ran in {} seconds'.format(dur))

        if sample_name in matrix:
            matrix[sample_name].append(avg_intensity)
        else:
            matrix[sample_name] = [avg_intensity]
    
    # create and return a pandas dataframe object
    avg_raw_df = pd.DataFrame(matrix, index=idx)
    return avg_raw_df

def store_raw_with_sampleMap_TMT(df, sampleMap_dict, non_avg_raw_df):
    # Because the same peptide will appear multiple times in the input file, we need to store the peptide intensities
    # for the same peptide sequence per sample in a list object to average after all the data has been appended.
    # Averaging happens in a separate function.
    for line in df[1:]:
        arr = line.split("\t")
        peptide_sequence = keepAminoAcidsOnly(arr[11])
        fn_parts = re.split('_', arr[0])
        fn_start = fn_parts[0]
        num_samples_by_file_name = len(sampleMap_dict[fn_start])
        for i in range(0, num_samples_by_file_name):
            sample_name = sampleMap_dict[fn_start][i]
            # Kind of specific to the CPTAC TMT10 psm output file in that position 22 is the first column that contains
            # peptide level intensity values. There are up to 10 values recorded per peptide (hence TMT10) as there are
            # ten total pooled samples profiled per fraction run, however, CPTAC sometimes added or did not add an 
            # internal reference control and/or a non patient specific tumor sample not a part of the original specimen collection
            intensity_position = 22 + i
            intensity = arr[intensity_position]
            intensity_parts = re.split('\/', intensity)
            if len(intensity_parts) == 2:
                intensity_val = intensity_parts[0]
            else:
                intensity_val = intensity
                
            k = peptide_sequence + "_" + sample_name
            intensity_val = float(intensity_val)
            if k in non_avg_raw_df:
                non_avg_raw_df[k].append(intensity_val)
            else:
                non_avg_raw_df[k] = [intensity_val]

    return non_avg_raw_df

def store_raw_with_sampleMap_DIA(df, sampleMap_dict, non_avg_raw_df):
    # similar function as the TMT version except the intensities mapping is based on the DIA data format
    # and we remap sample names from biognosys ids to Benchling ids to pull correct phenoData.
    # biognosys generated sample data starts from the 5th column.
    header = df[0].strip()
    header = header.split(",")
    col_num = len(header)
    biognosys_ids = []
    for i in range(5, col_num):
        biognosys_ids.append(header[i])

    for line in df[1:]:
        line = line.strip()
        if re.search(',\".+\",_', line):
            line = re.sub(',\".+\",_', ',,_', line)
        arr = line.split(',')
        peptide_sequence = re.sub('_', '', arr[4])
        if re.search('\\[Acetyl\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[Acetyl\\]', '', peptide_sequence)
        if re.search('\\[CAM\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[CAM\\]', '', peptide_sequence)
        if re.search('\\[OX\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[OX\\]', '', peptide_sequence)
        for j in range(5, col_num):
            biognosys_id = biognosys_ids[j - 5]
            sample_name = sampleMap_dict[biognosys_id]
            k = peptide_sequence + "_" + sample_name
            intensity_val = float(arr[j])
            if k in non_avg_raw_df:
                non_avg_raw_df[k].append(intensity_val)
            else:
                non_avg_raw_df[k] = [intensity_val]
    return non_avg_raw_df

def create_sampleMap_dict(ds_type, ds_sampleMap, sm_dict):
    dat = []
    fn = ds_sampleMap
    fi = open(fn, 'r')
    dat = fi.readlines()
    fi.close()
    if ds_type == 'TMT10':
        for line in dat[1:]:
            line = line.strip()
            arr = line.split("\t")
            # Only grab sample names that are in positions 1-8, otherwise, it refers to internal reference sample
            # Storage for TMT10 data uses the start of the file name containing intensity data as the key and the
            # value is a list of the actual sample names in the order reported depending on the mapping file
            if re.match('130C', arr[5]):
                next
            else: 
                name = arr[1]
                fn = arr[7]
                fn_parts = re.split('_', fn)
                fn_start = str(fn_parts[0])
                if fn_start in sm_dict:
                    sm_dict[fn_start].append(name)
                else:
                    sm_dict[fn_start] = [name]
    elif ds_type == 'DIA':
        for line in dat[1:]:
            line = line.strip()
            arr = line.split("\t")
            benchling = arr[0]
            biognosys = str(arr[1])
            sm_dict[biognosys] = benchling

    return sm_dict

def get_unique_peptide_sequences(ds_folder, ds_type, d_dict):
    if ds_folder and ds_type:
        dir = os.fsencode(ds_folder)
        for file in os.listdir(dir):
            filename = os.fsdecode(file)
            if filename.startswith("phenoData"):
                next
            elif filename.endswith(".psm") or filename.endswith(".csv"):
                print('Reading in data file ' + filename)
                df = []
                fn = ds_folder + "/" + filename
                fi = open(fn, 'r')
                df = fi.readlines()
                fi.close()
                # initially, we want to make a unique list of peptide sequences detected in the dataset.
                if ds_type == 'TMT10':
                    fill_detected_peptides_TMT10(df, d_dict)
                elif ds_type == 'DIA':
                    fill_detected_peptides_DIA(df, d_dict)
    return

def fill_detected_peptides_TMT10(rdata, dict):
    # TMT10 output by CPTAC in tab delimited .psm file format
    # filename and sample name required to deconvolute to patient id and clinical status (tumor or normal) 
    # via lookup with sample mapping file

    # TMT10 peptide sequence strings contain molecular weights for various static modifications due to the TMT reagents 
    # (+229.163 Da) added to lysines and N termini, carbamidomethyl (+57.021 Da) on cysteines and dynamic modifications
    # for oxidation of methionine residues (+15.9949 Da). We can use regex to remove non-amino acid characters from
    # the peptide sequence with confidence

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
    # DIA output by biognosys in comma delimited file
    # DIA peptide sequence strings contain [Acetyl], [CAM] or [OX] for specific positions that need to be removed.

    for line in rdata[1:]:
        line = line.strip()
        if re.search(',\".+\",_', line):
            line = re.sub(',\".+\",_', ',,_', line)
        arr = line.split(',')
        peptide_sequence = re.sub('_', '', arr[4])
        if re.search('\\[Acetyl\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[Acetyl\\]', '', peptide_sequence)
        if re.search('\\[CAM\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[CAM\\]', '', peptide_sequence)
        if re.search('\\[OX\\]', peptide_sequence):
            peptide_sequence = re.sub('\\[OX\\]', '', peptide_sequence)
        protein = arr[0]
        gene = arr[1]
        pro_gs = protein + ";" + gene
        dict[peptide_sequence] = pro_gs
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
