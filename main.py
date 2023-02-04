# Code to read in two different proteomics results files from two different runs of DIA based MS proteomics data
# and selectively map only peptides for proteins where the peptides detected are identical at the peptide sequence
# level to enable an apples to apples comparison.

#!/usr/bin/python3
import sys
import argparse
import os
import re
import pandas as pd
import multiprocessing as mp
from multiprocessing import Pool
from numpy import mean
from sortedcontainers import SortedList

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
                        required=True,
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
                        required=True,
                        help='sample annotation file to add to R ExpressionSet Object')
    parser.add_argument('-o',
                        action='store',
                        help='path to folder to store output results')
    return parser

def main(args=None):

    parser = get_parser()
    args = parser.parse_args(args)

    output_path = "./"
    if args.o:
        output_path = args.o
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            os.makedirs(output_path + "/Dataset_1")
            os.makedirs(output_path + "/Dataset_2")
            os.makedirs(output_path + "/Intersection")
    else:
        os.makedirs(output_path + "/Dataset_1")
        os.makedirs(output_path + "/Dataset_2")
        os.makedirs(output_path + "/Intersection")

    d1_dict = {}
    d2_dict = {}
    # Get the unique peptide sequences detected between the two datasets, 
    # then compare to get the intersection for our comparitive analysis
    print('Step 1:')
    print('Getting unique set of peptides detected by MS platform for dataset 1.')
    get_unique_peptide_sequences(args.d1, args.d1_type, d1_dict)
    print('Getting unique set of peptides detected by MS platform for dataset 2.')
    get_unique_peptide_sequences(args.d2, args.d2_type, d2_dict)

    d1_len = len(d1_dict.keys())
    d2_len = len(d2_dict.keys())
    print('Dataset 1: ' + str(d1_len) + ', Dataset 2: ' + str(d2_len))

    print('Step 2:')
    print('Identify intersection of unique peptide sequences between two datasets.')
    intersection = d1_dict.keys() & d2_dict.keys()
    intersection_list = list(intersection)
    ordered_intersection_list = SortedList(intersection_list)
    num_of_matches = len(intersection)
    print('There are ' + str(num_of_matches) + ' matches between the two datasets at the peptide level.')

    # Now that we have significant overlap in peptide sequences, let's create
    # dataframes that store the raw peptide intensities and normalize each
    # dataset using z-score
    print('Step 3:')
    print('Creating raw intensity matrix for dataset 1')
    d1_raw_df, d1_sorted_snames, d1_sorted_pseqs = store_raw_intensities(args.d1, args.d1_type, args.d1_sampleMap)
    print('Writing raw intensity matrix for dataset 1 as csv file')
    d1_raw_fn = output_path + "/Dataset_1/full_raw_" + args.d1_type + ".csv"
    d1_raw_df.to_csv(d1_raw_fn)
    print('Creating raw intensity matrix for dataset 2')
    d2_raw_df, d2_sorted_snames, d2_sorted_pseqs = store_raw_intensities(args.d2, args.d2_type, args.d2_sampleMap)
    print('Writing raw intensity matrix for dataset 1 as csv file')
    d2_raw_fn = output_path + "/Dataset_2/full_raw_" + args.d2_type + ".csv"
    d2_raw_df.to_csv(d2_raw_fn)
 
    print('Step 4:')
    print('Creating zscore matrix for dataset 1')
    d1_zscore_df = d1_raw_df.apply(z_score)
    print('Writing zscore matrix for dataset 1 to csv')
    d1_zscore_fn = output_path + "/Dataset_1/full_zscore_" + args.d1_type  + ".csv"
    d1_zscore_df.to_csv(d1_zscore_fn)
    print('Creating zscore matrix for dataset 2')
    d2_zscore_df = d2_raw_df.apply(z_score)
    print('Writing zscore matrix for dataset 2 to csv')
    d2_zscore_fn = output_path + "/Dataset_2/full_zscore_" + args.d2_type  + ".csv"
    d2_zscore_df.to_csv(d2_zscore_fn)

    # Make full featureData.csv file
    print('Step 5:')
    print('Writing full featureData file for dataset 1: R eSet object creation.')
    create_featureData("full", args.d1_type, output_path, "1", d1_dict, d1_sorted_pseqs)
    print('Writing full featureData file for dataset 2: R eSet object creation.')
    create_featureData("full", args.d2_type, output_path, "2", d2_dict, d2_sorted_pseqs)

    # Make full phenoData.csv file
    print('Step 6:')
    print('Writing full phenoData file for dataset 1: R eSet Object creation.')
    create_full_phenoData(args.d1_type, output_path, "1", d1_sorted_snames, args.d1_phenoData)
    print('Writing full phenoData file for dataset 2: R eSet Object creation.')
    create_full_phenoData(args.d1_type, output_path, "2", d2_sorted_snames, args.d2_phenoData)

    # Subset the raw intensities by the intersected peptide sequences for each dataframe then merge
    print('Step 7:')
    print('Creating intersection raw intensity matrix for dataset 1')
    subset_d1_raw_df = d1_raw_df.loc[ordered_intersection_list]
    print('Creating intersection raw intensity matrix for dataset 2')
    subset_d2_raw_df = d2_raw_df.loc[ordered_intersection_list]
    print('Merging two dataframe by index (by peptide sequence).')
    merged_subset_raw_df = subset_d1_raw_df.join(subset_d2_raw_df)
    print('Writing intersection raw intensity matrix for subset dataset as csv')
    mri_fn = output_path + "/Intersection/intersection_raw_ds1ds2.csv"
    merged_subset_raw_df.to_csv(mri_fn)

    # Subset the z-scores by the intersected peptide_sequences for each dataframe then merge
    print('Step 8:')
    print('Creating intersection zscore intensity matrix for dataset 1')
    subset_d1_zscore_df = d1_zscore_df.loc[ordered_intersection_list]
    print('Creating intersection zscore intensity matrix for dataset 2')
    subset_d2_zscore_df = d2_zscore_df.loc[ordered_intersection_list]
    print('Merging two zscore dataframes by index (by peptide sequence).')
    merged_subset_zscore_df = subset_d1_zscore_df.join(subset_d2_zscore_df)
    print('Writing intersection zscore intensity matrix for dataset 2 as csv')
    mzs_fn = output_path + "/Intersection/intersection_zscore_ds1ds2.csv"
    merged_subset_zscore_df.to_csv(mzs_fn)

    # Make subset featureData.txt file
    print('Step 9:')
    print('Create intersection featureData file for eSet building in R.')
    create_featureData("intersection", args.d2_type, output_path, 2, d2_dict, d2_sorted_pseqs)

    # Make subset phenoData.txt file
    print('Step 10:')
    print('Create intersection phenoData file for eSet building in R.')
    create_merged_phenoData(output_path, d1_sorted_snames, args.d1_phenoData, d2_sorted_snames, args.d2_phenoData)

def create_merged_phenoData(output_path, d1_sorted_snames, d1_phenoData, d2_sorted_snames, d2_phenoData):
    # Read in the d1 input phenoData file and parse.
    fi1 = open(d1_phenoData, 'r')
    df1 = fi1.readlines()
    fi1.close()
    # Read in the d2 input phenoData file and parse.
    fi2 = open(d2_phenoData, 'r')
    df2 = fi2.readlines()
    fi2.close()
    # iterate over d_sorted_snames to match, subset, and reorder if necessary
    pD_dict = {}
    for line in df1[1:]:
        line.strip()
        arr = line.split(",")
        pD_dict[arr[0]] = line

    for line in df2[1:]:
        line.strip()
        arr = line.split(",")
        pD_dict[arr[0]] = line

    fn = output_path + "/Intersection/intersection_phenoData_ds1ds2.csv"
    phenoData = open(fn, 'w')
    pD_header = df1[0].strip()
    phenoData.write(pD_header + "\n")
    for sname in d1_sorted_snames:
        phenoData.write(pD_dict[sname])
    phenoData.write("\n")
    for sname in d2_sorted_snames:
        phenoData.write(pD_dict[sname])
    phenoData.close()
    return

def create_full_phenoData(d_type, output_path, num, d_sorted_snames, phenoData_input):
    
    # Read in the input phenoData file and parse.
    fi = open(phenoData_input, 'r')
    df = fi.readlines()
    fi.close()
    # iterate over d_sorted_snames to match, subset, and reorder if necessary
    pD_dict = {}
    for line in df[1:]:
        line.strip()
        arr = line.split(",")
        pD_dict[arr[0]] = line

    fn = output_path + "/Dataset_" + str(num) + "/full_phenoData_" + d_type + ".csv"
    phenoData = open(fn, 'w')
    pD_header = df[0].strip()
    phenoData.write(pD_header + "\n")
    for sname in d_sorted_snames:
        phenoData.write(pD_dict[sname])
    phenoData.close()

    return

def create_featureData(fD_type, d_type, output_path, num, d_dict, d_sorted_pseqs):
    if fD_type == "full":
        fn = output_path + "/Dataset_" + str(num) + "/" + fD_type + "_featureData_" + d_type + ".csv" 
    elif fD_type == "intersection":
        fn = output_path + "Intersection/" + fD_type + "_featureData_ds1ds2.csv"

    if d_type == 'TMT10':
        featureData = open(fn, 'w')
        fD_header = "Peptide_Sequence,NCBI_Protein_ID\n"
        featureData.write(fD_header)
        for i in d_sorted_pseqs:
            v = d_dict[i]
            line = i + "," + v + "\n"
            featureData.write(line)
        featureData.close()
    elif d_type == 'DIA':
        featureData = open(fn, 'w')
        fD_header = "Peptide_Sequence,UniprotID,Gene_Symbol\n"
        featureData.write(fD_header)
        for i in d_sorted_pseqs:
            v = d_dict[i]
            pid, gs = re.split(',', v)
            line = i + "," + pid + "," + gs + "\n"
            featureData.write(line)
        featureData.close()
    return

def z_score(df): return (df-df.mean())/df.std(ddof=0)

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

    # Format dictionary to pandas DataFrame object
    # key is in the format peptide_sequence + "_" + sample_name and value is an average raw peptide intensity
    # This bit of code is to get a unique list of both peptide sequences and sample names so I can iterate using
    # a nested loop structure to build my matrix file.
    sorted_pseq = SortedList()
    sorted_sname = SortedList()
    for k, v in avg_raw_df.items():
        peptide_seq, sample_name = re.split("_", k)
        if sorted_pseq.__contains__(peptide_seq) == True:
            next
        else:
            sorted_pseq.add(peptide_seq)

        if sorted_sname.__contains__(sample_name) == True:
            next
        else:
            sorted_sname.add(sample_name)
    sorted_sname_as_list = list(sorted_sname.irange())
    sorted_pseq_as_list = list(sorted_pseq.irange())


    cpu_count = mp.cpu_count() - 1
    # Here let's use multiprocessing to improve the matrix creation time
    # Simplistic method to speed up is divide the number of unique pseqs by the number of available
    # cpus and launch that number of processes. Each process would tackle a smaller set of pseqs from 
    # the list to store the matrix in memory, then merge the multiple matrixes together before creating
    # the final matrix in pandas DataFrame format.

    print(str(cpu_count) + 'cpus. How many should we use?')
    total_u_pseqs = len(sorted_pseq_as_list)
    print(str(total_u_pseqs) + 'unique peptides. The numerator!')
    
    x_counter = round(total_u_pseqs / cpu_count)
    print('We need to process an average of ' + str(x_counter) + ' peptides per process.')
    list_of_pseq_indexes = [sorted_pseq_as_list[i * x_counter:(i + 1) * x_counter] 
                            for i in range((len(sorted_pseq_as_list) + x_counter - 1) // x_counter )]

    #print(list_of_pseq_indexes[0])
    # create tuples that contain the parameters for the function build_matrix
    tuples_list = []
    for a in range(len(list_of_pseq_indexes)):
        created_tuple = (a, list_of_pseq_indexes[a], sorted_sname_as_list, avg_raw_df)
        tuples_list.append(created_tuple)
    with Pool() as pool:
        matrices = pool.starmap(build_matrix, tuples_list)
    # starmap forces blocking so that the returns stack in the order they are submitted as defined by the tuples_list
    # just use pandas concat to merge the panda dataframes
    formatted_avg_raw_df = pd.concat(matrices)
    
    return formatted_avg_raw_df, sorted_sname_as_list, sorted_pseq_as_list
    
    
def build_matrix(order, list_of_pseq_indexes, sorted_sname_as_list, avg_raw_df):
    # nested for loops to build the matrix to create the pandas dataframe
    # call a function to do this and then assemble afterwards
    matrix = {}
    for i in list_of_pseq_indexes:
        for j in sorted_sname_as_list:
            built_key = i + "_" + j
            try:
                v = avg_raw_df[built_key]
                if j in matrix:
                    matrix[j].append(v)
                else:
                    matrix[j] = list()
                    matrix[j].append(v)
            except (KeyError):
                if j in matrix:
                    matrix[j].append(float('nan'))
                else:
                    matrix[j] = list()
                    matrix[j].append(float('nan'))
    
    formatted_avg_raw_df = pd.DataFrame(matrix, index=list_of_pseq_indexes)
    print('Finished building matrix ' + str(order))
    return formatted_avg_raw_df
    

def avg_raw_per_peptide_by_sample(raw_df):
    # This function is necesary to average the intensity values per peptide per sample before some
    # form of z-score based normalization can happen between different proteomics datasets. This is
    # especially true in the case of TMT style proteomics where multiple fractions are screened per
    # pool of samples and the same peptide may be identified across multiple fraction data files 
    # for the sample

    avg_raw_df = {}
    for k, v in raw_df.items():
        if len(v) == 1:
            avg_intensity = v[0]
        else:
            avg_intensity = mean(v)

        avg_raw_df[k] = avg_intensity
    
    return avg_raw_df

def store_raw_with_sampleMap_TMT(df, sampleMap_dict, non_avg_raw_df):
    # Because the same peptide will appear multiple times in the input file, we need to store the peptide intensities
    # for the same peptide sequence per sample in a list object to average after all the data has been appended.
    # Averaging happens in a separate function.

    if len(non_avg_raw_df.keys()) > 0:
        sorted_idx = SortedList(non_avg_raw_df.keys())
    else:
        sorted_idx = SortedList()

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

            # Need to speed up search by running binary on a sorted list of peptide sequences.
            # Using SortedList method in SortedContainers library for this! .add function allows
            # for inplace insertion of new values into the list while .__contains__ on the SortedList
            # returns True if found. If not found, then add the peptide sequence to the list.
            if sorted_idx.__contains__(k) == True:
                non_avg_raw_df[k].append(intensity_val)
            else:
                non_avg_raw_df[k] = [intensity_val]
                sorted_idx.add(k)

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
        pro_gs = protein + "," + gene
        dict[peptide_sequence] = pro_gs
    return

def keepAminoAcidsOnly(s):
    t = ''
    for i in s:
        if(i.isalpha()):
            t+=i
    return t

if __name__ == '__main__':
    main()