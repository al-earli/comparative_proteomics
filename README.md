# comparative_proteomics
Code to compare peptide sequence fragments between different datasets to use for comparative analytics

This was built because we wanted to be able to compare CPTAC generated proteomics data that used older
Tandem Mass Tagged data (specifically TMT10 - 10 pooled samples, 25 fractions [size based elutions] per pool)
to our internally generated Ultra deep DIA (discovery) proteomics data (using Biognosys).

Because of the differences in intensity values and the potential for issues as the peptides detected between
the two datasets are unlikely to be identical, we wanted to first identify a list of unique peptides based on
there peptide sequences reported that overlap between the two platforms and use that list as a mechanism to
perform an "apples to apples" comparison between the data generated.

Z-score scaled normalization per sample will also be required because the peptide intensity ranges between the
two mass spec based proteomics methodologies are vastly different.

Help to run the program can be executed with the -h flag which results in the following:

comparative_proteomics % ~/opt/anaconda3/bin/python main.py -h
usage: main.py [-h] -d1 D1 -d1_type {TMT10,DIA} -d1_sampleMap D1_SAMPLEMAP -d1_phenoData D1_PHENODATA -d2 D2 -d2_type {TMT10,DIA} -d2_sampleMap
               D2_SAMPLEMAP -d2_phenoData D2_PHENODATA [-o O]

Get to mapping and merging proteomics data!

optional arguments:
  -h, --help            show this help message and exit
  -d1 D1                folder path to first dataset
  -d1_type {TMT10,DIA}  choose the type of methodology used to generate peptide level intensity data
  -d1_sampleMap D1_SAMPLEMAP
                        sample mapping file that dereferences sample id to patient id
  -d1_phenoData D1_PHENODATA
                        sample annotation file to add to R ExpressionSet Object
  -d2 D2                folder path to second dataset
  -d2_type {TMT10,DIA}  choose the type of methodology used to generate peptide level intensity data
  -d2_sampleMap D2_SAMPLEMAP
                        sample mapping file that dereferences sample id to patient id
  -d2_phenoData D2_PHENODATA
                        sample annotation file to add to R ExpressionSet Object
  -o O                  path to folder to store output results


Command used to run and generate results from the example data included here:

~/opt/anaconda3/bin/python main.py 
    -d1 ./example_data_input/TMT_data 
    -d1_type TMT10 
    -d1_sampleMap ./example_data_input/TMT_data/CPTAC3_LUAD_sample_mapping.txt.txt 
    -d1_phenoData ./example_data_input/TMT_data/phenoData_CPTAC_TMT10.csv 
    -d2 ./example_data_input/DIA_data 
    -d2_type DIA 
    -d2_sampleMap ./example_data_input/DIA_data/Earli_to_Biognosis_ids.txt 
    -d2_phenoData ./example_data_input/DIA_data/phenoData_Earli_DIA.csv
    -o ./example_output