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