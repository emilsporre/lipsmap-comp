Contents

1. Version information
2. Data preparation
3. Running the script
4. Interpreting results

#### 1. Version information ####
This analysis was run on a linux system (Ubuntu 22.04) using R version 3.6.3 and MSstats version 3.18.5.


#### 2. Data preparation ####
This script is made to analyze LiP-SMap proteomics experiments with data that come in sets of 12 samples (3 conditions * 4 replicates) and has been converted to .mzml format and searched with the EncyclopeDIA search engine. The input data should have the columns "Peptide", "Protein", "numFragments" and one column per sample, named after the date, MS method file and well, e.g. HF_20210113_LSP_DIA_350to1000mz_40min_1to32_10mz_2uL_A1.mzML. The overall structure should look like so:

Peptide				Protein			numFragments	HF_20210113_LSP_DIA_350to1000mz_40min_1to32_10mz_2uL_A1.mzML	HF_20210113_LSP_DIA_350to1000mz_40min_1to32_10mz_2uL_A10.mzML
AAAAELAIPLYR			sp|P77972|ENO_SYNY3	5		82190064							76936552
AAAAGFEK			sp|P09193|PSBC_SYNY3	5		8337106								7604372.5
AAAEAEINDPEGDSPEEDFTGDTELGLK	tr|Q55726|Q55726_SYNY3	5		288031.4							377694.62
AAAEFEVGNLYINR			tr|P74275|P74275_SYNY3	5		649506.3							791845.3
...

In addition to the output from EncyclopeDIA, one needs a layout file that specifices which well belong to what sample. It should contain the columns "well" and "name" where the first column describes the samples position on the plate and the second column describes what metabolite was tested at what concentration and which technical replicate the sample represents. All experiments were run with 4 replicates containing no added metabolite (annotated with 0), 4 replicates containing a low concentration metabolite (annotated with 1) and 4 replicates with high concentration metabolite (annotated with 2). For example, "Suc_0_1" would indicate that the metabolite tested is sucrose, the concentration is zero (blank) and that it is technical replicate number one.

The layout file overall structure should look like so:

well	name
A1	Suc_0_1
A2	Suc_0_2
A3	Suc_0_3
...

A uniprot annotation file is also required, which is the proteome downloaded from UniProt as a .csv file with the following columns:  
"Entry", "Gene names  (ordered locus )", "Gene names  (primary )", "Gene names  (synonym )", "Protein names", "Function [CC]", "Catalytic activity", "Gene ontology (GO)", "Gene ontology (molecular function)", "Gene ontology (biological process)", "Gene ontology (cellular component)", "Subcellular location [CC]", "Pathway", "Cofactor", "Subunit structure [CC]", "EC number", "Length", "Activity regulation", "Kinetics", "Miscellaneous [CC]", "Status", "Protein existence", "Proteomes".


#### 3. Running the script ####
Before running the script the initial part titled "define input" must be updated to reflect the dataset being analyzed. The part of the script in need of updating looks like so:

#############################################################################################
### define input
#############################################################################################

pep_path <- list.files("/home/user/lipsmap-comp/data/pcc6803/raw_data/20210113/", pattern = ".txt")
pep_path <- lapply(pep_path, function(x){
  x <- paste("/home/user/lipsmap-comp/data/pcc6803/raw_data/20210113/", x, sep = "")
})
uniprot_path <- "/home/user/lipsmap-comp/data/pcc6803/annotations/uniprot_pcc6803_20220401.tab"
layout_path <- c("/home/user/lipsmap-comp/data/pcc6803/layout/layout_pcc6803_20210113.csv")
min_pep <- 3
mode <- "group" ### either "group" or "treatm", decides if the samples are grouped by group or treatment
save_dir <- "/home/user/lipsmap-comp/results/"
organism <- "pcc6803"

Here one needs to specify the paths to the EncyclopeDIA output, the uniprot annotation file and the layout file as well as specify what organism the data has been collected from. One can also adjust the minimum number of peptides required to be detected within each condition, the save directory and if the comparison is to be done across metabolite comparisons rather than within. None of these options were adjusted for any of our LiP-SMap experiments.


#### 4. Interpreting the data ####
The script outputs a table with one row per peptide and comparison for each metabolite experiment called pepcomp_organism_date_metabolite. It also outputs a file with all detected peptides prior to the statistical comparison, with the log2 intensities still in place. It also outputs a .png file that is a grid of 4 quality control plots. These can be used to assess the quality of the data and to identify any outliers. In the example data provided, there are a couple of outliers. Finally, a .png file of a volcano plot for each comparison is generated.






