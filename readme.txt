In this study, the transcriptome data, immune cell data and tissue section image data of gastrointestinal tumor micro were investigated, and the edge features of each tumor patient's tissue section image were extracted with the improved Canny operator, and the differential mRNA, miRNA and immune cell features of cancer patients and normal samples were screened and extracted, Finally, the clustering of gastrointestinal cancer samples was performed, and the hub genes of different cancer subtypes were explored to provide new theoretical methods for the clinical treatment of gastrointestinal tumors.

1, Gastrointestinal Cancer Tissue Image Data Analysis
Code:
data_extract.ipynb is the code to filter the images of gastric and colorectal cancer from inside the two image sets
canny.ipynb is the code to extract the edge features of each sample image by the improved Canny operator
Data:
pig_dimna.csv is the output of the edge features of the canny operator
Input:
1. mrna_diff_tumor.csv data for gastric and colorectal cancers, sample ids are required
2. name of the folder where the two image sets are located: msimut & mss
Output:
Save the screened images to the self-built folder

2. integrate mrna data and mirna data for gastrointestinal cancer data
Input:
combine_TCGA_data.R is the R code for integrating sample information.
mirbase.rds refers to the standard file to be used in the R file.
The data is downloaded from the TCGA official website, and the file miRNA_sample_sheet.tsv is downloaded for integrating the information. 
Output:
combined_miRNA_Expr.txt is the output of combined information, combined_miRNA_Expr.csv is the csv file after txt conversion.

3 Changing ids and grouping.R is the R code for changing gene ids and grouping samples
Firstly, the data corresponding to the two cancers (including patients and normal people) were input, and the mRNA data of the two cancers were converted into the ids of the genes using gencode.v22.annotation.gene.probeMap, and a total of 58,317 genes were retained, which corresponded to the stad_gene_414_58387.csv, crc_ gene_525_58387.csv data. Subsequently, the samples were grouped according to the 14th,15th result of the sample id (greater than 10 for normal samples, otherwise for patient samples).

4 all.R contains several sections, firstly, the data situation of the percentage of immune cells was calculated based on the mRNA data and 12 immune cell categories that differed significantly under both samples were filtered by drawing violin plots. Secondly, the mRNA and miRNA data were analyzed for differences to determine the differential characteristics of up-regulation and down-regulation.
After completing the above sections, intersections were taken for the samples and a total of 515 samples were retained for the four data categories. The names are “nextmrna.csv”, “nextmirna.csv”, “nextimmune.csv”, and “nextimage.csv”, and based on this data, the samples were calculated under the disTOM matrix under each set of features.

5 precluster&finally.py file firstly preclustered the samples under each histology data to get the optimal number of clusters and the corresponding sample labels, then the neighbor matrix and kernel matrix under each histology data were calculated, and finally, the categories to which each sample belonged were obtained by spectral clustering. And the differences in subtypes in different features were analyzed by Differences in subtypes after clustering.R

6: wgcna.R is the WGCNA code. Input mRNA data and sample clustering labels. The final hub genes of each subtype were obtained and analyzed by KEGG&GO using KEGG&GO.R Output: images of GO functional analysis.
