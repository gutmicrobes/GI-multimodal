01 : data_extract.ipynb: Filtering images of gastric and colorectal cancers from within two image sets; 
canny.ipynb: An improved Canny operator for extracting edge features from each sample image.
Input:
1.mrna_diff_tumor.csv: data of stomach and colorectal cancer, we need to use sample IDs.
2. The name of the folder where the two image sets are located: MSIMUT&MSS.
Output:
Save the screened images to the self-built folder; pig_dimna.csv is the output of canny operator edge features.
———————————————————————————————————————————————————————————
02 : In the study, the gene expression data and miRNA expression data of gastrointestinal cancer were obtained from the official website of TCGA, and the files miRNA_sample_sheet.tsv/RNAseq_sample_sheet.tsv were used to integrate the information.
For miRNA expression data, mirbase.rds refers to the standard files that need to be used in the R file to integrate the sample information using combine_TCGA_data.R. Combined_miRNA_Expr.txt is the output of the integrated information, and combined_miRNA_Expr.csv is the txt-converted csv file.
For gene expression information, according to moveFiles.pl extract the zip package of each sample and decompress the folders one by one, according to mRNA_merge.pl we can get mRNAmatrix.txt, use symbol.pl for ID conversion, and use data_candle.R for the fusion of the same probe. The data after probe fusion are unique.symbol.txt, symbol.txt is the file after probe ID conversion.
data_extract.R is the difference analysis code, and the corresponding results can be found in the "Results" folder.
———————————————————————————————————————————————————————————
03 : disTOM.R is the soft threshold selection code.
Input: pig_dimna.csv, finaly_mirna.csv, finaly_mrna.csv
Output: Disstom Infix File (Three Types of Data) - Disstom Matrix Data After Soft Threshold Determination 515*515
Cluster_new_SNF.ipynb is an improved clustering algorithm code.
Input: disstom infix file (three classes of data) ---- disstom matrix data after soft threshold determination (515*515)
Output: new_finally_label.csv.
———————————————————————————————————————————————————————————
04 : WGCNA.R is the WGCNA code.
Input: final_mrna.csv, new_finally_label.csv, sample clustering label.
All three classes correspond to labels first, second, and third, and the gene_spe suffixes are the hub genes obtained from the three classes of data.
KEGG&GO.R is the KEGG&GO code.
Input: The gene id of the hub genes of the three types of samples.
Output: A picture of GO functional analysis.
———————————————————————————————————————————————————————————
The appendice file is a supplement to the article.