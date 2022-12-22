# Investigating the role of CBX2 in breast cancer through analysis of differentially expressed genes

### By Sojeong Gwon, Fang Shi, Joseph Luper Tsenum, Cheng Zhang, Siming Zhao 

### Genomics and Applied Bioinformatics Project

### Contact Email: jtsenum3@gatech.edu

Engineered Biosystems Building (EBB)

Georgia Institute of Technology

950 Atlantic Drive Atlanta, GA 30332

### Introduction 

Proteins that regulate the gene expression process can be potential therapeutic targets for cancer. CBX2 is one of the important components of Polycomb Repressive Complex 1 (PRC1), which is involved in histone modification and chromatin remodeling (Jangal, Lebeau & Witcher, 2019). CBX2 expression is elevated in triple-negative breast cancer (TNBC), and its knockdown in TNBC models reduced cell numbers by inhibiting proliferation (Bilton et al., 2022). Similar effects have been observed in estrogen receptor (ER)-positive breast cancer (Bilton et al., 2022). Our chosen dataset is composed of sequenced mRNA profiles of two different types of breast cancer cells (MDA-MB-231 for TNBC and MCF-7 for ER-positive) where CBX2 has been knocked down with siRNA. Our project aims to analyze the effect of CBX2 knockdown on gene expression in MDA-MB-231 and MCF-7 cells. MCF-7 and MDA-MB-231 are both common breast cancer cell lines for ER-positive and TNBC, respectively. MCF-7 is estrogen- and EGF- (Epidermal Growth Factor) dependent, while MDA-MB-231 is hormone-independent. By analyzing the difference in the transcriptomic profiles of these two cell lines with CBX2 knockdown, we can further identify potential CBX2-regulated genes across two different breast cancer types.  


### Methods 

Datasets 

The processed data can be downloaded here: https://usegalaxy.org/u/joseph_luper_tsenum/h/all-combined-with-grch38

The raw RNA-Seq data can be found here: https://usegalaxy.org/u/joseph_luper_tsenum/h/raw-rna-seq-data

a. CBX2 Knockdown in MCF-7 cells 

GSE198418 

Platform: Illumina NovaSeq 6000 (Homo sapiens) 

Control: siSCR (non-silencing scrambled control siRNA) 

CBX2 Knockdown: siCBX2#3 (CBX2 targeting siRNA)

![now1](https://user-images.githubusercontent.com/58364462/208567266-6bbb8b71-6656-4031-b6a6-ea2236055d72.png)


b. CBX2 Knockdown in MDA-MB-231 cells 

GSE198417 

Platform: Illumina NovaSeq 6000 (Homo sapiens) 

Control: siSCR (non-silencing scrambled control siRNA) 

CBX2 Knockdown: siCBX2#3 (CBX2 targeting siRNA) 

![now2](https://user-images.githubusercontent.com/58364462/208567463-a4b0bfd2-4175-4317-9b22-d3e32530ada9.png)


### Overall Workflow

                                   Overall workflow of RNA-seq data analysis

<img width="539" alt="image" src="https://user-images.githubusercontent.com/58364462/208562768-7a9603dc-40d2-41de-a8f9-2f3e6b76de0b.png">

Kallisto Quant and DESeq2 were performed on Galaxy (https://usegalaxy.org), a graphical user interface providing various data analysis tools (Afgan et al., 2018).

### Quantifying Expression

Kallisto Quant is an alignment-free (pseudoalignment) expression estimation tool that quantifies the abundancies of RNA-seq transcripts (Bray, Pimentel, Melsted & Pachter, 2016). We first imported our datasets to Galaxy using Galaxy’s wrapper tool called “Faster Download and Extract Reads in FASTQ”, a tool specifically made for retrieving FASTQ files off SRA. Then, we performed Kallisto Quant in the paired-end format to generate TPMs (transcript per million) for each sample. All other settings were run at the default setting and no sequence bias correction or bootstrapping was performed. The reference file used was the prebuilt index file for Homo Sapiens provided by the Pachter lab at Caltech, which developed Kallisto. The Kallisto conversion was run in batches, and each file's tabular output was used for downstream processing. Kallisto version 0.46.2 was used in this analysis.

## Differential Gene Expression (DEG) Analysis and Principal Component Analysis (PCA)

DESeq2 (v1.34.0) is a differential gene expression analysis tool based on the negative binomial distribution (Love, Huber & Anders, 2014). The quantified results of each sample generated by Kallisto Quant in the tabular format were imported as input files. Tximport is the default method to input the data. The .gtf file from the Ensembl homo sapiens transcriptome shared by the Pachter lab was used as the annotation file. DESeq2 was performed to identify genes that are differentially expressed with CBX2 Knockdown relative to the controls on MCF-7 and MDA-MB-231 cell lines, respectively. All other settings were run by default.

Benjamini-Hochberg procedure was used for false discovery rate (FDR) control for our differential expression gene analysis. Genes significantly differentially expressed were defined as genes with an adjusted p-value of below 0.05 and a log2 fold-change of greater than 0.58, representing a minimum of a 1.5-fold change. We also performed principal components analysis (PCA) using covariance on MCF-7 datasets based on our differential expressed genes identified.

### Data import and summarization

We used Tximport to import TPMs, estimated counts, and transcript lengths generated from Kallisto Quant (Soneson, Love, & Robinson, 2015). Tximport also converted transcript ID into Gene ID using a table of transcript-to-gene data provided by the Kallisto developers, available on their Github (https://github.com/pachterlab/kallisto-transcriptome-indices/releases). All other Tximport settings were kept at defaults. The outputs were consolidated into a matrix table of both counts and abundances for all the genes covered by the reference genome used in Kallisto.  Only the abundance data were used for downstream analysis. Non-DEGs, according to DESeq2, were filtered out of the dataset. The R script for Tximport is provided in Supplementary Methods.

### Hierarchical Clustering

We performed agglomerative hierarchical clustering using the default complete-linkage on our differential expressed genes from MCF-7 datasets. First, data were scaled using R’s scale() function, which converts the original data in a column of a dataframe into z-scores, such that the column has a mean of 0 and a standard deviation of 1. This scaling function is a prerequisite for clustering using R’s hierarchical clustering function. 

Next, a Euclidean distance matrix was generated using dist(),  and the distance matrix was used as an input for hclust(), which clustered the data using the complete-linkage agglomerative method, where each element begins as an individual cluster and then is sequentially clustered until the entire dataset is consolidated. Clustering was performed along both the sample and gene axes. A heatmap of differential expression genes was also generated using the heatmap.2 function from R’s gplots package. R scripts used in this analysis is provided in Supplementary Methods.

### Pathway analysis

An online pathway analysis tool, Reactome, was used to perform pathway analysis (Gillespie et al., 2022). The dataset from Tximport (a matrix table of abundances for each sample) was supplied to Reactome, and the pathway analysis with the down-weighting of overlapping genes (PADOG) method was used to perform the differential expression analysis between two groups by a default setting. The control groups for the MCF-7 cell type were compared to the CBX2 knockdown treatment groups.

### Results

Differentially Expressed Gene (DEG) analysis using Kallisto and DESeq2

Using DESeq2 on the Galaxy platform, differentially expressed genes (DEGs) were identified for both cell lines. A summary of the total number of DEGs for each cell line is shown in the figure below, separated by up- and down-regulated genes.

                      Total numbers of differentially expressed genes for MCF-7 and MD-MBA-231
                                 
![image](https://user-images.githubusercontent.com/58364462/208563613-49a861e0-71f5-49fc-bc60-526cd83e07cc.png)

For the MCF-7 cell line, there were a total of 5024 differentially expressed genes, with 2560 up-regulated genes in the treatment cell line and 2464 down-regulated genes. For the MD-MBA-231 cell line, there were no differentially expressed genes identified.
   
A list of the top 10 up- and down-regulated DEGs for the MCF7 cell lines treated with siCBX2-knockdown is listed below:

                      Top up- and down-regulated genes in MCF-7 cells with CBX2 knockdown
      
![image](https://user-images.githubusercontent.com/58364462/208563743-e39139af-2827-47e5-943e-488b8441fbb8.png)

### Principal Component Analysis (PCA)

Using significantly differentially expressed genes, principal components were determined using the Galaxy DESeq2 tool’s optional visualization function, which generates plots using R. 

In the figure below, the six samples for MCF-7 cell line were plotted using the first two principal components. The replicates for the MCF-7 cells treated with the silent control siRNA are marked pink, while the replicates for the MCF-7 cell lines treated with siCBX2 are marked blue.

                                   PCA result for MCF-7 samples
                                   
![image](https://user-images.githubusercontent.com/58364462/208564084-79b11be5-eecc-4c9d-b6db-729d0ecbf94f.png)

PC1 accounted for 66% of the variance in samples, while PC2 accounted for 29% of the variance in samples. There is a clear separation between the control and treatment groups, with a small batch effect for replicate 3.

In the figure below, the six samples for MDA-MB-231 cell line were plotted using the first two principal components. Again, controls are marked in pink, while samples treated with the siRNA knockdown are marked blue.

                                   PCA result for MD-MBA-231 samples
                                   
![image](https://user-images.githubusercontent.com/58364462/208567905-d5ee25b9-ecc8-453e-80d4-1abeaaa097d1.png)

                                   
The first principal component accounts for 50% of the sample variance, while the second principal component accounts for 41. However, no meaningful difference is observed between control and siCBX2-treated samples, and each replicate of the treatment cell line is most similar to its corresponding control replicate. At this point, the MD-MBA-231 was omitted from further analysis, due to a failure to successfully identify DEGs using the computational methods detailed above.

#### Hierarchical clustering of MCF-7 cell line data by sample and gene

Two-way hierarchical clustering plots were generated using the heatmap.2 function from R’s gplots package. 

                                    Two-way hierarchical clustering plots for MCF-7 samples
                                    
![image](https://user-images.githubusercontent.com/58364462/208564385-40bb6b0f-be66-4b56-9e7c-4a452c655b4a.png)


The hierarchically clustered heatmap above shows the distribution of up- and down-regulated genes across all six samples. The green denotes genes where expression in the associated sample is higher relative to the normalized data, which indicates upregulation in the associated condition relative to its comparison. Red denotes downregulation relative to its comparison. A more intense color indicates a greater degree of up- or down-regulation. 

#### Pathway analysis for MCF-7 cell lines

Pathway analysis visualizations were generated using Reactome’s web application. In the figure below, pathways are grouped by biological function. In the figure below, the size of each node reflects the number of biological agents (including proteins, genes, and other molecules) involved with each pathway, while lines represent a connection between a pathway and its sub-pathways. 

 Reactome pathway analysis visualization. Pathways that are significantly upregulated are marked in yellow, and pathways that are downregulated are marked blue.
                                      
![image](https://user-images.githubusercontent.com/58364462/208564577-9b685ad5-a652-4481-8301-b391e2b043ca.png)

In the figure below, the same information is visualized using Reactome.org’s Reacfoam view. In this view, each “cell” represents one of Reactome’s designated pathways, and each embedded cell represents a sub-pathway.

![image](https://user-images.githubusercontent.com/58364462/208564657-ef22dc65-074f-4655-acd1-7ea46c1c7dfa.png)

Five of the most significant up- and down-regulated pathways are highlighted below for clarity.

                                      Top 5 significant up- and down-regulated pathways from Reactome analysis
                                      
![image](https://user-images.githubusercontent.com/58364462/208564779-ddd69398-935a-4d9d-bf73-e6dfb61d882f.png)

### Discussion

#### Discrepancies in gene expression analysis of the MD-MBA-231 cell line

While our workflow successfully generated a robust list of DEGs for the MCF-7 cell line using the 3 replicates per treatment condition, it identified no DEGs for the MD-MBA-231 cell line. This result was supported by a PCA analysis of the results, which showed that the majority of the variance in the dataset was between replicates, rather than between experimental conditions. The expression profiles of each treated sample were most similar to its corresponding replicate of the control sample.

This result differs from the originally published paper, which was successful in identifying up-regulated genes and pathways in MD-MBA-231 cells treated with siCBX2, compared to the controls. The analysis pipeline used by the authors of the original paper differed from our pipeline in several ways:
    
• For alignment, Bilton et al. used HiSAT2, whereas our pipeline used Kallisto. Previous studies (Liu et al., 2022) have indicated that Kallisto may lack sensitivity for lower-expression genes that would be detected using HiSAT2. Kallisto’s pseudoalignment method, while reducing computational needs, was demonstrated to sometimes produce a different number of DEGs in comparison with HiSAT2, which uses the burrows-wheeler transform to align reads to the reference. While our pipeline used a pre-built reference genome off of hg38 by Kallisto’s developers, HiSAT2 allows for a wider range of reference genome input formats and allows users to build their own reference genomes that may be more specific to the samples being aligned.
    
• For quantifying gene expression, Bilton et al. used fragments per kilobase of transcript per million fragments mapped (FPKM) whereas Kallisto measures abundance in transcripts per million (TPM). Again, prior studies have shown that TPM and FPKM, which use different approaches to normalize for transcript length, reads per sample and sequencing biases, do not necessarily produce comparable results in downstream analysis like hierarchical even when applied to the same dataset (Dillies et al., 2013; Zhao et al., 2021).
    
In the mTORC1 genes highlighted in the paper, we did not observe the same statistically significant changes in expression between control and treatment groups. Two genes highlighted in the paper were TSC1 and PRKAA2. In our DESeq2 results, TSC1 had a log-fold change of only 0.067758, with an adjusted p-value of 0.99. PRKAA2 had log-fold change of -0.05702, with an adjusted p-value of 0.99.

Due to the lack of any differentially expressed genes in the MD-MBA-231 dataset, we are not able to draw any biological conclusions about the effect of CBX2 knockdown on MD-MBA-231 cells or compare those results to the results of the MCF-7 cell line. However, our results underscore the importance of rigor and transparency in RNA-seq methodology. The use of different tools for the various steps of gene expression analysis can yield results that lead to significantly different conclusions about the biology and mechanism of action of a particular treatment. This suggests that it may be important to attempt to replicate RNA-seq analysis using additional tools or provide a clearer justification for the use of certain tools over others.

### Up-regulated pathways

Bilton et al. (2022) have reported that CBX2 knockdown in MDA‐MB‐231 cells upregulates expression of the mTORC1 inhibitors TSC1 and PRKAA2. Mammalian target of rapamycin (mTOR) is known to regulate cell proliferation, autophagy, and apoptosis, and studies have shown that mTORC1 is often activated in tumors and is predominantly associated with cell growth and metabolism (Zou et al., 2020 ). TSC1 forms a protein complex that inhibits Rheb, which is a crucial activator of mTORC1 signaling (Manning & Huang, 2008).  Bilton et al. (2022) observed similar effects of CBX2 knockdown in MCF-7 cells as well. In our current study, we corroborated the previous finding of upregulated TSC1 (but not PRKAA2) in MCF-7 cells that have CBX2 knocked down.

Another notable gene found in our study is TP53BP1, which is significantly upregulated in CBX2 knockdown cells compared to control. TP53BP1 participates in DNA repair pathway and maintains genomic stability (Mirza-Aghazadeh-Attari et al., 2019). Previous studies (Ward et al., 2005; Difilippantonio et al. 2008; Morales et al., 2006) suggest that TP53BP1 may act as a critical tumor suppressor because of its role in DNA repair.

### Down-regulated pathways

Reactome analysis identified some significant down-regulated pathways that may be significantly involved in breast cancer progression. For example, the most down-regulated pathway is R-HSA-1483101, the synthesis of phosphatidylserine. The presence of phosphatidylserine on the surface of cells is an engulfment signal that encourages phagocytosis in apoptotic cells (Yu et al., 2020).  Studies have suggested that PS exposure is significantly upregulated on the surface of tumor cells, and drugs that target PS have been developed and are being tested as part of drug cocktails to target various types of cancers (Chang et al., 2020).

Similarly, the second most down-regulated pathway is the macroautophagy pathway (R-HSA-1632852). Macroautophagy is the process by which portions of a cell’s cytoplasm are packaged and delivered to the lysosome for degradation. In normal cells, it recycles organelles and molecules into hydrolytic enzymes, providing new sources of energy while also removing intracellular “garbage”. In cancer cells, mutations along the autophagy pathway can have manifold effects, but are particularly impactful in breast cancer stem cells (BCSC). Wolf et al. (2013) demonstrated using RNAi screens that genes involved in regulating macroautophagy are critical for maintaining pluripotency in breast cancer-derived stem cells and also improved the tumorigenic potential of breast cancer cells in tissue graft experiments. Mechanism studies have suggested that this may have been due to the involvement of macroautophagy in the IL6/STAT3 and TGFB/SMAD pathways (Niklaus et al, 2021).

A gene of significant interest in breast cancer research that displayed statistically significant differential gene expression in our dataset is GAPDH. GAPDH expression has been shown to increase breast cancer cell proliferation and tumor aggressiveness (Guo et al. 2013). Prior studies have shown GAPDH is increased in a wide range of cell cancer types, but most relevant to our study, it has been demonstrated to be significantly upregulated in prior experiments with MCF7 cell lines (Revillion et al. 2000). In our dataset, GAPDH has a statistically significant logFC of -2.441. Although the mechanism of GAPDH’s association with cancer proliferation is not yet known, researchers have hypothesized that it may be involved in cell cycle regulation, because levels of GAPDH gene expression and protein quantity fluctuate in different stages of the cell cycle.

Another gene that has been associated with breast cancer proliferation that was down-regulated in our MCF-7 data set was ABCE1, a gene that has previously been shown to be overexpressed in breast cancer tissues. ABCE1 is an RNase L inhibitor, an interferon that targets and destroys intracellular DNA. It is involved in the 2-5A/RNase L pathway, which inhibits cell apoptosis via the modulation of cell metabolism. Huang et al. (2014) used siRNA to knock down ABCE1 in MCF-7 cells, which resulted in a decrease in cell viability and health, and a change in cell morphology, which they hypothesized was due to the increase in RNAse L expression as a consequence of the ABCE1 knockdown. In our data analysis, the statistically significant logFC of the ABCE1 gene in the treatment cell line was -0.6539.

he downregulation of these pathways and genes, all of which have been shown to be associated with greater progression of breast cancer when expressed at high levels, demonstrates the involvement of CBX2 in the expression of various genes that are directly involved in cancer cell proliferation, metastasis and development. The results of our analysis affirm the claim of the original authors that CBX2 is a potential target for breast cancer therapeutics, due to its influence on various significant biological pathways in breast cancers across multiple types.


References

1. Afgan, E., Baker, D., Batut, B., Van Den Beek, M., Bouvier, D., Čech, M., ... & Blankenberg, D. (2018). The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. Nucleic acids research, 46(W1), W537-W544.

2. Bilton, L. J., Warren, C., Humphries, R. M., Kalsi, S., Waters, E., Francis, T., ... & Wade, M. A. (2022). The Epigenetic Regulatory Protein CBX2 Promotes mTORC1 Signalling and Inhibits DREAM Complex Activity to Drive Breast Cancer Cell Growth. Cancers, 14(14), 3491

3. Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525-527.
Chang, W., Fa, H., Xiao, D., & Wang, J. (2020). Targeting phosphatidylserine for Cancer therapy: prospects and challenges. Theranostics, 10(20), 9214–9229. https://doi.org/10.7150/thno.45125

4. Difilippantonio, S., Gapud, E., Wong, N., Huang, C. Y., Mahowald, G., Chen, H. T., Kruhlak, M. J., Callen, E., Livak, F., Nussenzweig, M. C., Sleckman, B. P., & Nussenzweig, A. (2008). 53BP1 facilitates long-range DNA end-joining during V(D)J recombination. Nature, 456(7221), 529–533. https://doi.org/10.1038/nature07476

5. Dillies, M. A., Rau, A., Aubert, J., Hennequet-Antier, C., Jeanmougin, M., Servant, N., Keime, C., Marot, G., Castel, D., Estelle, J., Guernec, G., Jagla, B., Jouneau, L., Laloë, D., Le Gall, C., Schaëffer, B., Le Crom, S., Guedj, M., Jaffrézic, F., & French StatOmique Consortium (2013). A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis. Briefings in bioinformatics, 14(6), 671–683. https://doi.org/10.1093/bib/bbs046

6. Gillespie, M., Jassal, B., Stephan, R., Milacic, M., Rothfels, K., Senff-Ribeiro, A., ... & D’Eustachio, P. (2022). The reactome pathway knowledgebase 2022. Nucleic acids research, 50(D1), D687-D692.

7. Guo, C., Liu, S., & Sun, M. Z. (2013). Novel insight into the role of GAPDH playing in tumor. Clinical & translational oncology : official publication of the Federation of Spanish Oncology Societies and of the National Cancer Institute of Mexico, 15(3), 167–172. https://doi.org/10.1007/s12094-012-0924-x
Huang, B., Zhou, H., Lang, X., & Liu, Z. (2014). siRNA‑induced ABCE1 silencing inhibits proliferation and invasion of breast cancer cells. Molecular medicine reports, 10(4), 1685–1690. https://doi.org/10.3892/mmr.2014.2424

8. Huang, J., & Manning, B. D. (2008). The TSC1–TSC2 complex: a molecular switchboard controlling cell growth. Biochemical Journal, 412(2), 179-190.
Jangal, M., Lebeau, B., & Witcher, M. (2019). Beyond EZH2: is the polycomb protein CBX2 an emerging target for anti-cancer therapy? Expert opinion on therapeutic targets, 23(7), 565-578. 

9. Liu, X., Zhao, J., Xue, L., Zhao, T., Ding, W., Han, Y., & Ye, H. (2022). A comparison of transcriptome analysis methods with reference genome. BMC genomics, 23(1), 1-15.

10. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 1-21.

11. Mirza-Aghazadeh-Attari, M., Mohammadzadeh, A., Yousefi, B., Mihanfar, A., Karimian, A., & Majidinia, M. (2019). 53BP1: A key player of DNA damage response with critical functions in cancer. DNA repair, 73, 110–119. https://doi.org/10.1016/j.dnarep.2018.11.008

12. Morales, J. C., Franco, S., Murphy, M. M., Bassing, C. H., Mills, K. D., Adams, M. M., ... & Carpenter, P. B. (2006). 53BP1 and p53 synergize to suppress genomic instability and lymphomagenesis. Proceedings of the National Academy of Sciences, 103(9), 3310-3315.

13. Niklaus, N. J., Tokarchuk, I., Zbinden, M., Schläfli, A. M., Maycotte, P., & Tschan, M. P. (2021). The Multifaceted Functions of Autophagy in Breast Cancer Development and Treatment. Cells, 10(6), 1447. https://doi.org/10.3390/cells10061447

14. Revillion F, Pawlowski V, Hornez L et al (2000) Glyceraldehyde-3-phosphate dehydrogenase gene expression in human breast cancer. Eur J Cancer 36:1038–1042

15. Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4.

16. Ward, I. M., Difilippantonio, S., Minn, K., Mueller, M. D., Molina, J. R., Yu, X., Frisk, C. S., Ried, T., Nussenzweig, A., & Chen, J. (2005). 53BP1 cooperates with p53 and functions as a haploinsufficient tumor suppressor in mice. Molecular and cellular biology, 25(22), 10079–10086. https://doi.org/10.1128/MCB.25.22.10079-10086.2005

17. Wolf, J., Dewi, D. L., Fredebohm, J., Müller-Decker, K., Flechtenmacher, C., Hoheisel, J. D., & Boettcher, M. (2013). A mammosphere formation RNAi screen reveals that ATG4A promotes a breast cancer stem-like phenotype. Breast cancer research : BCR, 15(6), R109. https://doi.org/10.1186/bcr3576

18. Yu, M., Li, T., Li, B., Liu, Y., Wang, L., Zhang, J., ... & Shi, J. (2020). Phosphatidylserine-exposing blood cells, microparticles and neutrophil extracellular traps increase procoagulant activity in patients with pancreatic cancer. Thrombosis Research, 188, 5-16.

19. Zhao, Y., Li, M. C., Konaté, M. M., Chen, L., Das, B., Karlovich, C., ... & McShane, L. M. (2021). TPM, FPKM, or normalized counts? A comparative study of quantification measures for the analysis of RNA-seq data from the NCI patient-derived models repository. Journal of translational medicine, 19(1), 1-15.

20. Zou, Z., Tao, T., Li, H., & Zhu, X. (2020). mTOR signaling pathway and mTOR inhibitors in cancer: Progress and challenges. Cell & Bioscience, 10(1), 1-11.


### Supplementary Methods


tximport.R

<img width="800" alt="image" src="https://user-images.githubusercontent.com/58364462/208565751-545be4dc-ca9c-4f23-8e87-2c21cfe9b080.png">


clustering.R

![image](https://user-images.githubusercontent.com/58364462/208565950-a32717bf-fb02-4245-9095-213f95ab125c.png)


                                 


