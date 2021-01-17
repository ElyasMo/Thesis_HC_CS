# Thesis_HC_CS
The aim of this study is to find new genes and drugs which directly or indirectly play a role in Heparan sulfate and Chondroitin sulfate expansion at the surface of cancer cells

## Parsing the data
[Slinky R package](http://bioconductor.org/packages/slinky) was used to retreive the level three of LINCS L1000 gene expression data. Highest dose (10µm) and time points (24h) were considered.

```R
library(slinky)

setwd("Directory")
getwd()
key <- "personal_key"
gctx <- "Directory/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx"
info <- "Directory/GSE70138_Broad_LINCS_inst_info.txt"
sl <- Slinky(key, gctx, info)
col.ix <- which(metadata(sl)$cell_id =="HT29" & metadata(sl)$pert_iname == "DMSO" & metadata(sl)$pert_time == "24")
data <- readGCTX(sl[, col.ix])
write.table(data,"Directory/Name.txt", sep="\t")
```
The control (Cancer cell which is treated with DMSO) and treated (Cancer cell which is treated with various perturbages) data was retreived for MCF7 (human breast adenocarcinoma cell line; ATCC HTB-22), A549 (human non-small cell lung carcinoma cell line; ATCC CCL-185), HePG2 (human hepatocellular carcinoma cell line; ATCC HB-8065), HT29 (human colorectal adenocarcinoma cell line; ATCC HTB-38) cancer cell lines.
 
## Data manipulation
To comput the LogFoldChange (LFC) Between control and treated data, [NumPy library](https://numpy.org/) was used.

Example for computing LFC for one cancer cell line:
```Python
import pandas as pd
import numpy as np
import scipy
import math
import openpyxl
from openpyxl import Workbook
from scipy import stats

%cd Directory
%pwd

# Importing control and treated data
Treated=pd.read_csv("HEPG2_treated.csv")
Control=pd.read_csv("HEPG2_control.csv")

#Specifying the average value of gene expressions for replicates in the same batch in control group. In LINCS L1000 naming system, LJ00(number) specify the group of perturbagens which 
would be DMSO for all samples in control group and x(number) detrmines the batch.
LJP005_Average_x1=Control.loc[:, "LJP005_average_x1"]
LJP005_Average_x2=Control.loc[:, "LJP005_average_x2"]
LJP005_Average_x3=Control.loc[:, "LJP005_average_x3"]

LJP006_Average_x1=Control.loc[:, "LJP006_average_x1"]
LJP006_Average_x2=Control.loc[:, "LJP006_average_x2"]
LJP006_Average_x3=Control.loc[:, "LJP006_average_x3"]

LJP007_Average_x1=Control.loc[:, "LJP007_average_x1"]
LJP007_Average_x2=Control.loc[:, "LJP007_average_x2"]
LJP007_Average_x3=Control.loc[:, "LJP007_average_x3"]

LJP008_Average_x1=Control.loc[:, "LJP008_average_x1"]
LJP008_Average_x2=Control.loc[:, "LJP008_average_x2"]
LJP008_Average_x3=Control.loc[:, "LJP008_average_x3"]

LJP009_Average_x1=Control.loc[:, "LJP009_average_x1"]
LJP009_Average_x2=Control.loc[:, "LJP009_average_x2"]
LJP009_Average_x3=Control.loc[:, "LJP009_average_x3"]

#Categorizing the batchs in treated samples
LJP005_x1=Treated.iloc[:, np.r_[1:60]]
LJP005_x2=Treated.iloc[:, np.r_[60:117]]
LJP005_x3=Treated.iloc[:, np.r_[117:173]]

LJP006_x1=Treated.iloc[:, np.r_[173:231]]
LJP006_x2=Treated.iloc[:, np.r_[231:286]]
LJP006_x3=Treated.iloc[:, np.r_[286:344]]

LJP007_x1=Treated.iloc[:, np.r_[344:407]]
LJP007_x2=Treated.iloc[:, np.r_[407:469]]
LJP007_x3=Treated.iloc[:, np.r_[469:529]]

LJP008_x1=Treated.iloc[:, np.r_[529:590]]
LJP008_x2=Treated.iloc[:, np.r_[590:652]]
LJP008_x3=Treated.iloc[:, np.r_[652:713]]

LJP009_x1=Treated.iloc[:, np.r_[713:775]]
LJP009_x2=Treated.iloc[:, np.r_[775:835]]
LJP009_x3=Treated.iloc[:, np.r_[835:895]]

#retriving column names of treated samples forfurther steps
LJP005_x1_list=list(LJP005_x1.columns)
LJP005_x2_list=list(LJP005_x2.columns)
LJP005_x3_list=list(LJP005_x3.columns)

LJP006_x1_list=list(LJP006_x1.columns)
LJP006_x2_list=list(LJP006_x2.columns)
LJP006_x3_list=list(LJP006_x3.columns)

LJP007_x1_list=list(LJP007_x1.columns)
LJP007_x2_list=list(LJP007_x2.columns)
LJP007_x3_list=list(LJP007_x3.columns)

LJP008_x1_list=list(LJP008_x1.columns)
LJP008_x2_list=list(LJP008_x2.columns)
LJP008_x3_list=list(LJP008_x3.columns)

LJP009_x1_list=list(LJP009_x1.columns)
LJP009_x2_list=list(LJP009_x2.columns)
LJP009_x3_list=list(LJP009_x3.columns)

#LFC Calculation for x1-3 plates in LJP005
#this is an example for one batch and this computation should be done for all batches indivudally. 
m=0
m=0
LFC_LJP005_x1=pd.DataFrame()
LFC_LJP005_x2=pd.DataFrame()
LFC_LJP005_x3=pd.DataFrame()

for char in LJP005_x1_list:
    if m < len(LJP005_x1_list):
        d=LJP005_x1_list[m]
        LFC_x1=(np.log2(((LJP005_x1[d])+1)/((LJP005_Average_x1)+1)))
        LFC_LJP005_x1[d]=LFC_x1
        m=m+1
LFC_LJP005_x1.to_excel("LFC_LJP005_x1.xlsx")
```
Accordingly, we will have four matrixes for four cancer cell lines. In rows we can see the gene symbols and in columns we can see the information about purterbagens (LJ00), time point (24h), batches (x), and even wells in each plates as abbteviations (Figure 1).

![LFC](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/LFC_Example.PNG)
                                   **Figure 1: LFC matrix for a cancer cell line**
## Gene-Gene correlation
Various methods could be used to compute gene-gene correlation. In this study, we compared Spearman's rank correlation coefficient, Pearson correlation coefficient, and Kendall rank correlation coefficient. Based on the performance of these methods, one of them was considered as the gene-gene correlation reference method.
The [SciPy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) in python3 was used to calculate gene-gene correlations.

```Python
#an example of calculating gen-gene correlation for one cancer cell line
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import kendalltau

pro=pd.read_csv('A549_LFC_total_genenames.csv', index_col=0)
pro1=pd.read_csv('A549_LFC_total_genenames_filtered.csv')

#To correct the excel file regarding the changing gene names to dates
pro=pro.rename(index={"Sep-2":"Sptin-2","Mar-6":"MACHF6","Mar-7":"MACHF7","Sep-8":"SPTIN8","Mar-2":"MACHF2","Sep-4":"SPTIN4","Sep-10":"SPTIN10","Sep-7":"SPTIN7","Mar-3":"MACHF3","Sep-6":"SPTIN6","Mar-5":"MACHF5","Mar-1":"MTAC1","Dec-1":"ELEC1","Mar-2":"MTRC2","Mar-8":"MACHF8","Sep-9":"SPTIN9"})
#To skip NaNs
pro=pro.dropna()

symbol=pro1['gene symbol']
#preparingthe gene-gene symbol pairs

xInds = []
yInds = []
for i in range(len(pro.index)):
    for j in range(i+1, len(pro.index)):
        xInds.append(symbol[i])
        yInds.append(symbol[j]) 
symbol={'g1':xInds, 'g2':yInds}
symbol=pd.DataFrame(symbol)
$preparing the gene-gene indexes for the correlatio method
xInds = []
yInds = []
for i in range(len(pro.index)):
    for j in range(i+1, len(pro.index)):
        xInds.append(i)
        yInds.append(j) 
#the gene-gene correlation computation        
z=0
Rlist_sp = []
Plist_sp = []
Rlist_pe = []
Plist_pe = []
while z < len(xInds):
    b=xInds[z]
    c=yInds[z]    
    spR, spP = spearmanr(pro.iloc[b].values, pro.iloc[c].values)
    peR, peP = pearsonr(pro.iloc[b].values, pro.iloc[c].values) 
    keR, keP = kendalltau(pro.iloc[b].values, pro.iloc[c].values) 
    Rlist_sp.append(spR)
    Plist_sp.append(spP)
    Rlist_pe.append(peR)
    Plist_pe.append(peP) 
    Rlist_ke.append(keR)
    Plist_ke.append(keP)
    z=z+1    
G1=pd.DataFrame(xInds)
G2=pd.DataFrame(yInds)
R_sp=pd.DataFrame(Rlist_sp)
P_sp=pd.DataFrame(Plist_sp)
R_pe=pd.DataFrame(Rlist_pe)
P_pe=pd.DataFrame(Plist_pe) 
R_ke=pd.DataFrame(Rlist_ke)
P_ke=pd.DataFrame(Plist_ke)
Final=pd.concat([G1, G2, R_sp, R_pe, R_ke P_sp, P_pe, P_ke], axis=1)
np.savetxt('out_gg_A549.txt', Final.values, fmt='%s', delimiter='\t') 

#False discovery rate computation
#an example for calculating FDR based on one Pvalue (peP). The same procedure will be followed for other Pvalues
df_fdr=pd.DataFrame()
x=0
p_vals=Final['peP']
from scipy.stats import rankdata
ranked_p_values = rankdata(p_vals)
fdr = p_vals * len(p_vals) / ranked_p_values
fdr[fdr > 1] = 1
df_fdr=pd.DataFrame(fdr)
.
.
.
df_fdr= pd.concat([fdr_pe, fdr_sp, fdr_ke], axis=1, join='inner')
pe = pd.concat([symbol,Final, df_fdr], axis=1, join='inner')
```
## Functional analysis to decide which statistical method is the best.
According to [Kumari et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050411), the first 100 and 500 gene pairs (based on the lowest FDR) were chosen for functional analysis (both GO and KEGG pathway analysis). The aim was to determine which statistical method can extract more meaningful correlations. To address this, we used [Clustprofiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) R package. We investigated the method that can produce more enriched terms based on the first 100 and 500 correlated gene pairs.

```R
library(dbplyr)
library(org.Hs.eg.db)
library(AnnotationHub)
library(DOSE)


setwd("Directory")
DEGs="Directory/A549_500.csv"
DEG = as.data.frame(read.csv(DEGs,sep=',',stringsAsFactors = F,row.names = NULL))
library(clusterProfiler)

#GO analysis
ego3 <- enrichGO(gene         = DEG$top_500paires,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
head(summary(ego3))
dotplot(ego3, x='p.adjust')
barplot(ego3, showCategory=44)
heatplot(ego3, foldChange = NULL)

#Converting the gene symbols to to ENSEMBL ID
library(EnsDb.Hsapiens.v86)

hsens=EnsDb.Hsapiens.v86
my.symbols <- DEG$top_500paires
enterz<- select(hsens,  
                keys = my.symbols, 
                columns = c("ENTREZID", "SYMBOL", "GENEID"), 
                keytype = "SYMBOL")
enterz=enterz[complete.cases(enterz), ]
write.table(enterz,file='new_symbol_enterzid_500.txt',sep = '\t', na = '',row.names = T,col.names=NA)

#KEGG analysis
enterz = as.data.frame(read.csv('new_symbol_enterzid_500.txt',sep='\t',stringsAsFactors = F,row.names = NULL))

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

ego4 <- enrichKEGG(gene         = enterz$ENTREZID,
                   organism = 'hsa',
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.5,
                   qvalueCutoff  = 0.5)
head(summary(ego4))
barplot(ego4)
dotplot(ego4)
```
Functional analysis revlead that Pearson correlation coefficient outperform the other methods. According to **Figure 2**, the number of enriched terms in GO and KEGG analysis of first 100 gene pairs depict the better performance of Kendall rank correlation coefficient method while the number of enriched terms for the first 500 gene pairs showed the advantage of Pearson correlation coefficient and Spearman's rank correlation coefficient was in the third place out of three methods.

![correlation](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/Enriched%20terms.PNG)

**Figure 2. Number of enriched terms in GO and KEGG analysis of first 100 and 500 gene pairs retrived from Spearman's rank correlation coefficient, Pearson correlation coefficient, and Kendall rank correlation coefficient methods.**

Besides, according to the number of enriched terms and their level of significancy (**Figure 3** and **Figure4**), and also number of genes enrolled in enriched terms, Pearson correlation coefficient is the best statistical method to calculate gene-gene correlations in this study.

![Fig3](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/Dotplot.jpg)
**Figure 3. The enriched terms in y-axis and adjusted-pvalue in x-axis shows the higher amount of enriched terms in Pearson correlation coefficient merhod.**

![Fig4](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/Barplot.jpg)
**Figure 4. The enriched terms in y-axis and number of genes which were enrolled in the enriched terms in x-axis shows the advantage Pearson correlation coefficient merhod.**

![Fig5](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/Heatplots.jpg)
**Figure 5. The enriched terms are shown in y-axis and the number of genes which were enrolled in the enriched terms are placed in x-axis. According to the GO analysis, A) is the heatplot of Pearson correlation coefficient for the first five hundered gene-gene corelations. B) is the heatplot of Spearman's rank correlation coefficient for the first five hundered gene-gene corelations. C) is the heatplot of Kendall rank correlation coefficient for the C1) first one hundered gene-gene correlations and c2) first five hundered gene-gene correlation**

## Extracting HS and CS gene informations
After comingto conclusion that Pearson correlation coefficient is the best method to comput gene-gene correlations, this calculation was done for all four cancer cell lines.
Based on the experimentally aproved gene related to HS and CS which were obtained from literture reviewes, all co-expressed genes with these genes were extracted, were filtered based on the FDR<0.05, and were sorted based on the FDR values.

```python
g1=pe.loc[pe['g1'].isin(['SDC1','SDC2','SDC3','SDC4','GPC1','GPC2', 'GPC3', 'GPC4', 'GPC5', 'GPC6', 'PRCAN', 'AGRN',
                              'COL18A1','B3GAT3', 'EXTL2','EXT1','EXT2','NDST1','NDST2','NDST3','NDST4', 'GLCE', 'HS2ST1',
                             'HS6ST1','HS6ST2','HS6ST3', 'HS3ST1','HS3ST2','HS3ST3','HS3ST4','HS3ST5','HS3ST6', 'SULF1',
                             'SULF2','CSGALNACT1','CHSY1','CHPF','CHSY3','CHST11','CHST12','CHS14','CHST3','CHST7','CHS15',
                             'DSE','UST'])]
g2=pe.loc[pe['g2'].isin(['SDC1','SDC2','SDC3','SDC4','GPC1','GPC2', 'GPC3', 'GPC4', 'GPC5', 'GPC6', 'PRCAN', 'AGRN',
                              'COL18A1','B3GAT3', 'EXTL2','EXT1','EXT2','NDST1','NDST2','NDST3','NDST4', 'GLCE', 'HS2ST1',
                             'HS6ST1','HS6ST2','HS6ST3', 'HS3ST1','HS3ST2','HS3ST3','HS3ST4','HS3ST5','HS3ST6', 'SULF1',
                             'SULF2','CSGALNACT1','CHSY1','CHPF','CHSY3','CHST11','CHST12','CHS14','CHST3','CHST7','CHS15',
                             'DSE','UST'])]
frames=[g1,g2]
genes=pd.concat(frames)
genes=genes[genes['fdr_pe']<=0.05]
genes=genes.sort_values(by=["fdr_pe"])
genes.to_csv('genes.csv')
```
Top 200 gene pairs (with lowest FDR) for all available genes which were obtained from the literture was extracted from the "genes" matrix.

```python
list_genes=['SDC1','SDC2','SDC3','SDC4','GPC1','GPC2', 'GPC3', 'GPC4', 'GPC5', 'GPC6', 'PRCAN', 'AGRN',
                              'COL18A1','B3GAT3', 'EXTL2','EXT1','EXT2','NDST1','NDST2','NDST3','NDST4', 'GLCE', 'HS2ST1',
                             'HS6ST1','HS6ST2','HS6ST3', 'HS3ST1','HS3ST2','HS3ST3','HS3ST4','HS3ST5','HS3ST6', 'SULF1',
                             'SULF2','CSGALNACT1','CHSY1','CHPF','CHSY3','CHST11','CHST12','CHS14','CHST3','CHST7','CHS15',
                             'DSE','UST']
s=pd.DataFrame()
w=pd.DataFrame()
v=pd.DataFrame()
genes=pd.read_csv('genes.csv', sep=',', index_col=0)
genes=genes.sort_values(by=["fdr_pe"])
genes.columns=['g1', 'g2', 'peR', 'peP','fdr_pe']
z=0
while z<46:
    x=list_genes[z]
    y=genes.loc[genes['g1'].isin([x])]
    q=genes.loc[genes['g2'].isin([x])]
    y=y.head(100)
    q=q.head(100)
    s=pd.DataFrame(y)
    v=pd.DataFrame(q)
    frames = [w, s, v]
    w=pd.concat(frames)
    z=z+1
w.to_csv('top100_each_genes.csv')
```
Also, top gene-gene correlations for each experimentally aproved gene was extracted seperately in a dataframe.

```python
list_genes=['SDC1','SDC2','SDC3','SDC4','GPC1','GPC2', 'GPC3', 'GPC4', 'GPC5', 'GPC6', 'PRCAN', 'AGRN',
                              'COL18A1','B3GAT3', 'EXTL2','EXT1','EXT2','NDST1','NDST2','NDST3','NDST4', 'GLCE', 'HS2ST1',
                             'HS6ST1','HS6ST2','HS6ST3', 'HS3ST1','HS3ST2','HS3ST3','HS3ST4','HS3ST5','HS3ST6', 'SULF1',
                             'SULF2','CSGALNACT1','CHSY1','CHPF','CHSY3','CHST11','CHST12','CHS14','CHST3','CHST7','CHS15',
                             'DSE','UST']
list_csv=['SDC1.csv','SDC2,csv','SDC3.csv','SDC4.csv','GPC1.csv','GPC2.csv', 'GPC3.csv', 'GPC4.csv', 'GPC5.csv', 'GPC6.csv', 'PRCAN.csv', 'AGRN.csv',
                              'COL18A1.csv','B3GAT3.csv', 'EXTL2.csv','EXT1.csv','EXT2.csv','NDST1.csv','NDST2.csv','NDST3.csv','NDST4.csv', 'GLCE.csv', 'HS2ST1.csv',
                             'HS6ST1.csv','HS6ST2.csv','HS6ST3.csv', 'HS3ST1.csv','HS3ST2.csv','HS3ST3.csv','HS3ST4.csv','HS3ST5.csv','HS3ST6.csv', 'SULF1.csv',
                             'SULF2.csv','CSGALNACT1.csv','CHSY1.csv','CHPF.csv','CHSY3.csv','CHST11.csv','CHST12.csv','CHS14.csv','CHST3.csv','CHST7.csv','CHS15.csv',
                             'DSE.csv','UST.csv']
s=pd.DataFrame()
w=pd.DataFrame()
v=pd.DataFrame()
genes=pd.read_csv('genes.csv', sep=',', index_col=0)
genes=genes.sort_values(by=["fdr_pe"])
genes.columns=['g1', 'g2', 'peR', 'peP','fdr_pe']
z=0
while z<46:
    x=list_genes[z]
    y=genes.loc[genes['g1'].isin([x])]
    q=genes.loc[genes['g2'].isin([x])]
    y=y.head(500)
    q=q.head(500)
    s=pd.DataFrame(y)
    v=pd.DataFrame(q)
    frames = [s, v]
    w=pd.concat(frames)
    w.to_csv(list_csv[z])
    z=z+1
```
In order to discover common co-expressed genes with predefined genes in all four cancer cell lines, it is necessary to have co-expressed genes with each predefined gene for all four cancer cell lines in a dataframe. Accordingly, the 32 out of 46 HS and CS experimentally aproved genes which are available in our dataset was considered for this analysis and 32 dataframes (one for each gene) were prepared which includes four coulumns for significantly co-expressed genes in four cancer cell lines (sorted basedon FDR) and four related FDR columns.

```python
list_csv=['SDC2.csv','SDC3.csv','SDC4.csv','GPC1.csv', 'GPC3.csv', 'GPC4.csv', 'GPC5.csv', 'AGRN.csv',
                              'COL18A1.csv','B3GAT3.csv', 'EXTL2.csv','EXT1.csv','EXT2.csv','NDST1.csv','NDST2.csv','NDST3.csv','NDST4.csv', 'GLCE.csv', 'HS2ST1.csv',
                             'HS6ST1.csv', 'HS3ST1.csv','HS3ST2.csv', 'SULF1.csv',
                             'CSGALNACT1.csv','CHSY1.csv','CHPF.csv','CHST11.csv','CHST12.csv','CHST3.csv','CHST7.csv',
                             'DSE.csv','UST.csv']
z=0
while z<32:
    %cd "D:\P.H.D\Thesis\new\Matrixes\MCF7"
    x=list_csv[z]
    MCF7=pd.read_csv(x,  usecols=range(1,6))
    MCF7=MCF7[['g1','g2','fdr_pe']]
    MCF7=MCF7.sort_values(by=["fdr_pe"])
    MCF7_1=MCF7[['g1','fdr_pe']]
    MCF7_1.columns=['g_MCF7','fdr_pe']
    MCF7_2=MCF7[['g2', 'fdr_pe']]
    MCF7_2.columns=['g_MCF7','fdr_pe']
    frame=MCF7_1.append(MCF7_2, ignore_index=True)
    frame=frame.sort_values(by=["fdr_pe"])
    frame=frame.drop_duplicates(subset='g_MCF7', keep="first")
    frame_MCF7=frame.sort_index(ignore_index=True)
    %cd "D:\P.H.D\Thesis\new\Matrixes\HT29"
    HT29=pd.read_csv(x,  usecols=range(1,6))
    HT29=HT29[['g1','g2','fdr_pe']]
    HT29=HT29.sort_values(by=["fdr_pe"])
    HT29_1=HT29[['g1','fdr_pe']]
    HT29_1.columns=['g_HT29','fdr_pe']
    HT29_2=HT29[['g2', 'fdr_pe']]
    HT29_2.columns=['g_HT29','fdr_pe']
    frame=HT29_1.append(HT29_2, ignore_index=True)
    frame=frame.sort_values(by=["fdr_pe"])
    frame=frame.drop_duplicates(subset='g_HT29', keep="first")
    frame_HT29=frame.sort_index(ignore_index=True)
    %cd "D:\P.H.D\Thesis\new\Matrixes\A549" 
    A549=pd.read_csv(x,  usecols=range(1,6))
    A549=A549[['g1','g2','fdr_pe']]
    A549=A549.sort_values(by=["fdr_pe"])
    A549_1=A549[['g1','fdr_pe']]
    A549_1.columns=['g_A549','fdr_pe']
    A549_2=A549[['g2', 'fdr_pe']]
    A549_2.columns=['g_A549','fdr_pe']
    frame=A549_1.append(A549_2, ignore_index=True)
    frame=frame.sort_values(by=["fdr_pe"])
    frame=frame.drop_duplicates(subset='g_A549', keep="first")
    frame_A549=frame.sort_index(ignore_index=True)
    %cd "D:\P.H.D\Thesis\new\Matrixes\HEPG2" 
    HEPG2=pd.read_csv(x,  usecols=range(1,6))
    HEPG2=HEPG2[['g1','g2','fdr_pe']]
    HEPG2=HEPG2.sort_values(by=["fdr_pe"])
    HEPG2_1=HEPG2[['g1','fdr_pe']]
    HEPG2_1.columns=['g_HEPG2','fdr_pe']
    HEPG2_2=HEPG2[['g2', 'fdr_pe']]
    HEPG2_2.columns=['g_HEPG2','fdr_pe']
    frame=HEPG2_1.append(HEPG2_2, ignore_index=True)
    frame=frame.sort_values(by=["fdr_pe"])
    frame=frame.drop_duplicates(subset='g_HEPG2', keep="first")
    frame_HEPG2=frame.sort_index(ignore_index=True)
    w=pd.concat([frame_MCF7,frame_HT29,frame_A549,frame_HEPG2], axis=1, join='inner')
    w.columns=['g_MCF7','fdr_MCF7','g_HT29','fdr_HT29','g_A549','fdr_A549','g_HEPG2','fdr_HEPG2']
    %cd "D:\P.H.D\Thesis\new\Matrixes\Merged\Pairs\NEW" 
    w.to_csv(list_csv[z])
    z=z+1
```
To investigate the common co-expressed genes for each HS and CS genes in all four cancer cell lines a [Venn diagram visualisation tool](https://bioinfogp.cnb.csic.es/tools/venny/) was used.
