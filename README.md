# Thesis_HC_CS
The aim of this study is to find new genes and drugs which directly or indirectly play a role in Heparan sulfate and Chondroitin sulfate expansion at the surface of cancer cells

## Parsing the data
[Slinky R package](http://bioconductor.org/packages/slinky) was used to retreive the level three of LINCS L1000 gene expression data. highest dose (10µm) and time points (24h) were considered.

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
The control (Cancer cell which is treated with DMSO) and treated (Cancer cell which is treated with various perturbages) data was retreived for A549 (adenocarcinomic human alveolar basal epithelial
 cells:Lung cancer), MCF7 (Breast cancer), HEPG2 (Liver cancer), and HT-29 (Colon cancer) cancer cell lines.
 
## Data manipulation
To comput the LogFoldChange (LFC) Between control and treated data, [NumPy library](https://numpy.org/) was used.

Example for one cancer cell line:
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
Accordingly, we will have four matrixes for four cancer cell lines. In rows we can see the gene symbols and in columns we can see the information about purterbagens (LJ00), time point (24h), batches (x), and even wells in each plates as abbteviations.

![LFC](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/LFC_Example.PNG)

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
According to [Kumari et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050411), the first 100 and 500 gene pairs (based on the lowest FDR) were chosen for functional analysis. The aim was to determine which statistical method can extract more meaningful correlations. To address this, we used [Clustprofiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) R package. We investigated the method that can produce more enriched terms based on the first 100 and 500 correlated gene pairs.

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
dotplot(ego3, x='p.adjust', showCategory=44)
barplot(ego3, showCategory=44)
heatplot(ego3, showCategory = 44, foldChange = NULL)
