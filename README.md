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
Various methods could be used to compute gene-gene correlation. In this study, we have used Spearman's rank correlation coefficient, Pearson correlation coefficient, and Kendall rank correlation coefficient. Based on the performance of these methods, one of them was used as the gene-gene correlation method.
the [SciPy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) in python3 was used to calculate gene-gene correlations.

