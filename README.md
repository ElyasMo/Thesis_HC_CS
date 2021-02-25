# Pipeline development to propose repositionable drugs which increase the functionality of anti-cancer peptides ##

### In this study we inquire the possibility of escalating anti-cancer peptides (ACPs) functionality by diminishing Heparan Sulfate (HS), Heparan Sulfate ProteoGlycans (HSPGs), and Chondroitin Sulfate (CS) branchs at the surface of cancer cells. The reason is that these cell surface components, act as an obstacle on way of ACPs lytic effect #





![Graphical Abstract](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/abstract.jpg)
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
![Fig 6](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/Example.png)
**Figure 6. An instance on how to find the common co-expressed genes with each HC and CS defined gene in all four cancer cell lines. Accordingly, 2 and 137 common co-expressed genes with AGRN, and B3GAT3 genes in all four cancer cell lines can be seen.**

Next step is to choose top 10 coexpressed genes for each experimentally aproved gene (if available).
In order to visualize the pattern of gene expression and changes against various perturbagens a heatmap was provided for each cancer cell line which have genes as rows and perturbagens as columns. Prior to plotting the heatmap, all gene expressions were sorted in rows to distinguish between up and downregulated expression patterns against various perturbagens.
To do so, first we should retrive the expression profile of the HS and CS genese and their coexpressed genes from the LFC matrixes for all 4 cancer cell lines.

```python
%cd "D:\P.H.D\Thesis\new\Matrixes\main_matrixes"
A549=pd.read_csv('A549_LFC_total.csv')
HT29=pd.read_csv('HT29_LFC_total.csv')
HEPG2=pd.read_csv('HEPG2_LFC_total.csv')
MCF7=pd.read_csv('MCF7_LFC_total.csv')

%cd "D:\P.H.D\Thesis\new\Matrixes\results"
alls=pd.read_csv('all_in_a_column.csv')
all_list=alls['all_top_10'].tolist()
expr_list=alls['expr_apr'].tolist()

A549_expr=A549.loc[A549['gene symbol'].isin(expr_list)]
A549_expr=A549_expr.sort_index(ignore_index=True)
A549_all=A549.loc[A549['gene symbol'].isin(all_list)]
A549_all=A549_all.sort_index(ignore_index=True)

HEPG2_expr=HEPG2.loc[HEPG2['gene symbol'].isin(expr_list)]
HEPG2_expr=HEPG2_expr.sort_index(ignore_index=True)
HEPG2_all=HEPG2.loc[HEPG2['gene symbol'].isin(all_list)]
HEPG2_all=HEPG2_all.sort_index(ignore_index=True)

HT29_expr=HT29.loc[HT29['gene symbol'].isin(expr_list)]
HT29_expr=HT29_expr.sort_index(ignore_index=True)
HT29_all=HT29.loc[HT29['gene symbol'].isin(all_list)]
HT29_all=HT29_all.sort_index(ignore_index=True)

MCF7_expr=MCF7.loc[MCF7['gene symbol'].isin(expr_list)]
MCF7_expr=MCF7_expr.sort_index(ignore_index=True)
MCF7_all=MCF7.loc[MCF7['gene symbol'].isin(all_list)]
MCF7_all=MCF7_all.sort_index(ignore_index=True)

A549_expr.to_csv('A549_expr.csv')
A549_all.to_csv('A549_all.csv')
HEPG2_expr.to_csv('HEPG2_expr.csv')
HEPG2_all.to_csv('HEPG2_all.csv')
HT29_expr.to_csv('HT29_expr.csv')
HT29_all.to_csv('HT29_all.csv')
MCF7_expr.to_csv('MCF7_expr.csv')
MCF7_all.to_csv('MCF7_all.csv')
```
In order to plot the heatmap the gplot package in R was used.

``` R
setwd('Directory')

HS_CS_genes=read.csv("HT29_all_sorted.csv", sep=",", row.names=1) # I import it from the option on up right of the rstudio
matrix=as.matrix(HS_CS_genes)

library(gplots)

yb <-colorRampPalette(c("gold", "black", "blue"))
heatmap.2(matrix, col=yb, trace = "none", margins = c(6,10), cexCol =0.1,cexRow = 0.3, 
          Rowv = FALSE, Colv = FALSE, scale="row", key = TRUE, key.title = "Range"
          ,key.xlab = "LogFoldChange", key.ylab = "Down", keysize = 1, densadj = 0.25, 
          density.info="none", key.par=list(mgp=c(1, 0.5, 0),mar=c(1, 3, 4, 0))) #, key.xtickfun=FALSE
```

Fig 7 represents the effect of all perturbagens on HC and CS genes. Considering that the rows are sorted based on the value of LFC, the perturbagens which cause downregulation are located at the left and the ones which cause upregulation are placed at the right side of the heatmap and accordingly, they can be easily extracted.

![Fig 7](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/all_heatmap%20(1).jpg)
**Figure 7. The heatmap plot of the effect of perturbagens on HC and CS gene expression.**

In order to find the common perturbagens which cause up or down regulations for HS and CS genes, the most effective drugs were extracted for all four cancer cell lines (Fig 8).
![Fig8](https://github.com/ElyasMo/Thesis_HC_CS/blob/main/heatmap1.jpg)
**Figure 8. all important drugs and/or chemichals which cause up or down regulations in HS and CS genes.**

Once again, [Venn diagram visualisation tool](https://bioinfogp.cnb.csic.es/tools/venny/) can be used to find common chemichals which cause the same LFC changes. Accordingly, five common perturbagens caused downregulation and 16 made upregulations in all four cancer cell lines.
