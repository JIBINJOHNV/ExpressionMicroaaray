




#https://www.biostars.org/p/9466751/
#https://www.biostars.org/p/388949/#388960
#https://www.biostars.org/p/281922/#282138
#https://www.biostars.org/p/401319/
#https://www.biostars.org/p/9485414/
#https://bioconductor.org/packages/release/bioc/html/limma.html
#Reference https://www.biostars.org/p/388949/#388960
#




library(biomaRt)
library(limma)


require(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup2 <- getBM(
  mart = mart,
  attributes = c(
    'agilent_gpl6848',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype',
    'external_gene_name'))


write.table(
  annotLookup2,
  paste0('agilent_gpl6848_', gsub("-", "_", as.character(Sys.Date())), '.tsv'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE)


# convert the data to an EListRaw object, which is a data object for single channel data
# specify green.only = TRUE for Agilent
# retain information about background via gIsWellAboveBG


targetinfo <- readTargets("/Users/jibinjohn/Downloads/Microarray_dataanalysis/SG13134300_257236337155_S002_Sample_Details.tsv", sep = '\t')


project <- read.maimages(
                targetinfo, 
                source = 'agilent.median',
                green.only = TRUE,
                other.columns = 'gIsWellAboveBG')

colnames(project) <- gsub('raw\\/', '', colnames(project))


# generate QC plots of raw intensities
dir.create('QC/')
# histograms, box, and density plots - use your own code


##Probe details downloaded from https://www.ebi.ac.uk/biostudies/files/A-GEOD-21185/A-GEOD-21185_comments.txt



# annotate the probes
annotLookup <- read.csv("A-GEOD-21185_comments.txt",header = TRUE,sep = '\t',stringsAsFactors = FALSE)

colnames(annotLookup)[1] <- 'AgilentID'
annotLookup <- annotLookup[which(annotLookup$AgilentID %in% project$genes$ProbeName),]
annotLookup <- annotLookup[match(project$genes$ProbeName, annotLookup$AgilentID),]
table(project$genes$ProbeName == annotLookup$AgilentID) # check that annots are aligned
project$genes$AgilentID <- annotLookup$AgilentID
project$genes$ensembl_gene_id <- annotLookup$ensembl_gene_id
project$genes$entrezgene <- annotLookup$entrezgene
project$genes$gene_biotype <- annotLookup$gene_biotype
project$genes$external_gene_name <- annotLookup$external_gene_name
#project$genes$GeneName <- annotLookup$GeneName


# perform background correction on the fluorescent intensities
project.bgcorrect <- backgroundCorrect(project, method = 'normexp')

# normalize the data with the 'quantile' method
project.bgcorrect.norm <- normalizeBetweenArrays(project.bgcorrect, method = 'quantile')


# filter out control probes, those with no symbol, and those that fail:
Control <- project.bgcorrect.norm$genes$ControlType==1L

NoSymbol <- is.na(project.bgcorrect.norm$genes$GeneName)
IsExpr <- rowSums(project.bgcorrect.norm$other$gIsWellAboveBG>0)>= 2
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & IsExpr, ]


# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'ensembl_gene_id','entrezgene','gene_biotype','GeneName','external_gene_name',"ProbeName","AgilentID")]
head(project.bgcorrect.norm.filt$genes)




expr<-as.data.frame(project.bgcorrect.norm$E)
genes<-as.data.frame(project.bgcorrect.norm$genes)


row.names(expr) <- NULL
row.names(genes) <- NULL
gene_expr<-cbind(genes,expr)


write.csv(gene_expr,"Gene_Expreession.csv", row.names = FALSE)


############################-------------------------PYthon--------------------------------------------------------------

from gprofiler import GProfiler
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns


##https://www.biostars.org/p/209790/

df=pd.read_csv("Gene_Expreession.csv")

#Inveerse logtransorm

for col in [ x for x in df.columns  if "SG13134300" in x]:
  df[col]=2**df[col]


df.rename(columns={'SG13134300_257236337155_S002_GE1_1105_Oct12_1_2':"MS0_02",
'SG13134300_257236337155_S002_GE1_1105_Oct12_1_3':"MS0_03", 
'SG13134300_257236337155_S002_GE1_1105_Oct12_1_4':"MS0_04",
'SG13134300_257236337155_S002_GE1_1105_Oct12_2_1':"MS0_05",
'SG13134300_257236337155_S002_GE1_1105_Oct12_2_2':"MS0_06",
'SG13134300_257236337155_S002_GE1_1105_Oct12_2_3':"MS0_07"},inplace=True)



ComparisonDict={"SetA":"MS0_03:MS0_04",
                "SetB1":"MS0_03:MS0_06",
                "SetB2":"MS0_03:MS0_07",
                "SetC":"MS0_02:MS0_05",
                "SetD":"MS0_02:MS0_03"}
FC=3

for Str in ComparisonDict.items():
  Cont=Str[1].split(":")[0]
  Disease=Str[1].split(":")[1]
  Newcol1=Str[0]+"_l2FC_"+str(FC)
  Newcol2=Str[0]+"_updown_"+str(FC)
  df[Newcol1] = np.log2(np.nan_to_num(np.divide(df[Disease],df[Cont]), nan=1))
  conditions=[(df[Newcol1]>FC),(df[Newcol1]<-FC)]
  choices=["Up","Down"]
  df[Newcol2]=np.select(conditions,choices, default="NoChange")

df['SetB_updown_'+str(FC)]=np.where(df['SetB1_updown_'+str(FC)]==df['SetB2_updown_'+str(FC)],df['SetB2_updown_'+str(FC)],"NoChange")


updowncols=[x for x in df.columns if re.search("Set[A,B,C,D]_[updown]",x)]
lfccols=[x for x in df.columns if re.search("Set[A,B,C,D]_[l2FC]",x)]
othercols=[x for x in df.columns if x not in updowncols+lfccols]


df=df[othercols+lfccols+updowncols]
df["NonSignificant"]=df[[x for x in df.columns if re.search("Set[A,B,C,D]_updown_",x)]].isin(["NoChange"]).sum(axis=1)

df.to_csv("Diffrentially_Expressed_3lFC.csv",index=None)


df2=df.drop(['ensembl_gene_id','entrezgene','gene_biotype','external_gene_name','ProbeName', 'AgilentID'],axis=1)
df2.to_csv("Diffrentially_Expressed_3FC_Results.csv",index=None)

df2[[ x for x  in df2.columns if "updown" in  x]].apply(pd.Series.value_counts).drop(["SetB1_updown_3","SetB2_updown_3"],axis=1)[["SetA_updown_3","SetB_updown_3","SetC_updown_3","SetD_updown_3"]]


df[((df['SetA_updown_3']=="Down") | (df['SetB_updown_3']=="Down")) &  ((df['SetC_updown_3']=="Up") | (df['SetD_updown_3']=="Up"))  ]
df[((df['SetA_updown_3']=="Up") | (df['SetB_updown_3']=="Up")) &  ((df['SetC_updown_3']=="Down") | (df['SetD_updown_3']=="Down"))  ]

df[((df['SetA_updown_3']=="Down") & (df['SetB_updown_3']=="Down")) &  ((df['SetC_updown_3']=="Up") & (df['SetD_updown_3']=="Up"))  ]
df[((df['SetA_updown_3']=="Up") & (df['SetB_updown_3']=="Up")) &  ((df['SetC_updown_3']=="Down") & (df['SetD_updown_3']=="Down"))  ]


for Set in [x for x in df.columns if re.search("Set[A,B,C,D]_[updown]",x)]:
  gp = GProfiler(return_dataframe=True)
  UP=list(df[df[Set]=="Up"]["GeneName"].unique())
  gp = GProfiler(return_dataframe=True)
  test=gp.profile(organism='hsapiens',query=UP)
  Outname=Set.replace("_updown","")+"_Up_Regulated_Gene_Enrichment.csv"
  test.to_csv(Outname,index=None)
  
  gp = GProfiler(return_dataframe=True)
  DOWN=list(df[df[Set]=="Down"]["GeneName"].unique())
  gp = GProfiler(return_dataframe=True)
  test=gp.profile(organism='hsapiens',query=DOWN)
  Outname=Set.replace("_updown","")+"_Down_Regulated_Gene_Enrichment.csv"
  test.to_csv(Outname,index=None)
  
  gp = GProfiler(return_dataframe=True)
  DOWN_UP=list(df[df[Set]!="NoChange"]["GeneName"].unique())
  gp = GProfiler(return_dataframe=True)
  test=gp.profile(organism='hsapiens',query=DOWN_UP)
  Outname=Set.replace("_updown","")+"_UP_Down_Regulated_Gene_Enrichment.csv"
  test.to_csv(Outname,index=None)


for file in glob.glob("*Gene_Enrichment.csv"):
  test=pd.read_csv(file)
  print(file,test[test['significant']==True].shape[0])















###https://www.biostars.org/p/9485414/
targets <- project.bgcorrect.norm$targets
targets$group <- factor(targets$SetName, levels = c('set1','set2'))
expr <- project.bgcorrect.norm$E

design <- model.matrix(~ 0 + targets$group)
colnames(design) <- c('set1', 'set2')

# Fit the linear model on the study's data
project.fitmodel <- lmFit(expr,design)


output <- topTable(project.fitmodel, coef="set2-set1",
genelist=fit$genes,number=Inf)


# Applying the empirical Bayes method to the fitted values
# Acts as an extra normalisation step and aims to bring the different
#   probe-wise variances to common values
project.fitmodel.eBayes <- eBayes(project.fitmodel)\
names(project.fitmodel.eBayes)


# Make individual contrasts: cancer vs normal
res <- makeContrasts(res = 'set2-set1', levels = design)
res.fitmodel <- contrasts.fit(project.fitmodel.eBayes, res)
res.fitmodel.eBayes <- eBayes(res.fitmodel)




toptable <- topTable(res.fitmodel.eBayes,adjust = 'BH',coef = 'res',number = Inf,
              p.value = 1)


head(toptable)





not_ncbi=expression_outer2_biomart_lnc_ncbi[expression_outer2_biomart_lnc_ncbi['NCBI_GeneID'].isna()]

ncbi_complete=expression_outer2_biomart_lnc_ncbi[~expression_outer2_biomart_lnc_ncbi['NCBI_GeneID'].isna()]
ensembl_complete=not_ncbi[~not_ncbi['ensembl_gene_id'].isna()]

ens_transcript_complete= not_ens_trascript[~not_ens_trascript['Gene stable ID'].isna()][['ProbeName', 'GeneName', 'Comment[GENE_SYMBOL]','Gene stable ID', 'Transcript stable ID', 'Gene name','Gene type', 'NCBI gene (formerly Entrezgene) ID']]

not_intranscript=not_ens_trascript[not_ens_trascript['Gene stable ID'].isna()][["ProbeName","Comment[GENE_SYMBOL]"]]


ncbi_complete2=ncbi_complete[['ProbeName','GeneName','Comment[GENE_SYMBOL]',"NCBI_Symbol","NCBI_GeneID","NCBI_Ensembl","gene_biotype"]]
ensembl_complete2=ensembl_complete[["ProbeName","GeneName",'Comment[GENE_SYMBOL]',"external_gene_name","entrezgene_id","ensembl_gene_id","gene_biotype"]]
ens_transcript_complete2=ens_transcript_complete[['ProbeName', 'GeneName','Comment[GENE_SYMBOL]',"Gene stable ID","Gene name","Gene type","NCBI gene (formerly Entrezgene) ID"]]

not_intranscript



ncbi_complete2[["ProbeName","NCBI_Symbol","NCBI_GeneID","NCBI_Ensembl","gene_biotype"]].columns=['ProbeName',"Gene_synbol","NCBI_ID","Ensembl_ID","Biotype"]

ensembl_complete2.columns=['ProbeName',"Gene_synbol","NCBI_ID","Ensembl_ID","Biotype"]
ens_transcript_complete2.columns=['ProbeName',"Gene_synbol","NCBI_ID","Ensembl_ID","Biotype"]
not_intranscript.columns=['ProbeName',"Gene_synbol","NCBI_ID","Ensembl_ID","Biotype"]



['ProbeName',"Gene_synbol","NCBI_ID","Ensembl_ID","Biotype"]



gcard_cosmic,
onkokb_cancermine

metastasis_group


onkokb_cancermine_metastasis=pd.merge(metastasis_group,onkokb_cancermine,on=["Genesymbol","NCBI_ID"],how="outer")


onkokb_cancermine_metastasis_cosmic=pd.merge(onkokb_cancermine_metastasis,cosmic,on="NCBI_ID",how="outer")




expression=pd.read_csv("Diffrentially_Expressed_3lFC_with_Annotation.csv")
cancergene=pd.read_csv("All_cancer_Gene_List_drug.csv")



cancergene_None=cancergene[cancergene['NCBI_ID']=="None"]
cancergene=cancergene[cancergene['NCBI_ID']!="None"]

expression_None=expression[expression['NCBI_ID'].isna()]
expression=expression[~expression['NCBI_ID'].isna()]

expression_cancergene=pd.merge(expression,cancergene,on="NCBI_ID",how="left")
expression_cancergene






annot,
annot_none

cancer_None=cancer[cancer['NCBI_ID'].isna()]
cancer=cancer[~cancer['NCBI_ID'].isna()]


annot_cancer=pd.merge(annot,cancer,on="NCBI_ID",how="right")

