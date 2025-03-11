library(tidyverse)
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(R.utils)
library(grid)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(fastDummies)
library(viridis)
library(ggVennDiagram)
library(UpSetR)
library(pROC)

data_wd<-"/Users/erv894/OneDrive - Harvard University/PostDoc/Data/"
plot_wd<-"~/OneDrive - Harvard University/PostDoc/Data/variants/"
sourceDirectory(paste0(data_wd,"/miamiplot-master/R/"),modifiedOnly = FALSE)

load(paste0(data_wd,"genes_biomaRt_all.RData"))
genes_combine$mid_point<-ceiling((genes_combine$start_position+genes_combine$end_position)/2)
genes_combine$gene_length<-genes_combine$end_position-genes_combine$start_position + 1
genes_combine<-genes_combine%>%filter(gene_biotype=="protein_coding")
gene_loc<-genes_combine%>%filter(gene_biotype=="protein_coding",hgnc_symbol!="")%>%dplyr::select(hgnc_symbol,chromosome_name,gene_length,mid_point,start_position,end_position)
colnames(gene_loc)<-c("gene","chr","gene_length","mid_point","start","end")

ENCODE_cCREs_V3<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/TOPMED-cCREs/GRCh38-cCREs.bed",col_names = FALSE)
colnames(ENCODE_cCREs_V3)<-c("chr","start","end","other_id","cCRE_accession","element_classification")
ENCODE_cCREs_V3$ENCODE_element_type<-NA
ENCODE_cCREs_V3$ENCODE_element_type[grepl("dELS",ENCODE_cCREs_V3$element_classification)]<-"dELS"
ENCODE_cCREs_V3$ENCODE_element_type[grepl("pELS",ENCODE_cCREs_V3$element_classification)]<-"pELS"
ENCODE_cCREs_V3$ENCODE_element_type[grepl("PLS",ENCODE_cCREs_V3$element_classification)]<-"PLS"
ENCODE_cCREs_V3<-ENCODE_cCREs_V3%>%filter(!is.na(ENCODE_element_type))
ENCODE_cCREs_V3$chr[ENCODE_cCREs_V3$chr=="chrX"]<-"chr23"
ENCODE_cCREs_V3<-ENCODE_cCREs_V3%>%filter(!chr=="chrY")
ENCODE_cCREs_V3$chr<-as.numeric(gsub("chr","",ENCODE_cCREs_V3$chr))


Harvard_MGH_variants<-read_csv(paste0(data_wd,"/variants/0922list_annotations.csv"))
colnames(Harvard_MGH_variants)<-gsub(" ","_",colnames(Harvard_MGH_variants))
colnames(Harvard_MGH_variants)<-gsub("-","_",colnames(Harvard_MGH_variants))
colnames(Harvard_MGH_variants)[1:2]<-c("variant","chr")
Harvard_MGH_variants<-Harvard_MGH_variants%>%filter(!rsID=="rsID")
Harvard_MGH_variants$chr[Harvard_MGH_variants$chr=="X"]<-"23"
Harvard_MGH_variants$chr<-as.numeric(Harvard_MGH_variants$chr)
Harvard_MGH_variants$Funseq_Description[is.na(Harvard_MGH_variants$Funseq_Description)]<-"indel"



# get_genes<-function(position,chromosome,max_dist){
#   g<-unique(gene_loc%>%filter(chr==chromosome,abs(position-start)<max_dist|abs(position-end)<max_dist)%>%pull(gene))
#   return(paste0(g,collapse=","))
# }
# system.time({
#   g<-character(length=nrow(Harvard_MGH_variants))
#   for(i in 1:nrow(Harvard_MGH_variants)){
#     g[i]<-get_genes(Harvard_MGH_variants$Position[i],Harvard_MGH_variants$chr[i],500000)
#     if(i%%1000==0){print(i)}
#   }
# })
# Harvard_MGH_variants_genes_within_500000<-g
# save(Harvard_MGH_variants_genes_within_500000,file=paste0(data_wd,"/variants/Harvard_MGH_variants_genes_within_500000.RData"))

load(paste0(data_wd,"/variantsHarvard_MGH_variants_genes_within_10000.RData"))
load(paste0(data_wd,"/variantsHarvard_MGH_variants_genes_within_100000.RData"))
load(paste0(data_wd,"/variants/Harvard_MGH_variants_genes_within_500000.RData"))
Harvard_MGH_variants$genes_within_10000<-Harvard_MGH_variants_genes_within_10000
Harvard_MGH_variants$genes_within_100000<-Harvard_MGH_variants_genes_within_100000
Harvard_MGH_variants$genes_within_500000<-Harvard_MGH_variants_genes_within_500000

#This could/should probably be implemented using gRanges for faster code
# But this seemed much easier to understand and verify at the time
# get_ENCODE_element<-function(position,chromosome){
#   temp<-unique(ENCODE_cCREs_V3%>%filter(chr==chromosome,start<=position & position<=end)%>%pull(ENCODE_element_type))
#   #browser()
#   if(length(temp)==1){
#     return(paste0("ENCODE ",temp))
#   }else{
#     if(length(temp)==0){
#       return("Not In ENCODE cCRE")
#     }else{
#       #this would happen if overlapping elements in ENCODE
#       # there is only one, and it is both dELS,
#       return("Multiple Rows?")
#     }
# 
#   }
# }
# e<-character(length=nrow(Harvard_MGH_variants))
# all_pos<-Harvard_MGH_variants$Position
# all_chr<-Harvard_MGH_variants$chr
# for(i in 1:nrow(Harvard_MGH_variants)){
#   e[i]<-get_ENCODE_element(all_pos[i],all_chr[i])
#   if(i%%1000==0){print(i)}
# }
# Harvard_MGH_variants_ENCODE_elements<-e
# save(Harvard_MGH_variants_ENCODE_elements,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_variants_ENCODE_elements.RData")

load(paste0(data_wd,"/variantsHarvard_MGH_variants_ENCODE_elements.RData"))
Harvard_MGH_variants$ENCODE_element<-Harvard_MGH_variants_ENCODE_elements
Harvard_MGH_variants$ENCODE_element[Harvard_MGH_variants$ENCODE_element%in%c("ENCODE dELS","ENCODE pELS")]<-"ENCODE Enhancer"
Harvard_MGH_variants$ENCODE_element[Harvard_MGH_variants$ENCODE_element%in%c("ENCODE PLS")]<-"ENCODE Promoter"

Harvard_MGH_variants$ENCODE_element<-factor(Harvard_MGH_variants$ENCODE_element,levels=c("Not In ENCODE cCRE","ENCODE Enhancer","ENCODE Promoter"))
Harvard_MGH_variants$ENCODE_element[is.na(Harvard_MGH_variants$ENCODE_element)]<-"Not In ENCODE cCRE"
# Evaluated two ways to classify a variant as coding or not
# one way is to look for "exon" in Genecode_Comprehensive_Category
# other is to see if Genecode_Comprehensive_Exonic_Info has something in it
# ncRNA_exonic variants should be considered non-coding
# in the end, start w non-coding, consider coding based on two conditions
# and change back only if Gencode explicitly classifies as non-coding

Harvard_MGH_variants$variant_category<-"noncoding"
Harvard_MGH_variants$variant_category[grepl("exon",Harvard_MGH_variants$Genecode_Comprehensive_Category)
                              |!is.na(Harvard_MGH_variants$Genecode_Comprehensive_Exonic_Info) ]<-"coding"
Harvard_MGH_variants$variant_category[grepl("ncRNA_exonic",Harvard_MGH_variants$Genecode_Comprehensive_Category)]<-"noncoding"

Harvard_MGH_variants$nearest_gene_gencode<-Harvard_MGH_variants$Genecode_Comprehensive_Info
Harvard_MGH_variants$nearest_gene_gencode<-gsub("\\(.*","",Harvard_MGH_variants$nearest_gene_gencode)
Harvard_MGH_variants$nearest_gene_gencode<-gsub("\\;.*","",Harvard_MGH_variants$nearest_gene_gencode)
Harvard_MGH_variants$nearest_gene_gencode<-gsub("\\,.*","",Harvard_MGH_variants$nearest_gene_gencode)

Harvard_MGH_variants$CAGE_noncoding_category<-NA
Harvard_MGH_variants$CAGE_noncoding_category[Harvard_MGH_variants$variant_category=="noncoding"]<-"Not CAGE Prom. or Enh."
Harvard_MGH_variants$CAGE_noncoding_category[!is.na(Harvard_MGH_variants$CAGE_Promoter)]<-"CAGE Promoter"
Harvard_MGH_variants$CAGE_noncoding_category[!is.na(Harvard_MGH_variants$CAGE_Enhancer)]<-"CAGE Enhancer"

Harvard_MGH_variants$CAGE_noncoding_category2<-"Not CAGE Prom. or Enh."
#Harvard_MGH_variants$CAGE_noncoding_category2[Harvard_MGH_variants$variant_category=="noncoding"]<-"other_noncoding"
Harvard_MGH_variants$CAGE_noncoding_category2[!is.na(Harvard_MGH_variants$CAGE_Promoter)]<-"CAGE Promoter"
Harvard_MGH_variants$CAGE_noncoding_category2[!is.na(Harvard_MGH_variants$CAGE_Enhancer)]<-"CAGE Enhancer"

#Harvard_MGH_variants$CAGE_noncoding_category2<-factor(Harvard_MGH_variants$CAGE_noncoding_category2,levels=c("CAGE Promoter","CAGE Enhancer"),exclude = NULL )

Harvard_MGH_variants<-Harvard_MGH_variants%>%distinct()

# Add in MACIE scores

Harvard_MGH_MACIE_hg19<-read_tsv(paste0(data_wd,"/variants/Harvard_MGH_MACIE_scores.tsv"))
Harvard_MGH_MACIE_hg19<-Harvard_MGH_MACIE_hg19%>%dplyr::select(-chr,-pos,-ref,-alt)%>%distinct()
prob_PHRED<-function(x){return(pmax(-10*log10(1-x+.0001),0))}

scores<-data.frame(lapply(Harvard_MGH_MACIE_hg19%>%dplyr::select(-variant,-Consequence),prob_PHRED))
#colnames(scores)<-paste0(colnames(scores),"_PHRED")
Harvard_MGH_MACIE_hg19_PHRED<-bind_cols(Harvard_MGH_MACIE_hg19%>%dplyr::select(variant,Consequence),
                                scores)

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,Harvard_MGH_MACIE_hg19_PHRED,by="variant")

Harvard_MGH_variants<-Harvard_MGH_variants%>%mutate(MACIE_max_noncoding=pmax(MACIE_conserved,MACIE_anyclass,MACIE_regulatory,na.rm=TRUE))
Harvard_MGH_variants<-Harvard_MGH_variants%>%mutate(MACIE_max_coding=pmax(MACIE_conserved,MACIE_anyclass,MACIE_protein,na.rm=TRUE))



#Harvard_MGH_variants$Variant_Category<-Harvard_MGH_variants$CAGE_noncoding_category
#index<-!Harvard_MGH_variants$Genecode_Comprehensive_Exonic_Category=="noncoding"
#Harvard_MGH_variants$Variant_Category[index]<-Harvard_MGH_variants$Genecode_Comprehensive_Exonic_Category[index]

#https://stackoverflow.com/questions/49572416/r-convert-to-factor-with-order-of-levels-same-with-case-when
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}
create_aPC_categories<-function(obj,col_name,new_var_name=paste0(col_name,"_category")){
  obj<-obj%>%mutate(temp=fct_case_when((is.na(!!sym(col_name)) | !!sym(col_name)< 20) ~"(0, 20)",
                                                        !!sym(col_name)< 30 ~ "[20, 30)",
                                                        !!sym(col_name)> 30 ~ ">30"))
                                                        #!!sym(col_name)> 40 ~ ">40"))
  colnames(obj)[colnames(obj)=="temp"]<-new_var_name
  return(obj)
}

Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Protein_Function")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Transcription_Factor")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Conservation")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Epigenetics_Active")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Epigenetics_Repressed")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="aPC_Epigenetics_Transcription")

Harvard_MGH_variants<-Harvard_MGH_variants%>%mutate(TOPMed_AF_category=fct_case_when(is.na(TOPMed_Bravo_AF)~"unobserved",
                                  TOPMed_Bravo_AF<.01 ~"(0,.01)",
                                  TOPMed_Bravo_AF<.05 ~ "[.01,.05)",
                                  TOPMed_Bravo_AF<.25 ~ "[.05,.25)",
                                  TOPMed_Bravo_AF>=.25 ~ "[.25,1)"))

Harvard_MGH_variants$Drug_Response<-NA
Harvard_MGH_variants$Drug_Response[Harvard_MGH_variants$Clinical_Significance%in%c("drug_response")]<-"Drug Response"

Harvard_MGH_variants$Pathogenic<-NA
Harvard_MGH_variants$Pathogenic[Harvard_MGH_variants$Clinical_Significance%in%c("Pathogenic",
                                                                "Pathogenic/Likely_pathogenic"
                                                                ,"Likely_pathogenic")]<-"Pathogenic"

Harvard_MGH_variants$SIFTcat[is.na(Harvard_MGH_variants$SIFTcat)]<-"Unknown"
colnames(Harvard_MGH_variants)[colnames(Harvard_MGH_variants)=="SIFTcat"]<-"SIFT_cat"

Harvard_MGH_variants$near_gene_CAGE<-NA
index<-!is.na(Harvard_MGH_variants$CAGE_noncoding_category2) & !is.na(Harvard_MGH_variants$nearest_gene_gencode)
Harvard_MGH_variants$near_gene_CAGE[index]<-paste0(Harvard_MGH_variants$nearest_gene_gencode[index],": ",gsub("_"," ",Harvard_MGH_variants$CAGE_noncoding_category2[index]))

Harvard_MGH_variants$near_gene_ENCODE<-NA
index<-!is.na(Harvard_MGH_variants$ENCODE_element) & !is.na(Harvard_MGH_variants$nearest_gene_gencode)
Harvard_MGH_variants$near_gene_ENCODE[index]<-paste0(Harvard_MGH_variants$nearest_gene_gencode[index],": ",gsub("_"," ",Harvard_MGH_variants$ENCODE_element[index]))


compute_PHRED_scale<-function(x){
  return(10*-log10(rank(-x)/length(x)))
}

Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_conserved")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_anyclass")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_protein")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_regulatory")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_max_coding")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="MACIE_max_noncoding")

Harvard_MGH_variants$MACIE_conserved_g.t._20<-Harvard_MGH_variants$MACIE_conserved>20
Harvard_MGH_variants$MACIE_conserved_g.t._20[is.na(Harvard_MGH_variants$MACIE_conserved_g.t._20)]<-FALSE
Harvard_MGH_variants$MACIE_protein_g.t._20<-Harvard_MGH_variants$MACIE_protein>20
Harvard_MGH_variants$MACIE_protein_g.t._20[is.na(Harvard_MGH_variants$MACIE_protein_g.t._20)]<-FALSE
Harvard_MGH_variants$MACIE_anyclass_g.t._20<-Harvard_MGH_variants$MACIE_anyclass>20
Harvard_MGH_variants$MACIE_anyclass_g.t._20[is.na(Harvard_MGH_variants$MACIE_anyclass_g.t._20)]<-FALSE
Harvard_MGH_variants$MACIE_regulatory_g.t._20<-Harvard_MGH_variants$MACIE_regulatory>20
Harvard_MGH_variants$MACIE_regulatory_g.t._20[is.na(Harvard_MGH_variants$MACIE_regulatory_g.t._20)]<-FALSE
# 
Harvard_MGH_variants$aPC_Conservation_g.t._20<-Harvard_MGH_variants$aPC_Conservation>20
Harvard_MGH_variants$aPC_Conservation_g.t._20[is.na(Harvard_MGH_variants$aPC_Conservation>20)]<-FALSE

Harvard_MGH_variants$aPC_Epigenetics_Active_g.t._20<-Harvard_MGH_variants$aPC_Epigenetics_Active>20
Harvard_MGH_variants$aPC_Epigenetics_Active_g.t._20[is.na(Harvard_MGH_variants$aPC_Epigenetics_Active>20)]<-FALSE

Harvard_MGH_variants$aPC_Epigenetics_Transcription_g.t._20<-Harvard_MGH_variants$aPC_Epigenetics_Transcription>20
Harvard_MGH_variants$aPC_Epigenetics_Transcription_g.t._20[is.na(Harvard_MGH_variants$aPC_Epigenetics_Transcription>20)]<-FALSE

Harvard_MGH_variants$aPC_Transcription_Factor_g.t._20<-Harvard_MGH_variants$aPC_Transcription_Factor>20
Harvard_MGH_variants$aPC_Transcription_Factor_g.t._20[is.na(Harvard_MGH_variants$aPC_Transcription_Factor>20)]<-FALSE

Harvard_MGH_variants$aPC_Protein_Function_g.t._20<-Harvard_MGH_variants$aPC_Protein_Function>20
Harvard_MGH_variants$aPC_Protein_Function_g.t._20[is.na(Harvard_MGH_variants$aPC_Protein_Function>20)]<-FALSE
Harvard_MGH_variants<-Harvard_MGH_variants%>%
  mutate(max_aPC=pmax(aPC_Conservation,aPC_Epigenetics_Active,aPC_Epigenetics_Repressed
                      ,aPC_Epigenetics_Transcription
                      ,aPC_Protein_Function,aPC_Transcription_Factor))


Harvard_MGH_variants$DNase_PHRED<-compute_PHRED_scale(Harvard_MGH_variants$DNase)
Harvard_MGH_variants$H3K27ac_PHRED<-compute_PHRED_scale(Harvard_MGH_variants$H3K27ac)
Harvard_MGH_variants$H3K4me3_PHRED<-compute_PHRED_scale(Harvard_MGH_variants$H3K4me3)

Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="DNase_PHRED")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="H3K27ac_PHRED")
Harvard_MGH_variants<-create_aPC_categories(Harvard_MGH_variants,col_name="H3K4me3_PHRED")

index<-Harvard_MGH_variants$Genecode_Comprehensive_Exonic_Category=="nonsynonymous SNV"&Harvard_MGH_variants$MetaSVM_Score==TRUE
Harvard_MGH_variants$disruptive_missense<-NA
Harvard_MGH_variants$disruptive_missense[index]<-"Disruptive Missense"
#stop("260")
load(paste0("liver_ASE_sig_variants.RData"))
liver_ASE_sig_variants$liver_Allele_Specific_Evidence<-TRUE

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,liver_ASE_sig_variants,by="variant")
Harvard_MGH_variants$liver_Allele_Specific_Evidence[is.na(Harvard_MGH_variants$liver_Allele_Specific_Evidence)]<-FALSE

load(paste0(data_wd,"/variants/liver_ASE_sig_variants_noFDR.RData"))
liver_ASE_sig_variants_noFDR$liver_Allele_Specific_Evidence_noFDR<-TRUE

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,liver_ASE_sig_variants_noFDR,by="variant")
Harvard_MGH_variants$liver_Allele_Specific_Evidence_noFDR[is.na(Harvard_MGH_variants$liver_Allele_Specific_Evidence_noFDR)]<-FALSE
#stop("272")
load(paste0(data_wd,"/variants/Harvard_MGH_AS_elements.RData"))
load(paste0(data_wd,"/variants/Harvard_MGH_AS_information.RData"))
# 
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,Harvard_MGH_AS_elements,by="variant")
Harvard_MGH_variants$AS_PHRED_max[is.na(Harvard_MGH_variants$AS_PHRED_max)]<-0
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,Harvard_MGH_AS_information,by="variant")
# 
# # table<-tibble(thresh=numeric(),num_AS_variants=numeric())
# # j<-0
# # for(thresh in c(5E-2,5E-3,5E-4,5E-5,5E-6,5E-7,5E-8)){
# #   table<-bind_rows(table,c(thresh=thresh,num_AS_variants=sum(Harvard_MGH_variants$overall_min_p_ASE<thresh,na.rm=TRUE)))
# # }
# # table[,1]<-formatC(table%>%pull(1),format = "e", digits = 2)
# 
Harvard_MGH_variants$Allele_Specific_Evidence<-Harvard_MGH_variants$overall_min_p_ASE<1E-4
Harvard_MGH_variants$Allele_Specific_Evidence[is.na(Harvard_MGH_variants$Allele_Specific_Evidence)]<-FALSE
#stop("283")


load(paste0(data_wd,"/variants/Harvard_MGH_variants_scores_ENCFF536RJV.RData"))
to_merge<-Harvard_MGH_variants_scores_ENCFF536RJV%>%dplyr::select(variant,in_peak_ENCFF536RJV)
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,to_merge,by="variant")

load(paste0(data_wd,"/variants/Harvard_MGH_variants_scores_ENCFF262URW.RData"))
to_merge<-Harvard_MGH_variants_scores_ENCFF262URW%>%dplyr::select(variant,score_ENCFF262URW,PHRED_score_ENCFF262URW)
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,to_merge,by="variant")
###### 
#Catlas
#######

load(paste0(data_wd,"/variants/Harvard_MGH_variants_scores_Hepatocyte.RData"))
to_merge<-Harvard_MGH_variants_scores_Hepatocyte%>%dplyr::select(variant,in_peak_Hepatocyte)
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,to_merge,by="variant")
colnames(Harvard_MGH_variants)[colnames(Harvard_MGH_variants)=="in_peak_Hepatocyte"]<-"in_CATlas_peak_Hepatocyte"


unique_chrs<-unique(Harvard_MGH_variants$chr)
new_df<-tibble()

for(chr_in in unique_chrs)
{
  print(chr_in)
  cV2F_scores<-read_tsv(paste0("~/OneDrive - Harvard University/PostDoc/Data/ukbb_liver_annotation_chrmbpnet_LIVER/ukbb_0.9_0.01_liver_annotation_chrmbpnet_LIVER_lava_ld_ld_maf.",chr_in,".cv2f.txt"))%>%dplyr::select(-CM)%>%rename(chr=CHR,Position=BP,rsID=SNP,liver_cV2F=cV2F)
  df1<-Harvard_MGH_variants%>%filter(chr==chr_in)
  df1<-left_join(df1,cV2F_scores,by=c("chr","Position","rsID"))
  
  cV2F_scores<-read_tsv(paste0("~/OneDrive - Harvard University/PostDoc/Data/ukbb_all_baseline_annotation_chrmbpnet/ukbb_0.9_0.01_all_baseline_annotation_chrmbpnet_lava_ld_ld_maf.",chr_in,".cv2f.txt"))%>%dplyr::select(-CM)%>%rename(chr=CHR,Position=BP,rsID=SNP,cV2F=cV2F)
  df1<-left_join(df1,cV2F_scores,by=c("chr","Position","rsID"))
  new_df<-bind_rows(new_df,df1)
}
#stop("308")
Harvard_MGH_variants<-new_df
### Add in cV2F scores


# # Add in Kushal scores
# 
# LDL_Kushal<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/variants/LDL_TOPMed_9genes_0.8_sentinels_prox_finemapping.txt")
# HDL_Kushal<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/variants/HDL_TOPMed_9genes_0.8_sentinels_prox_finemapping.txt")
# all_Kushal<-bind_rows(LDL_Kushal,HDL_Kushal)%>%distinct()
# colnames(all_Kushal)[1]<-"variant"
# all_Kushal<-all_Kushal%>%dplyr::select(variant,cV2F.prob,cV2F.binary,Liver.cV2F.prob
#                                 ,Liver.cV2F.binary,LDL_PoPS_cS2G,LDL_PoPS_cS2G.rank
#                                 ,HDL_PoPS_cS2G,HDL_PoPS_cS2G.rank,regulomedb_ranking
#                                 ,regulomedb_prob)%>%
#                                 filter(!variant=="Variant (VCF)",!variant=="-",!cV2F.prob=="-")
# 
# cols<-c("cV2F.prob","cV2F.binary","Liver.cV2F.prob"
#         ,"Liver.cV2F.binary","LDL_PoPS_cS2G","LDL_PoPS_cS2G.rank"
#         ,"HDL_PoPS_cS2G","HDL_PoPS_cS2G.rank","regulomedb_ranking"
#         ,"regulomedb_prob")
# all_Kushal[,cols]<-apply(all_Kushal[,cols],2,FUN=as.numeric)
# 
# all_Kushal_original<-all_Kushal
# 
# scores<-data.frame(lapply(all_Kushal%>%dplyr::select(cV2F.prob,Liver.cV2F.prob),prob_PHRED))
# colnames(scores)<-gsub(".prob","",colnames(scores))
# #colnames(scores)<-paste0(colnames(scores),"_PHRED")
# all_Kushal<-bind_cols(all_Kushal%>%dplyr::select(-cV2F.prob,-Liver.cV2F.prob),
#                                 scores)
# 
# # dup<-all_Kushal2%>%filter(duplicated(variant))%>%pull(variant)
# # dup2<-all_Kushal%>%filter(variant%in%dup)%>%arrange(variant)
# # index<-which(!dup2[1,]==dup2[2,])
# # data.frame(dup2[1,index])
# # data.frame(dup2[2,index])
# 
# Harvard_MGH_variants<-left_join(Harvard_MGH_variants,all_Kushal,by="variant")
# #rm(all_Kushal)


# logical<-((!is.na(Harvard_MGH_variants$aPC_Conservation)&Harvard_MGH_variants$aPC_Conservation>20)|
#             (!is.na(Harvard_MGH_variants$aPC_Protein_Function)&Harvard_MGH_variants$aPC_Protein_Function>20)|
#             (!is.na(Harvard_MGH_variants$aPC_Epigenetics_Active)&Harvard_MGH_variants$aPC_Epigenetics_Active>20)|
#             (!is.na(Harvard_MGH_variants$aPC_Transcription_Factor)&Harvard_MGH_variants$aPC_Transcription_Factor>20)|
#             #(!is.na(Harvard_MGH_variants$MACIE_conserved_PHRED)&Harvard_MGH_variants$MACIE_conserved_PHRED>20)|
#             (!is.na(Harvard_MGH_variants$MACIE_anyclass_PHRED)&Harvard_MGH_variants$MACIE_anyclass_PHRED>20)|
#             #(!is.na(Harvard_MGH_variants$MACIE_regulatory_PHRED)&Harvard_MGH_variants$MACIE_regulatory_PHRED>20)|
#             #(!is.na(Harvard_MGH_variants$MACIE_protein_PHRED)&Harvard_MGH_variants$MACIE_protein_PHRED>20)|
#             (!is.na(Harvard_MGH_variants$Clinical_Significance) 
#              & Harvard_MGH_variants$Clinical_Significance%in%c("Pathogenic",
#                                                        "Pathogenic/Likely_pathogenic"
#                                                        ,"Likely_pathogenic"
#                                                        ,"drug_response"))|
#           !is.na(Harvard_MGH_variants$disruptive_missense)&Harvard_MGH_variants$disruptive_missense=="Disruptive Missense")

#chromBPnet_prediction<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/variants/chromBPnet_predictions.tsv")
#chromBPnet_prediction<-chromBPnet_prediction%>%mutate(variant=paste0(chr,"-",pos,"-",allele1,"-",allele2))%>%mutate(variant=gsub("chr","",variant))%>%dplyr::select(variant,rsid,logfc.mean)%>%dplyr::rename(rsID=rsid,chromBPnet_logfc_mu=logfc.mean)

#chromBPnet_prediction$abs_chromBPnet_logfc<-abs(chromBPnet_prediction$chromBPnet_logfc_mu)

loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }
  
}
#stop("382")
all_chromBPnet<-tibble()
for(file in list.files(path=paste0(data_wd,"/variants/ChromBPNet_Liver_and_HepG2/",pattern="chromBPnet_prediction3"))){
  t1<-loadRData(paste0(paste0(data_wd,"/variants/ChromBPNet_Liver_and_HepG2/",file)))
  all_chromBPnet<-bind_rows(all_chromBPnet,t1)
}
all_chromBPnet<-all_chromBPnet%>%filter(!grepl("\\+",variant1))%>%filter(!grepl("\\+",variant2))
full_chromBPnet_data<-all_chromBPnet
all_chromBPnet<-all_chromBPnet%>%filter(chromBPnet_new_func_defn==TRUE)%>%group_by(variant1,variant2)%>%mutate(all_chromBPnet_samples=paste0(sample_name,collapse=","))%>%dplyr::select(-sample_name,-jsd.mean.pval,-max_percentile.mean,-logfc.mean.pval)%>%distinct()%>%ungroup()

temp1<-inner_join(Harvard_MGH_variants,all_chromBPnet%>%dplyr::select(-variant2),by=c("variant"="variant1"))
temp2<-left_join(Harvard_MGH_variants,all_chromBPnet%>%dplyr::select(-variant1),by=c("variant"="variant2"))
temp2<-temp2%>%filter(!variant%in%temp1$variant)
#stop("394")

lol<-read_xlsx(paste0(data_wd,"/variants/1224_IGVF18lipidlocus_datasummary.xlsx"),sheet=2)

lol<-lol%>%mutate(variant=paste0(hg38Chromosome,"-",hg38Position,"-",hg38Ref,"-",hg38Alt))

 out_var<-Harvard_MGH_variants%>%filter(rsID%in%lol$rsID)%>%dplyr::select(rsID,all_of(contains("GENECODE")))%>%distinct()

colnames(out_var)<-gsub("Genecode","Gencode",colnames(out_var))

write_tsv(out_var,file=paste0(data_wd,"/variants/44_loci_coding_status.tsv"))

lol$rsID[!lol$rsID%in%out_var$rsID]

Harvard_MGH_variants<-bind_rows(temp1,temp2)
out_var2<-Harvard_MGH_variants%>%filter(rsID%in%lol$rsID)%>%distinct()


logical_aPC<-((!is.na(Harvard_MGH_variants$aPC_Conservation)&Harvard_MGH_variants$aPC_Conservation>20)|
            (!is.na(Harvard_MGH_variants$aPC_Protein_Function)&Harvard_MGH_variants$aPC_Protein_Function>20)|
            (!is.na(Harvard_MGH_variants$aPC_Epigenetics_Active)&Harvard_MGH_variants$aPC_Epigenetics_Active>20)|
            (!is.na(Harvard_MGH_variants$aPC_Epigenetics_Repressed)&Harvard_MGH_variants$aPC_Epigenetics_Repressed>20)|
            (!is.na(Harvard_MGH_variants$aPC_Epigenetics_Transcription)&Harvard_MGH_variants$aPC_Epigenetics_Transcription>20)|
            (!is.na(Harvard_MGH_variants$aPC_Transcription_Factor) &Harvard_MGH_variants$aPC_Transcription_Factor>20))
      
logical_aPC[is.na(logical_aPC)]<-FALSE
    
logical_MACIE<-(!is.na(Harvard_MGH_variants$MACIE_anyclass)&Harvard_MGH_variants$MACIE_anyclass>20)   
logical_MACIE[is.na(logical_MACIE)]<-FALSE
### TLand 

TLand_scores<-read_tsv(paste0(data_wd,"/variants/Harvard_MGH_variants_hg38_SNVs_TLand.tsv"))
colnames(TLand_scores)<-gsub(" ","_",colnames(TLand_scores))
TLand_scores<-TLand_scores%>%dplyr::select(-contains("light"))%>%mutate(variant=paste0(chrom,"-",end,"-",ref,"-",alt))
TLand_scores$variant<-gsub("chr","",TLand_scores$variant)

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,TLand_scores
                                ,by="variant")

Harvard_MGH_variants$liver_TLand_top5pct<-Harvard_MGH_variants$liver_TLand>=quantile(Harvard_MGH_variants$liver_TLand,.95,na.rm=TRUE)

  
logical_ClinVar<-(!is.na(Harvard_MGH_variants$Clinical_Significance) 
             & Harvard_MGH_variants$Clinical_Significance%in%c("Pathogenic",
                                                       "Pathogenic/Likely_pathogenic"
                                                       ,"Likely_pathogenic"
                                                       ,"drug_response"))|
            (!is.na(Harvard_MGH_variants$disruptive_missense) &Harvard_MGH_variants$disruptive_missense=="Disruptive Missense"
)
logical_ClinVar[is.na(logical_ClinVar)]<-FALSE

logical_ASE<-Harvard_MGH_variants$Allele_Specific_Evidence
logical_ASE[is.na(logical_ASE)]<-FALSE

logical_liver_ASE<-Harvard_MGH_variants$liver_Allele_Specific_Evidence
logical_liver_ASE[is.na(logical_liver_ASE)]<-FALSE

#logical_chromBPnet<-Harvard_MGH_variants$abs_chromBPnet_logfc_top10pct
#logical_chromBPnet[is.na(logical_chromBPnet)]<-FALSE

logical_chromBPnet<-Harvard_MGH_variants$chromBPnet_new_func_defn
logical_chromBPnet[is.na(logical_chromBPnet)]<-FALSE

#Harvard_MGH_variants%>%group_by(variant)%>%filter(n()>1)%>%ungroup()%>%dplyr::select(variant,max_percentile.mean,jsd.mean.pval,logfc.mean.pval)

#stop("455")
logical_liver_TLand<-Harvard_MGH_variants$liver_TLand_top5pct
logical_liver_TLand[is.na(logical_liver_TLand)]<-FALSE
#stop("411")
logical_liver_cV2F<-Harvard_MGH_variants$liver_cV2F>.75
logical_liver_cV2F[is.na(Harvard_MGH_variants$liver_cV2F)]<-FALSE

logical_cV2F<-Harvard_MGH_variants$cV2F>.75
logical_cV2F[is.na(Harvard_MGH_variants$cV2F)]<-FALSE

logical<-logical_aPC|logical_ClinVar|logical_MACIE|logical_liver_ASE|logical_chromBPnet|logical_liver_TLand|logical_liver_cV2F|logical_cV2F

Harvard_MGH_variants$Predicted_Functional<-logical
Harvard_MGH_variants$Predicted_Functional_aPC<-logical_aPC
Harvard_MGH_variants$Predicted_Functional_MACIE<-logical_MACIE
Harvard_MGH_variants$Predicted_Functional_ClinVar<-logical_ClinVar
Harvard_MGH_variants$Predicted_Functional_ASE<-logical_ASE
Harvard_MGH_variants$Predicted_Functional_liver_ASE<-logical_liver_ASE
Harvard_MGH_variants$Predicted_Functional_chromBPnet<-logical_chromBPnet
Harvard_MGH_variants$Predicted_Functional_liver_TLand<-logical_liver_TLand
Harvard_MGH_variants$Predicted_Functional_liver_cV2F<-logical_liver_cV2F
Harvard_MGH_variants$Predicted_Functional_cV2F<-logical_cV2F

Harvard_MGH_variants<-Harvard_MGH_variants%>%rowwise()%>%mutate(sum_pred_func=sum(Predicted_Functional_aPC,Predicted_Functional_MACIE,Predicted_Functional_liver_ASE,Predicted_Functional_chromBPnet,Predicted_Functional_ClinVar,Predicted_Functional_liver_TLand,Predicted_Functional_liver_cV2F,na.rm=TRUE))%>%ungroup()

#table(Harvard_MGH_variants$Predicted_Functional,Harvard_MGH_variants$abs_chromBPnet_logfc_mu_top10pct)

#table(Harvard_MGH_variants$Predicted_Functional,Harvard_MGH_variants$abs_chromBPnet_logfc_mu_top1pct)

#1+"e"

load(paste0(data_wd,"/variants/Harvard_MGH_DHS_elements.RData"))
Harvard_MGH_DHS_elements<-Harvard_MGH_DHS_elements%>%dplyr::select(-promoter_linked_gene,-enhancer_linked_gene)%>%distinct()
Harvard_MGH_DHS_elements2<-Harvard_MGH_DHS_elements%>%group_by(variant)%>%mutate(DHS_promoter=max(DHS_promoter)
                                                                ,DHS_enhancer=max(DHS_enhancer))%>%distinct()

Harvard_MGH_DHS_elements2$DHS_category<-"None"
Harvard_MGH_DHS_elements2$DHS_category[Harvard_MGH_DHS_elements2$DHS_enhancer==1]<-"DHS enhancer"
Harvard_MGH_DHS_elements2$DHS_category[Harvard_MGH_DHS_elements2$DHS_promoter==1]<-"DHS promoter"

table(Harvard_MGH_DHS_elements2$DHS_category)

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,Harvard_MGH_DHS_elements2%>%dplyr::select(-chr,-Position),by="variant")
Harvard_MGH_variants$DHS_category[is.na(Harvard_MGH_variants$DHS_category)]<-"Not DHS Prom. or Enh."

Harvard_MGH_variants$DHS_category<-factor(Harvard_MGH_variants$DHS_category,levels=c("Not DHS Prom. or Enh.","DHS enhancer"
                                                                     ,"DHS promoter"))


table(Harvard_MGH_variants$in_peak_ENCFF536RJV,
      Harvard_MGH_variants$Predicted_Functional)

cor(Harvard_MGH_variants$in_peak_ENCFF536RJV,
      Harvard_MGH_variants$PHRED_score_ENCFF262URW)
Harvard_MGH_variants%>%group_by(in_peak_ENCFF536RJV)%>%
  summarise(mean_PHRED_score=mean(PHRED_score_ENCFF262URW),
            mean_score=mean(score_ENCFF262URW))

IGV_out_MGH<-Harvard_MGH_variants%>%dplyr::select(chr,Position)%>%
  dplyr::rename(start=Position)
IGV_out_MGH$end<-IGV_out_MGH$start+1
IGV_out_MGH$chr<-paste0("chr",IGV_out_MGH$chr)
IGV_out_MGH$name<-"."
IGV_out_MGH$score<-1000
IGV_out_MGH$strand<-"."
write_tsv(IGV_out_MGH,file=paste0(data_wd,"/variants/IGV_out_MGH.bed",col_names=FALSE))
colnames(Harvard_MGH_variants)[colnames(Harvard_MGH_variants)=="in_peak_ENCFF536RJV"]<-"in_HepG2_ATAC-seq_peak"

Harvard_MGH_variants<-Harvard_MGH_variants%>%mutate(functional_category=as.factor(ifelse(Predicted_Functional,"Predicted Functional","Predicted Non-Functional")))

Harvard_MGH_variants<-Harvard_MGH_variants%>%ungroup()%>%distinct()



Harvard_MGH_variants_merge<-Harvard_MGH_variants%>%separate(variant,c("chr","start","ref","alt"))%>%
  mutate(prioritizedVariantID=paste0(chr,"_",Position,"_hg38_",ref,"_",alt))%>%dplyr::select(prioritizedVariantID,rsID,Predicted_Functional,`in_HepG2_ATAC-seq_peak`)


midyear_variants_in<-read_xlsx(paste0(data_wd,"/variants/IGVF_18_lipidlocus_variants_041123.xlsx"),sheet="18locus_variants")%>%dplyr::select(-"...31")


midyear_variants_out<-left_join(midyear_variants_in,Harvard_MGH_variants_merge,by=c("prioritizedVariantID","rsID"))

midyear_variants_out$Predicted_Functional[is.na(midyear_variants_out$Predicted_Functional)]<-FALSE

midyear_variants_out$`in_HepG2_ATAC-seq_peak`[is.na(midyear_variants_out$`in_HepG2_ATAC-seq_peak`)]<-FALSE

midyear_variants_out<-midyear_variants_out%>%dplyr::rename(FAVOR_Predicted_Functional=Predicted_Functional)

write_xlsx(midyear_variants_out,path=paste0(data_wd,"/variants/IGVF_midyear_variants_out.xlsx"))

tissue_cCREs<-read_csv(paste0(data_wd,"/variants/0922list_annotations_updatedwtissuecCREs.csv")
                       ,col_names=TRUE)

colnames(tissue_cCREs)<-gsub(" ","_",colnames(tissue_cCREs))
colnames(tissue_cCREs)<-gsub("\\.","_",colnames(tissue_cCREs))

tissue_cCREs<-tissue_cCREs%>%
  dplyr::rename(in_ENCODE_cts_adipose_cCRE=adipose_cCRE
                ,in_ENCODE_cts_liver_cCRE=liver_cCRE
                ,in_ENCODE_cts_pancreas_cCRE=pancreas_cCRE
                ,in_ENCODE_cts_small_intestine_cCRE=small_intestine_cCRE
                ,variant=`Variant_(VCF)`)%>%
  dplyr::select(in_ENCODE_cts_adipose_cCRE,in_ENCODE_cts_liver_cCRE
                ,in_ENCODE_cts_pancreas_cCRE,in_ENCODE_cts_small_intestine_cCRE
                ,variant)
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,tissue_cCREs
                                ,by="variant")
Harvard_MGH_variants<-Harvard_MGH_variants%>%ungroup()%>%rowwise()%>%mutate(in_Hep_ATAC_peak = in_CATlas_peak_Hepatocyte|`in_HepG2_ATAC-seq_peak`)%>%ungroup()
#stop("568")
###########
# caQTL data
###########
#caQTL_variants<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/variants/caQTL_variants_alreadyInLipidStudy.txt",col_names = FALSE)%>%pull(1)
#caQTL_variants<-gsub("\\:.*","",caQTL_variants)
caQTL_variants<-read_tsv(paste0(data_wd,"/variants/caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed"),col_names = TRUE)
colnames(caQTL_variants)[1]<-"chr"

Harvard_MGH_variants$liver_caQTL<-FALSE
Harvard_MGH_variants$liver_caQTL[Harvard_MGH_variants$rsID%in%c(caQTL_variants$proxy_rsID,caQTL_variants$lead_rsID)]<-TRUE

caQTL_eQTL_coloc<-read_xlsx(paste0(data_wd,"/variants/tableS15_caQTL_eQTL_colocalization_liver_caQTL.xlsx"),skip=2)

sum(Harvard_MGH_variants$rsID%in%caQTL_eQTL_coloc$eQTL_rsID)
sum(Harvard_MGH_variants$rsID%in%caQTL_eQTL_coloc$caQTL_rsID)
Harvard_MGH_variants$caQTL_eQTL_coloc<-FALSE
Harvard_MGH_variants$caQTL_eQTL_coloc[Harvard_MGH_variants$rsID%in%c(caQTL_eQTL_coloc$eQTL_rsID,caQTL_eQTL_coloc$caQTL_rsID)]<-TRUE

#stop("528")
########
# UKB Results
########

load(paste0(data_wd,"/variants/Harvard_MGH_variants_UKB_results.RData"))

Harvard_MGH_variants<-left_join(Harvard_MGH_variants,Harvard_MGH_variants_UKB_results,by="variant")

Harvard_MGH_variants$sig_pvalue_LDL_UKB_200K<-Harvard_MGH_variants$UKB_200K_LDL_pvalue<5E-8
Harvard_MGH_variants$sig_pvalue_LDL_UKB_200K[is.na(Harvard_MGH_variants$sig_pvalue_LDL_UKB_200K)]<-FALSE
#stop("539")
########
# TOPMed F8 Results
#########
#TOPMed_Lipids_results<-read_xlsx("~/OneDrive - Harvard University/PostDoc/Data/variants/TOPMed_Lipids_supp_table.xlsx",sheet = "Supplementary Data 4",skip=3)
#colnames(TOPMed_Lipids_results)[5]<-c("rsID")
TOPMed_LDL_sig_results<-read_delim(paste0(data_wd,"/variants/LDL_SigPos_Full_ForEric.txt.gz",delim=" "))
TOPMed_LDL_merge<-TOPMed_LDL_sig_results
TOPMed_LDL_merge$CHR<-as.numeric(gsub("chr","",TOPMed_LDL_merge$CHR))
TOPMed_LDL_merge<-TOPMed_LDL_merge%>%dplyr::select(CHR,POS,Allele1,Allele2,p.value)
colnames(TOPMed_LDL_merge)<-c("chr","Position","TOPMed_Allele1","TOPMed_Allele2","TOPMed_F8_p.value_LDL")
Harvard_MGH_variants<-left_join(Harvard_MGH_variants,TOPMed_LDL_merge
                                 ,by=c("chr","Position"))

#TOPMed_Lipids_results<-TOPMed_Lipids_results[!is.na(TOPMed_Lipids_results$CHR),]
#TOPMed_Lipids_results<-TOPMed_Lipids_results[!is.na(TOPMed_Lipids_results$rsID),]
#TOPMed_Lipids_results$p.value<-as.numeric(TOPMed_Lipids_results$p.value)
#TOPMed_Lipids_results$BETA<-as.numeric(TOPMed_Lipids_results$BETA)
#TOPMed_Lipids_results$SE<-as.numeric(TOPMed_Lipids_results$SE)
#TOPMed_Lipids_results<-TOPMed_Lipids_results%>%dplyr::select(rsID,BETA,SE,p.value,MAF,pheno)


#TOPMed_Lipids_results_LDL<-TOPMed_Lipids_results%>%filter(pheno=="LDL")
#TOPMed_Lipids_results_LDL<-TOPMed_Lipids_results_LDL%>%dplyr::select(rsID,p.value,BETA,SE)%>%rename(TOPMed_F8_p.value_LDL=p.value,TOPMed_F8_BETA_LDL=BETA,TOPMed_F8_SE_LDL=SE)

#TOPMed_Lipids_results_HDL<-TOPMed_Lipids_results%>%filter(pheno=="HDL")
#TOPMed_Lipids_results_HDL<-TOPMed_Lipids_results_HDL%>%dplyr::select(rsID,p.value,BETA,SE)%>%rename(TOPMed_F8_p.value_HDL=p.value,TOPMed_F8_BETA_HDL=BETA,TOPMed_F8_SE_HDL=SE)

#Harvard_MGH_variants<-left_join(Harvard_MGH_variants,TOPMed_Lipids_results_LDL,by=c("rsID"))


Harvard_MGH_variants$sig_pvalue_LDL_TOPMed_F8<-Harvard_MGH_variants$TOPMed_F8_p.value_LDL<5E-8
Harvard_MGH_variants$sig_pvalue_LDL_TOPMed_F8[is.na(Harvard_MGH_variants$sig_pvalue_LDL_TOPMed_F8)]<-FALSE

#Harvard_MGH_variants$sig_pvalue_HDL_TOPMed_F8<-Harvard_MGH_variants$TOPMed_F8_p.value_HDL<5E-8
#Harvard_MGH_variants$sig_pvalue_HDL_TOPMed_F8[is.na(Harvard_MGH_variants$sig_pvalue_HDL_TOPMed_F8)]<-FALSE
#stop("575")


########
#Load in preliminary LDL results
########
uptake_results<-read_xlsx(paste0(data_wd,"/variants/062923_bean_element_result.LDLuptake_simple.xlsx"))%>%filter(!target_type%in%c("PosControl","NegControl"))%>%
  select(target_id,mu_z_adj,prioritizedVariantID
         ,`Gene/chrom`,`exon/position`,libraryName,finemappingscore_UKB_SUSIE)%>%
  dplyr::rename(LDL_uptake_mu_z_adj=mu_z_adj,chr=`Gene/chrom`,position=`exon/position`)%>%mutate(LDL_uptake_p=pnorm(2*-1*abs(LDL_uptake_mu_z_adj)))

efflux_results<-read_xlsx(paste0(data_wd,"/variants/062923_bean_element_result.TopFluor_simple.xlsx"))%>%filter(!target_type%in%c("PosControl","NegControl"))%>%
  select(target_id,mu_z_adj,prioritizedVariantID,`Gene/chrom`,`exon/position`)%>%
  dplyr::rename(LDL_efflux_mu_z_adj=mu_z_adj,chr=`Gene/chrom`,position=`exon/position`)%>%mutate(LDL_efflux_p=pnorm(2*-1*abs(LDL_efflux_mu_z_adj)))

experiment_results<-full_join(efflux_results,uptake_results,by=c("target_id","prioritizedVariantID","chr","position"))
experiment_results$position<-as.numeric(experiment_results$position)
experiment_results$chr<-as.numeric(experiment_results$chr)
experiment_results$LDL_efflux_p_FDR<-p.adjust(experiment_results$LDL_efflux_p,'fdr')
experiment_results$LDL_uptake_p_FDR<-p.adjust(experiment_results$LDL_uptake_p,'fdr')
experiment_results$finemappingscore_UKB_SUSIE<-as.numeric(experiment_results$finemappingscore_UKB_SUSIE)
sum(experiment_results$finemappingscore_UKB_SUSIE>.9,na.rm=TRUE)

Harvard_MGH_variants_w_results<-right_join(Harvard_MGH_variants,experiment_results,by=c("Position"="position"
                                                             ,"chr"="chr"))
Harvard_MGH_variants_w_results$LDL_uptake_p_FDR_lt_.05<-Harvard_MGH_variants_w_results$LDL_uptake_p_FDR<.05

Harvard_MGH_variants_w_results$LDL_efflux_p_FDR_lt_.05<-Harvard_MGH_variants_w_results$LDL_efflux_p_FDR<.05

Harvard_MGH_variants_w_results$neg_log10_LDL_efflux_p_FDR<--1*log10(Harvard_MGH_variants_w_results$LDL_efflux_p_FDR)

Harvard_MGH_variants_w_results$neg_log10_LDL_uptake_p_FDR<--1*log10(Harvard_MGH_variants_w_results$LDL_uptake_p_FDR)

# Variants that are not in original files?
# Some may be indels
temp2<-Harvard_MGH_variants_w_results%>%filter(!is.na(Predicted_Functional))

temp3<-temp2%>%filter(LDL_uptake_p_FDR<.05 | LDL_efflux_p_FDR<.05)
temp4<-temp2%>%filter(!(LDL_uptake_p_FDR<.05 | LDL_efflux_p_FDR<.05))
temp5<-temp2%>%filter(`in_HepG2_ATAC-seq_peak`)

round(prop.table(table(temp2$LDL_uptake_p_FDR<.05,temp2$LDL_efflux_p_FDR<.05)),3)

round(prop.table(table(temp2$LDL_uptake_p_FDR<.05,temp2$Predicted_Functional)),2)
round(prop.table(table(temp2$LDL_efflux_p_FDR<.05,temp2$Predicted_Functional)),2)

round(prop.table(table(temp5$LDL_uptake_p_FDR<.05,temp5$Predicted_Functional)),2)
round(prop.table(table(temp5$LDL_efflux_p_FDR<.05,temp5$Predicted_Functional)),2)

round(prop.table(table(temp5$LDL_efflux_p_FDR<.05,temp5$Predicted_Functional)),2)

temp3<-temp2%>%filter(Predicted_Functional==TRUE)

round(prop.table(table(temp3$LDL_efflux_p_FDR<.05|temp3$LDL_uptake_p_FDR<.05,temp3$Predicted_Functional)),2)

#prop.table(table(temp2$LDL_uptake_p_FDR<.05,temp2$abs_chromBPnet_logfc_top10pct))
#prop.table(table(temp2$LDL_efflux_p_FDR<.05,temp2$abs_chromBPnet_logfc_top10pct))

table(temp2%>%arrange(LDL_efflux_p_FDR)%>%dplyr::slice(1:50)%>%pull(Predicted_Functional))
table(temp2%>%arrange(LDL_uptake_p_FDR)%>%dplyr::slice(1:50)%>%pull(Predicted_Functional))


out_var2<-Harvard_MGH_variants%>%filter(rsID%in%lol$rsID)%>%distinct()
#stop("692")
table(Harvard_MGH_variants_w_results$sum_pred_func,(Harvard_MGH_variants_w_results$LDL_efflux_p_FDR_lt_.05&Harvard_MGH_variants_w_results$LDL_uptake_p_FDR_lt_.05))



out_var3<-Harvard_MGH_variants_w_results%>%filter(rsID%in%lol$rsID)%>%dplyr::select(variant,rsID,LDL_efflux_p_FDR,LDL_uptake_p_FDR,all_of(contains("Predicted_Functional")),sum_pred_func)
#stop("698")
write_tsv(out_var3,file = paste0(data_wd,"/variants/44_loci_functional_data.tsv"))
#stop("700")
out_df<-Harvard_MGH_variants_w_results%>%dplyr::select(variant,rsID,LDL_efflux_p,LDL_efflux_mu_z_adj,LDL_uptake_mu_z_adj,LDL_uptake_p,sum_pred_func,all_of(contains("Predicted_Functional")))

colnames(out_df)[grepl("Predicted_Functional",colnames(out_df))]<-paste0(colnames(out_df)[grepl("Predicted_Functional",colnames(out_df))],"_BE")
#stop("704")
write_tsv(out_df,file=paste0(data_wd,"/variants/BE_data_full_results.tsv"))

rofl<-full_chromBPnet_data%>%filter(variant1%in%lol$variant|variant2%in%lol$variant)
vars<-unique(c(rofl$variant1,rofl$variant2))

sum(vars%in%lol$variant) # 40 variants match

to_merge<-lol%>%select(variant,rsID)

temp<-left_join(rofl,to_merge,by=c("variant1"="variant"))
temp<-left_join(temp,to_merge,by=c("variant2"="variant"))
temp$rsID<-temp$rsID.x
temp$rsID[is.na(temp$rsID)&!is.na(temp$rsID.y)]<-temp$rsID.y[is.na(temp$rsID)&!is.na(temp$rsID.y)]
temp<-temp%>%dplyr::select(-rsID.x,-rsID.y)
temp<-temp%>%dplyr::select(variant1,variant2,rsID,everything())
temp<-temp%>%rename(ENCODE_sample_name=sample_name)

write_tsv(temp,file = paste0(data_wd,"/variants/44_loci_chromBPnet_data.tsv"))

lol$rsID[!lol$rsID%in%temp$rsID]

rs<-unique(c(lol$rsID))

out<-Harvard_MGH_variants_w_results%>%select(variant,all_of(contains("Predicted_Functional")))%>%summarise(n())



#stop("722")

all_var<-read_xlsx(paste0(data_wd,"/variants/0922_Lipid_Topmed_22gene_08sentinels_prox_finemapping_eQTL.xlsx"))%>%dplyr::select(hg38Chromosome,hg38Position,hg38Ref,hg38Alt,prioritizedVariantID,rsID)
colnames(all_var)<-c("chr","pos","ref","alt","variant","rsID")
write_csv(all_var,file=paste0(data_wd,"/variants/all_Harvard_MGH_variants.csv"))

# var1<-Harvard_MGH_variants_w_results%>%dplyr::select(chr,Position,ref.x,alt.x,variant,rsID,Predicted_Functional,LDL_efflux_p_FDR_lt_.05,LDL_uptake_p_FDR_lt_.05)
# colnames(var1)[1:6]<-c("chr","position","ref","alt","variant","rsID")
# write_csv(var1,file="~/OneDrive - Harvard University/PostDoc/Data/variants/all_Harvard_MGH_variants.csv")

var1<-Harvard_MGH_variants_w_results%>%filter((LDL_efflux_p_FDR_lt_.05 | LDL_uptake_p_FDR_lt_.05)&!Predicted_Functional)%>%dplyr::select(chr,Position,ref.x,alt.x,variant,rsID,Predicted_Functional,LDL_efflux_p_FDR_lt_.05,LDL_uptake_p_FDR_lt_.05)
colnames(var1)[1:6]<-c("chr","position","ref","alt","variant","rsID")
write_csv(var1,file=paste0(data_wd,"/variants/Harvard_MGH_variants_sig_experimental_result_no_functional_evidence.csv"))

var1<-Harvard_MGH_variants_w_results%>%filter((LDL_efflux_p_FDR_lt_.05 | LDL_uptake_p_FDR_lt_.05)&Predicted_Functional)%>%dplyr::select(chr,Position,ref.x,alt.x,variant,rsID,Predicted_Functional,LDL_efflux_p_FDR_lt_.05,LDL_uptake_p_FDR_lt_.05)
colnames(var1)[1:6]<-c("chr","position","ref","alt","variant","rsID")
write_csv(var1,file=paste0(data_wd,"/variants/Harvard_MGH_variants_sig_experimental_result_and_functional_evidence.csv"))


var1<-Harvard_MGH_variants_w_results%>%filter(!LDL_efflux_p_FDR_lt_.05 & !LDL_uptake_p_FDR_lt_.05&Predicted_Functional)%>%dplyr::select(chr,Position,ref.x,alt.x,variant,rsID,Predicted_Functional,LDL_efflux_p_FDR_lt_.05,LDL_uptake_p_FDR_lt_.05)
colnames(var1)[1:6]<-c("chr","position","ref","alt","variant","rsID")
write_csv(var1,file=paste0(data_wd,"/variants/Harvard_MGH_variants_no_sig_experimental_result_with_functional_evidence.csv"))

save(Harvard_MGH_variants_w_results,file=paste0(data_wd,"/variants/Harvard_MGH_variants_w_results.RData"))

#stop("739")
Harvard_MGH_variants_w_results<-Harvard_MGH_variants_w_results%>%filter(!is.na(Predicted_Functional))
Harvard_MGH_variants_w_results$either_sig<-Harvard_MGH_variants_w_results$LDL_efflux_p_FDR_lt_.05|Harvard_MGH_variants_w_results$LDL_uptake_p_FDR_lt_.05

recall_df<-Harvard_MGH_variants_w_results%>%
  mutate(TruePositive_Overall=(either_sig&Predicted_Functional)
         ,TruePositive_chromBPnet=(either_sig&Predicted_Functional_chromBPnet)
         ,TruePositive_aPC=(either_sig&Predicted_Functional_aPC)
         ,TruePositive_MACIE=(either_sig&Predicted_Functional_MACIE)
         ,TruePositive_cV2F=(either_sig&Predicted_Functional_cV2F)
         ,TruePositive_liver_cV2F=(either_sig&Predicted_Functional_liver_cV2F)
         ,TruePositive_liver_ASE=(either_sig&Predicted_Functional_liver_ASE)
         ,TruePositive_liver_TLand=(either_sig&Predicted_Functional_liver_TLand))%>%
  mutate(FalseNegative_Overall=(either_sig&(!Predicted_Functional))
         ,FalseNegative_chromBPnet=(either_sig&(!Predicted_Functional_chromBPnet))
         ,FalseNegative_liver_ASE=(either_sig&(!Predicted_Functional_liver_ASE))
         ,FalseNegative_aPC=(either_sig&(!Predicted_Functional_aPC))
         ,FalseNegative_MACIE=(either_sig&(!Predicted_Functional_MACIE))
         ,FalseNegative_cV2F=(either_sig&(!Predicted_Functional_cV2F))
         ,FalseNegative_liver_cV2F=(either_sig&(!Predicted_Functional_liver_cV2F))
         ,FalseNegative_liver_TLand=(either_sig&(!Predicted_Functional_liver_TLand)))%>%
  mutate(FalsePositive_Overall=((!either_sig)&Predicted_Functional)
         ,FalsePositive_chromBPnet=((!either_sig)&Predicted_Functional_chromBPnet)
         ,FalsePositive_liver_ASE=((!either_sig)&Predicted_Functional_liver_ASE)
         ,FalsePositive_aPC=((!either_sig)&Predicted_Functional_aPC)
         ,FalsePositive_MACIE=((!either_sig)&Predicted_Functional_MACIE)
         ,FalsePositive_cV2F=((!either_sig)&Predicted_Functional_cV2F)
         ,FalsePositive_liver_cV2F=((!either_sig)&Predicted_Functional_liver_cV2F)
         ,FalsePositive_liver_TLand=((!either_sig)&Predicted_Functional_liver_TLand))

for(category in c("Overall","aPC","MACIE","chromBPnet"
                  ,"cV2F","liver_cV2F","liver_TLand","liver_ASE")){
  tp<-recall_df[,paste0("TruePositive_",category)]%>%pull(1)
  fp<-recall_df[,paste0("FalsePositive_",category)]%>%pull(1)
  fn<-recall_df[,paste0("FalseNegative_",category)]%>%pull(1)
  precision<-sum(tp)/(sum(tp)+sum(fp))
  recall<-sum(tp)/(sum(tp)+sum(fn))
  precision_boot<-c()
  recall_boot<-c()
  for(nboot in 1:1000){
    if(nboot%%100==0){print(nboot)}
    idx = sample(1:length(tp), length(tp), replace = T)
    tp_boot<-tp[idx]
    fp_boot<-fp[idx]
    fn_boot<-fn[idx]
    
    precision_boot<-c(precision_boot,sum(tp_boot)/(sum(tp_boot)+sum(fp_boot)))
    recall_boot<-c(recall_boot,sum(tp_boot)/(sum(tp_boot)+sum(fn_boot)))
  }
  recall_df[,paste0("Recall_",category)]<-recall
  recall_df[,paste0("Precision_",category)]<-precision
  
  #Will have NA if precision_boot or recall_boot are missing
  #This is most likely to happen when denominator is 0
  recall_df[,paste0("Recall_",category,"_SE")]<-sd(recall_boot,na.rm=TRUE)
  recall_df[,paste0("Precision_",category,"_SE")]<-sd(precision_boot,na.rm=TRUE)
  
}
recall_by_group<-recall_df%>%dplyr::select(all_of(contains(c("Recall","Precision"))))%>%distinct()

recall_by_group<-recall_by_group%>%dplyr::select(all_of(contains(c("Recall","Precision"))))%>%distinct()

recall_df<-recall_by_group%>%dplyr::select(contains("Recall"))%>%mutate(id=1)%>%
  pivot_longer(cols = starts_with("Recall"),names_to = "Method",values_to = "Recall")%>%dplyr::select(-id)

recall_df$Method<-gsub("Recall_","",recall_df$Method)
recall_SE<-recall_df%>%filter(grepl("_SE",Method))
recall_SE$Method<-gsub("_SE","",recall_SE$Method)
colnames(recall_SE)[2]<-"Recall_SE"

recall_df<-recall_df%>%filter(!grepl("_SE",Method))
recall_df<-left_join(recall_df,recall_SE,by="Method")

precision_df<-recall_by_group%>%dplyr::select(contains("Precision"))%>%mutate(id=1)%>%
  pivot_longer(cols = starts_with("Precision"),names_to = "Method",values_to = "Precision")%>%dplyr::select(-id)

precision_df$Method<-gsub("Precision_","",precision_df$Method)

precision_SE<-precision_df%>%filter(grepl("_SE",Method))
precision_SE$Method<-gsub("_SE","",precision_SE$Method)
colnames(precision_SE)[2]<-"Precision_SE"

precision_df<-precision_df%>%filter(!grepl("_SE",Method))
precision_df<-left_join(precision_df,precision_SE,by="Method")
#stop("860")
recall_df<-inner_join(recall_df,precision_df,by="Method")

recall_df<-recall_df%>%
  mutate(Recall_CI_lower=Recall-Recall_SE,
         Recall_CI_upper=Recall+Recall_SE,
         Precision_CI_lower=Precision-Precision_SE,
         Precision_CI_upper=Precision+Precision_SE)

recall_df<-recall_df%>%
  mutate(Recall_CI_lower=Recall-Recall_SE,
         Recall_CI_upper=Recall+Recall_SE,
         Precision_CI_lower=Precision-Precision_SE,
         Precision_CI_upper=Precision+Precision_SE)

load(paste0(data_wd,"/variants/IGVF_10M_variants_out_w_predictions.RData"))

IGVF_10M_variants_out_w_predictions$chr<-as.numeric(gsub("-.*","",IGVF_10M_variants_out_w_predictions$variant))

mean_table<-IGVF_10M_variants_out_w_predictions%>%
  summarize(across(starts_with("Predicted"), list(mean = mean), .names = "{.col}_{.fn}"))%>%pivot_longer(names_to="Annotation",values_to="Mean",cols=everything())

mean_table$Annotation<-gsub("_mean","",mean_table$Annotation)
mean_table$Annotation<-gsub("Predicted_Functional_","",mean_table$Annotation)
mean_table$Annotation[mean_table$Annotation=="Predicted_Functional"]<-"Overall"
mean_table<-mean_table%>%arrange(Annotation)
mean_table<-mean_table[c(4,1:3,5:nrow(mean_table)),]

format_number <- function(x) {
  ifelse(x < 0.01, formatC(x, format = "e", digits = 3), formatC(x, format = "f", digits = 3))
}


mean_table$Mean<-sapply(mean_table$Mean, format_number)
write_tsv(mean_table,file="~/OneDrive - Harvard University/PostDoc/Docs/annotate_variants/IGVF_10M_mean_anno_table.tsv")
#stop("889")
Harvard_MGH_variants_w_results_w_background<-bind_rows(Harvard_MGH_variants_w_results%>%filter((LDL_efflux_p_FDR_lt_.05|LDL_uptake_p_FDR_lt_.05)),IGVF_10M_variants_out_w_predictions%>%filter(!variant%in%Harvard_MGH_variants_w_results$variant)%>%mutate(LDL_efflux_p_FDR_lt_.05=FALSE,LDL_uptake_p_FDR_lt_.05=FALSE))

Harvard_MGH_variants_w_results_w_background$either_sig<-Harvard_MGH_variants_w_results_w_background$LDL_efflux_p_FDR_lt_.05|Harvard_MGH_variants_w_results_w_background$LDL_uptake_p_FDR_lt_.05

Harvard_MGH_variants_w_results$Predicted_Functional_Overall<-Harvard_MGH_variants_w_results$Predicted_Functional

Harvard_MGH_variants_w_results_w_background$Predicted_Functional_Overall<-Harvard_MGH_variants_w_results_w_background$Predicted_Functional

for(dataset in c("Harvard_MGH_variants_w_results"
                 ,"Harvard_MGH_variants_w_results_w_background")){
  for(cts_status in c("binary")){
    if(cts_status=="cts"){
      logistic_model<-glm(either_sig~aPC_Conservation_g.t._20+aPC_Protein_Function_g.t._20+aPC_Transcription_Factor_g.t._20+aPC_Epigenetics_Active_g.t._20,aPC_Epigenetics_Repressed_g.t._20+aPC_Transcription_Factor_g.t._20
                          ,data=get(dataset),family="binomial")
    }else{
      logistic_model<-glm(either_sig~Predicted_Functional_aPC+Predicted_Functional_chromBPnet+Predicted_Functional_MACIE+Predicted_Functional_cV2F+Predicted_Functional_liver_cV2F
                          +Predicted_Functional_liver_ASE
                          #+Predicted_Functional_ClinVar
                          +Predicted_Functional_liver_TLand
                          #+liver_caQTL
                          +Predicted_Functional_Overall
                          ,data=get(dataset),family="binomial")
    }
    
    coef_summary <- summary(logistic_model)$coefficients
    odds_ratios <- exp(coef_summary[, "Estimate"])
    ci_lower <- exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"])
    ci_upper <- exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
    
    # Step 2: Create a data frame for plotting
    plot_data <- data.frame(
      Term = rownames(coef_summary),
      Estimate = coef_summary[, "Estimate"],
      OddsRatio = odds_ratios,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    )
    plot_data$Term<-gsub("TRUE","",plot_data$Term)
    plot_data<-plot_data%>%filter(!Term=="(Intercept)")
    plot_data$Term<-gsub("Predicted_Functional_","",plot_data$Term)
    # Reorder terms for plotting
    plot_data$Term <- factor(plot_data$Term, levels = plot_data$Term)
    
    plot_data<-left_join(plot_data,recall_df,by=c("Term"="Method"))
    #stop("940")
    # Step 3: Create the plot using ggplot2
    p <- ggplot(plot_data, aes(x = Term, y = OddsRatio)) +
      geom_point() +
      geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
      coord_flip() +  # Flip coordinates to make it horizontal
      theme_cowplot() +
      geom_hline(yintercept =1,lty=2,lwd=1,color="red")+
      labs(x = "Annotation",
           y = "Odds Ratio for Significant Experimental Result (95% CI)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)
            ,title=element_text(hjust=.5,size=rel(1.5)),axis.text=element_text(size=rel(1.3))
            ,axis.title=element_text(size=rel(1)))
    
    plot_data$Method<-plot_data$Term
    p2 <- ggplot(plot_data, aes(x = OddsRatio, y = Recall,color=Method,shape=Method)) +
      geom_point(size=3) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.005) +
      geom_errorbar(aes(ymin = Recall_CI_lower, ymax = Recall_CI_upper), width = 0.005) +
      theme_cowplot() +
      labs(x = "Odds Ratio for Significant Experimental Result",
           y = "Recall",
           color="Method") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)
            ,title=element_text(hjust=.5,size=rel(1.5)),axis.text=element_text(size=rel(1.3))
            ,axis.title=element_text(size=rel(1)))+scale_color_viridis(discrete = TRUE)+scale_shape_manual(values = c(1,2,3,4,5,6,7,8))
    p2a <- ggplot(plot_data, aes(x = OddsRatio, y = Precision,color=Method,shape=Method)) +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.005) +
      geom_errorbar(aes(ymin = Precision_CI_lower, ymax = Precision_CI_upper), width = 0.005) +
      geom_point(size=3) +
      theme_cowplot() +
      labs(x = "Odds Ratio for Significant Experimental Result",
           y = "Precision",
           color="Method") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)
            ,title=element_text(hjust=.5,size=rel(1.5)),axis.text=element_text(size=rel(1.3))
            ,axis.title=element_text(size=rel(1)))+scale_color_viridis(discrete = TRUE)+scale_shape_manual(values = c(1, 2, 3, 4,5,6,7,8))
    
    p3<-ggplot(recall_df, aes(x = Recall, y = Precision, color = Method,shape=Method)) +
      geom_line(size = 1) +        # Line for each method
      geom_point(size = 3) +       # Points at each recall-precision pair
      geom_errorbarh(aes(xmin = Recall_CI_lower, xmax = Recall_CI_upper), height = 0.005) +
      geom_errorbar(aes(ymin = Precision_CI_lower, ymax = Precision_CI_upper), width = 0.005) +
      labs(
        title = "Precision-Recall Plot",
        x = "Recall",
        y = "Precision",
      )+scale_shape_manual(values = c(1, 2, 3, 4,5,6,7,8))+
      theme_cowplot()+scale_color_viridis(discrete = TRUE)+theme(title=element_text(hjust=.5,size=rel(1.5)))
    
    # Print the plot
    print(p)
    pdf_name<-dataset
    pdf(file=paste0(data_wd,"/variants/",pdf_name,"_logistic_model_",cts_status,".pdf"),width=10,height=6)
    print(p+ggtitle("Base-Editing Screen (784 significant variants)"))
    print(p2+ggtitle("Base-Editing Screen (784 significant variants)"))
    print(p2a+ggtitle("Base-Editing Screen (784 significant variants)"))
    print(p3+ggtitle("Base-Editing Screen (784 significant variants)"))
    dev.off()
  }
}

stop("1004")

logistic_model2<-glm((LDL_efflux_p_FDR_lt_.05|LDL_uptake_p_FDR_lt_.05)~Predicted_Functional_aPC+Predicted_Functional_chromBPnet+Predicted_Functional_MACIE+Predicted_Functional_cV2F+Predicted_Functional_liver_cV2F#+Predicted_Functional_liver_TLand
                     ,data=Harvard_MGH_variants_w_results_w_background)

# Save the plot as an image file
ggsave(paste0(data_wd,"/variants/baseedit_logistic_model_coefficients.png"), plot = p, width = 10, height = 6)


table((Harvard_MGH_variants_w_results$LDL_efflux_p_FDR_lt_.05|Harvard_MGH_variants_w_results$LDL_uptake_p_FDR_lt_.05),Harvard_MGH_variants_w_results$Predicted_Functional)

ggplot(data=temp2,aes(x=-1*log10(LDL_efflux_p_FDR),y=-1*log10(LDL_uptake_p_FDR)))+geom_hex()
#stop("790")

var1<-Harvard_MGH_variants_w_results%>%dplyr::select(chr,Position,ref.x,alt.x,variant,rsID,Predicted_Functional,contains("Predicted_Functional"),aPC_Conservation,aPC_Epigenetics_Active,aPC_Epigenetics_Repressed,aPC_Epigenetics_Transcription,aPC_Transcription_Factor,aPC_Protein_Function,MACIE_anyclass,Clinical_Significance,Allele_Specific_Evidence,liver_TLand,liver_cV2F,cV2F,LDL_efflux_p_FDR,LDL_uptake_p_FDR)
colnames(var1)[1:6]<-c("chr","position","ref","alt","variant","rsID")
write_csv(var1,file=paste0(data_wd,"/variants/Harvard_MGH_variants_all_predictions.csv"))

miami_plot_anno<-function(obj,file_name,max.overlaps=100
                          ,only_upper_plot=FALSE
                          ,color_var
                          ,split_var="class",plot.title
                          ,split_at_var="dELS",variants_always_label_lower="",
                          variants_always_label_upper=""
                          ,upper_ylab="",filtered_chrs=c(1,19)
                          ,lower_ylab="",data_col,genome_line=NA,suggestive_line=20
                          ,only_filtered_plots=FALSE,including_coding_variants=TRUE
                          ,label_genes_manual=FALSE,genes_to_label_manual=NULL
                          ,diff_y_limits=FALSE,suggestive_line_upper=NULL,suggestive_line_lower=NULL){
  #browser()
  pdf(file=file_name,width=10,height=6)
  if(only_filtered_plots==FALSE){
    plot_data_EVB<-prep_miami_data_nonp(data = obj, pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    #upper_labs<-""
    #lower_labs<-""
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    #browser()
    ggmiami_nonp(obj,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                 ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                 color_var=color_var,
                 chr_colors = viridis_pal()(length(unique(obj[,color_var]))),
                 pos="Position",max.overlaps = max.overlaps
                 ,title_name = plot.title,genome_line =genome_line
                 ,suggestive_line=suggestive_line
                 ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                 ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    #theme(plot.margin=unit(c(1,100,1.5,100),"pt"))
    
    plot_data_EVB<-prep_miami_data_nonp(data = obj%>%filter(!!sym(color_var)=="coding"), pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    #upper_labs<-""
    #lower_labs<-""
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    
    obj2<-obj%>%filter(!!sym(color_var)=="coding")
    obj2$Genecode_Comprehensive_Exonic_Category[is.na(obj2$Genecode_Comprehensive_Exonic_Category)]<-"unknown"
    ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                 ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                 color_var="Genecode_Comprehensive_Exonic_Category",
                 chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                 pos="Position",max.overlaps = max.overlaps
                 ,title_name = paste0(plot.title,": "," Coding Variants Only (Genecode_Comprehensive_Exonic_Category)"),genome_line =genome_line
                 ,suggestive_line=suggestive_line
                 ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                 ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    #browser()
    plot_data_EVB<-prep_miami_data_nonp(data = obj%>%filter(!is.na(CAGE_noncoding_category)), pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    
    obj2<-obj%>%filter(!is.na(CAGE_noncoding_category))
    obj2$CAGE_noncoding_cat<-factor(obj2$CAGE_noncoding_category,levels=c("Not CAGE Prom. or Enh.","CAGE Enhancer","CAGE Promoter"))
    
    ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                 ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                 color_var="CAGE_noncoding_cat",
                 chr_colors=viridis_pal()(length(unique(obj2$CAGE_noncoding_cat))),
                 pos="Position",max.overlaps = max.overlaps
                 ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                 ,suggestive_line=suggestive_line
                 ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                 ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    
    ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                 ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                 color_var="ENCODE_element",
                 chr_colors=viridis_pal()(length(unique(obj2$CAGE_noncoding_cat))),
                 pos="Position",max.overlaps = max.overlaps
                 ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                 ,suggestive_line=suggestive_line
                 ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                 ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
  }
  
  #browser()
  if(!is.null(filtered_chrs)){
    #browser()
    for(filtered_chr in filtered_chrs){
      plot.title.filt<-paste0(plot.title," in Chr. ",filtered_chr)
      obj2<-obj%>%filter(chr%in%filtered_chr)
      if(length(unique(obj2$Position))<2){
        obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
      }else{
        obj2<-obj2%>%mutate(Position_category=cut_interval(Position,2))
      }
      
      upper_variants_label<-obj2%>%filter(!!sym(split_var)==split_at_var)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      lower_variants_label<-obj2%>%filter(!!sym(split_var)!=split_at_var)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      if(including_coding_variants==TRUE){
        
        #Overall plot showing both coding and noncoding
        plot_data_EVB<-prep_miami_data_nonp(data = obj2, pos="Position",chr="chr"
                                            ,split_by = split_var, split_at = split_at_var, p = data_col)
        
        upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
        lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
        
        ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var
                     ,p=data_col
                     ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                     color_var=color_var,
                     #chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                     chr_colors = viridis_pal()(length(levels(obj2[,color_var]))),
                     pos="Position",max.overlaps = max.overlaps
                     ,title_name = paste0(plot.title)
                     ,genome_line =genome_line
                     ,suggestive_line=suggestive_line
                     ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                     ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual)
        
        ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var
                     ,p=data_col
                     ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                     color_var=color_var,
                     #chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                     chr_colors = viridis_pal()(length(levels(obj2[,color_var]))),
                     pos="Position",max.overlaps = max.overlaps
                     ,title_name = paste0(plot.title)
                     ,genome_line =genome_line
                     ,suggestive_line=suggestive_line
                     ,upper_labels_df = NULL,lower_labels_df = NULL
                     ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual)
        
        obj2<-obj%>%filter(chr%in%filtered_chr,!!sym(color_var)=="coding")
        if(nrow(obj2)>0){
          if(length(unique(obj2$Position))<2){
            obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
          }else{
            obj2<-obj2%>%mutate(Position_category=cut_interval(Position,2))
          }
          
          upper_variants_label<-obj2%>%filter(!!sym(split_var)==split_at_var)%>%
            group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
          
          lower_variants_label<-obj2%>%filter(!!sym(split_var)!=split_at_var)%>%
            group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
          
          #Just coding plot
          
          plot_data_EVB<-prep_miami_data_nonp(data = obj2, pos="Position",chr="chr"
                                              ,split_by = split_var, split_at = split_at_var, p = data_col)
          upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
          lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
          #browser()
          
          obj2$Genecode_Comprehensive_Exonic_Category[is.na(obj2$Genecode_Comprehensive_Exonic_Category)]<-"unknown"
          ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var="Genecode_Comprehensive_Exonic_Category",
                       chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = paste0(plot.title,": ","Coding Variants Only (Genecode_Comprehensive_Exonic_Category)")
                       ,genome_line =genome_line
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                       ,diff_y_limits=diff_y_lim,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
          ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var="Genecode_Comprehensive_Exonic_Category",
                       chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = paste0(plot.title,": ","Coding Variants Only (Genecode_Comprehensive_Exonic_Category)")
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = NULL,lower_labels_df = NULL
                       ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                       ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower) 
        }
      }
      
      obj2<-obj%>%filter(chr%in%filtered_chr,!is.na(CAGE_noncoding_category))
      if(length(unique(obj2$Position))<2){
        obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
      }else{
        obj2<-obj2%>%mutate(Position_category=cut_interval(Position,3))
      }
      #browser()
      upper_variants_label<-obj2%>%filter((!!sym(split_var))==split_at_var,value>suggestive_line)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      lower_variants_label<-obj2%>%filter((!!sym(split_var))!=split_at_var,value>suggestive_line)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      plot_data_EVB<-prep_miami_data_nonp(data = obj2%>%filter(!is.na(CAGE_noncoding_category)), pos="Position",chr="chr"
                                          ,split_by = split_var, split_at = split_at_var, p = data_col)
      
      upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
      lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
      
      obj2$CAGE_noncoding_cat<-factor(obj2$CAGE_noncoding_category,levels=c("Not CAGE Prom. or Enh.","CAGE Enhancer","CAGE Promoter"))
      #browser()
      
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="CAGE_noncoding_cat",
                   chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line
                   ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      ### Duplicate without Labels
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="CAGE_noncoding_cat",
                   chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line
                   ,upper_labels_df = NULL,lower_labels_df = NULL
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      # Lower Labs Only
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="CAGE_noncoding_cat",
                   chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line
                   ,upper_labels_df = NULL,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      
      #browser()
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="ENCODE_element",
                   chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = upper_labs,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      ##Duplicate without Labels
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="ENCODE_element",
                   chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = NULL,lower_labels_df = NULL
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      #Lower Labels only
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="ENCODE_element",
                   chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = NULL,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      ###############
      # Begin DHS category plots
      ###############    
      #browser()
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="DHS_category",
                   chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = upper_labs,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      ##Duplicate without Labels
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="DHS_category",
                   chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = NULL,lower_labels_df = NULL
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      #Lower Labels only
      ggmiami_nonp(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                   ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                   color_var="DHS_category",
                   chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
                   pos="Position",max.overlaps = max.overlaps
                   ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                   ,suggestive_line=suggestive_line,
                   upper_labels_df = NULL,lower_labels_df = lower_labs
                   ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                   ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      ######################    # end DHS plots
      
      
    }
  }
  
  dev.off()
}

miami_plot_anno_shape<-function(obj,file_name,max.overlaps=100
                                ,only_upper_plot=FALSE
                                ,color_var
                                ,shape_var
                                ,shape_var_levels=c("Not In ENCODE cCRE","ENCODE Enhancer","ENCODE Promoter")
                                ,split_var="class",plot.title
                                ,split_at_var="dELS",variants_always_label_lower="",
                                variants_always_label_upper=""
                                ,upper_ylab="",filtered_chrs=c(1,19)
                                ,lower_ylab="",data_col,genome_line=NA,suggestive_line=20
                                ,only_filtered_plots=FALSE,including_coding_variants=TRUE
                                ,label_genes_manual=FALSE,genes_to_label_manual=NULL
                                ,diff_y_limits=FALSE,suggestive_line_upper=NULL,suggestive_line_lower=NULL){
  #browser()
  pdf(file=file_name,width=11,height=6)
  if(only_filtered_plots==FALSE){
    plot_data_EVB<-prep_miami_data_nonp(data = obj, pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    #upper_labs<-""
    #lower_labs<-""
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    #browser()
    ggmiami_nonp_shape(obj,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var=color_var,shape_var=shape_var,
                       chr_colors = viridis_pal()(max(length(levels(obj[,color_var])))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = plot.title,genome_line =genome_line
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                       ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    #theme(plot.margin=unit(c(1,100,1.5,100),"pt"))
    
    plot_data_EVB<-prep_miami_data_nonp(data = obj%>%filter(!!sym(color_var)=="coding"), pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    #upper_labs<-""
    #lower_labs<-""
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    
    obj2<-obj%>%filter(!!sym(color_var)=="coding")
    obj2$Genecode_Comprehensive_Exonic_Category[is.na(obj2$Genecode_Comprehensive_Exonic_Category)]<-"unknown"
    ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var="Genecode_Comprehensive_Exonic_Category",shape_var=shape_var,
                       chr_colors=viridis_pal()(max(2,length(levels(obj2$Genecode_Comprehensive_Exonic_Category)))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = paste0(plot.title,": "," Coding Variants Only (Genecode_Comprehensive_Exonic_Category)"),genome_line =genome_line
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                       ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    #browser()
    plot_data_EVB<-prep_miami_data_nonp(data = obj%>%filter(!is.na(CAGE_noncoding_category)), pos="Position",chr="chr"
                                        ,split_by = split_var, split_at = split_at_var, p = data_col)
    
    upper_labs<-plot_data_EVB$upper%>%filter(variant%in%variants_always_label_upper)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    lower_labs<-plot_data_EVB$lower%>%filter(variant%in%variants_always_label_lower)%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
    
    obj2<-obj%>%filter(!is.na(CAGE_noncoding_category))
    obj2$CAGE_noncoding_cat<-factor(obj2$CAGE_noncoding_category,levels=c("Not CAGE Prom. or Enh.","CAGE Enhancer","CAGE Promoter"))
    
    ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                       ,shape_var=shape_var
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var="CAGE_noncoding_cat",
                       chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                       ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)#+
    
    ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col,shape_var=shape_var
                       ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                       color_var="ENCODE_element",
                       chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                       pos="Position",max.overlaps = max.overlaps
                       ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                       ,suggestive_line=suggestive_line
                       ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                       ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
  }
  
  #browser()
  if(!is.null(filtered_chrs)){
    #browser()
    for(filtered_chr in filtered_chrs){
      plot.title.filt<-paste0(plot.title," in Chr. ",filtered_chr)
      obj2<-obj%>%filter(chr%in%filtered_chr)
      if(length(unique(obj2$Position))<2){
        obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
      }else{
        obj2<-obj2%>%mutate(Position_category=cut_interval(Position,2))
      }
      
      upper_variants_label<-obj2%>%filter(!!sym(split_var)==split_at_var)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      lower_variants_label<-obj2%>%filter(!!sym(split_var)!=split_at_var)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      if(including_coding_variants==TRUE){
        
        #Overall plot showing both coding and noncoding
        plot_data_EVB<-prep_miami_data_nonp(data = obj2, pos="Position",chr="chr"
                                            ,split_by = split_var, split_at = split_at_var, p = data_col)
        
        upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
        lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
        
        ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var
                           ,shape_var=shape_var
                           ,p=data_col
                           ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                           color_var=color_var,
                           #chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                           chr_colors = viridis_pal()(length(levels(obj2[,color_var]))),
                           pos="Position",max.overlaps = max.overlaps
                           ,title_name = paste0(plot.title)
                           ,genome_line =genome_line
                           ,suggestive_line=suggestive_line
                           ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                           ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual)
        
        ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var
                           ,shape_var=shape_var
                           ,p=data_col
                           ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                           color_var=color_var,
                           #chr_colors=viridis_pal()(length(unique(obj2$Genecode_Comprehensive_Exonic_Category))),
                           chr_colors = viridis_pal()(length(levels(obj2[,color_var]))),
                           pos="Position",max.overlaps = max.overlaps
                           ,title_name = paste0(plot.title)
                           ,genome_line =genome_line
                           ,suggestive_line=suggestive_line
                           ,upper_labels_df = NULL,lower_labels_df = NULL
                           ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual)
        
        obj2<-obj%>%filter(chr%in%filtered_chr,!!sym(color_var)=="coding")
        if(nrow(obj2)>0){
          if(length(unique(obj2$Position))<2){
            obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
          }else{
            obj2<-obj2%>%mutate(Position_category=cut_interval(Position,2))
          }
          
          upper_variants_label<-obj2%>%filter(!!sym(split_var)==split_at_var)%>%
            group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
          
          lower_variants_label<-obj2%>%filter(!!sym(split_var)!=split_at_var)%>%
            group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
          
          #Just coding plot
          
          plot_data_EVB<-prep_miami_data_nonp(data = obj2, pos="Position",chr="chr"
                                              ,split_by = split_var, split_at = split_at_var, p = data_col)
          upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
          lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
          #browser()
          
          obj2$Genecode_Comprehensive_Exonic_Category[is.na(obj2$Genecode_Comprehensive_Exonic_Category)]<-"unknown"
          ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col,shape_var=shape_var
                             ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",
                             color_var="Genecode_Comprehensive_Exonic_Category",
                             chr_colors=viridis_pal()(length(levels(obj2$Genecode_Comprehensive_Exonic_Category))),
                             pos="Position",max.overlaps = max.overlaps
                             ,title_name = paste0(plot.title,": ","Coding Variants Only (Genecode_Comprehensive_Exonic_Category)")
                             ,genome_line =genome_line
                             ,suggestive_line=suggestive_line
                             ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                             ,diff_y_limits=diff_y_lim,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
          ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                             ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                             color_var="Genecode_Comprehensive_Exonic_Category",
                             chr_colors=viridis_pal()(length(levels(obj2$Genecode_Comprehensive_Exonic_Category))),
                             pos="Position",max.overlaps = max.overlaps
                             ,title_name = paste0(plot.title,": ","Coding Variants Only (Genecode_Comprehensive_Exonic_Category)")
                             ,suggestive_line=suggestive_line
                             ,upper_labels_df = NULL,lower_labels_df = NULL
                             ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                             ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower) 
        }
      }
      
      obj2<-obj%>%filter(chr%in%filtered_chr,!is.na(CAGE_noncoding_category))
      if(length(unique(obj2$Position))<2){
        obj2<-obj2%>%mutate(Position_category=factor(unique(obj2$Position)))
      }else{
        obj2<-obj2%>%mutate(Position_category=cut_interval(Position,3))
      }
      #browser()
      upper_variants_label<-obj2%>%filter((!!sym(split_var))==split_at_var,value>suggestive_line)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      lower_variants_label<-obj2%>%filter((!!sym(split_var))!=split_at_var,value>suggestive_line)%>%
        group_by(Position_category)%>%arrange(desc(value))%>%dplyr::slice(1:5)%>%pull(variant)
      
      plot_data_EVB<-prep_miami_data_nonp(data = obj2%>%filter(!is.na(CAGE_noncoding_category)), pos="Position",chr="chr"
                                          ,split_by = split_var, split_at = split_at_var, p = data_col)
      
      upper_labs<-plot_data_EVB$upper%>%filter(variant%in%c(variants_always_label_upper,upper_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
      lower_labs<-plot_data_EVB$lower%>%filter(variant%in%c(variants_always_label_lower,lower_variants_label))%>%group_by(chr)%>%arrange(!!sym(data_col))%>%ungroup()%>%mutate(label=paste0(variant,"\nAlt. AF: ",round(100*TOPMed_Bravo_AF,2),"%","\nClosest Gene: ",nearest_gene_gencode))%>%dplyr::select(rel_pos,value,label)
      
      obj2$CAGE_noncoding_cat<-factor(obj2$CAGE_noncoding_category,levels=c("Not CAGE Prom. or Enh.","CAGE Enhancer","CAGE Promoter"))
      #browser()
      
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="CAGE_noncoding_cat",
                         chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line
                         ,upper_labels_df = upper_labs,lower_labels_df = lower_labs
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      ### Duplicate without Labels
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="CAGE_noncoding_cat",
                         chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line
                         ,upper_labels_df = NULL,lower_labels_df = NULL
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      # Lower Labs Only
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="CAGE_noncoding_cat",
                         chr_colors=viridis_pal()(length(levels(obj2$CAGE_noncoding_cat))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line
                         ,upper_labels_df = NULL,lower_labels_df = lower_labs
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      
      #browser()
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="ENCODE_element",
                         chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line,
                         upper_labels_df = upper_labs,lower_labels_df = lower_labs
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      ##Duplicate without Labels
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="ENCODE_element",
                         chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line,
                         upper_labels_df = NULL,lower_labels_df = NULL
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      #Lower Labels only
      ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
                         ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
                         color_var="ENCODE_element",
                         chr_colors=viridis_pal()(length(levels(obj2$ENCODE_element))),
                         pos="Position",max.overlaps = max.overlaps
                         ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
                         ,suggestive_line=suggestive_line,
                         upper_labels_df = NULL,lower_labels_df = lower_labs
                         ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
                         ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      ###############
      # Begin DHS category plots
      ###############    
      #browser()
      # ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
      #              ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
      #              color_var="DHS_category",
      #              chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
      #              pos="Position",max.overlaps = max.overlaps
      #              ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
      #              ,suggestive_line=suggestive_line,
      #              upper_labels_df = upper_labs,lower_labels_df = lower_labs
      #              ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
      #              ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      # ##Duplicate without Labels
      # ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
      #              ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
      #              color_var="DHS_category",
      #              chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
      #              pos="Position",max.overlaps = max.overlaps
      #              ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
      #              ,suggestive_line=suggestive_line,
      #              upper_labels_df = NULL,lower_labels_df = NULL
      #              ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
      #              ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      # #Lower Labels only
      # ggmiami_nonp_shape(obj2,only_upper_plot=only_upper_plot,split_by=split_var,split_at=split_at_var,p=data_col
      #              ,upper_ylab=upper_ylab,lower_ylab = lower_ylab,chr="chr",shape_var=shape_var,
      #              color_var="DHS_category",
      #              chr_colors=viridis_pal()(length(levels(obj2$DHS_category))),
      #              pos="Position",max.overlaps = max.overlaps
      #              ,title_name = paste0(plot.title,": "," Non-Coding Variants Only"),genome_line =genome_line
      #              ,suggestive_line=suggestive_line,
      #              upper_labels_df = NULL,lower_labels_df = lower_labs
      #              ,genes_to_label_manual = genes_to_label_manual,label_genes_manual = label_genes_manual
      #              ,diff_y_limits=diff_y_limits,suggestive_line_upper=suggestive_line_upper,suggestive_line_lower=suggestive_line_lower)
      
      ######################    # end DHS plots
      
      
    }
  }
  
  dev.off()
}
#stop("1170")
Harvard_MGH_variants_w_results$ENCODE_liver_cCRE<-as.character(Harvard_MGH_variants_w_results$in_ENCODE_cts_liver_cCRE)
Harvard_MGH_variants_w_results$ENCODE_liver_cCRE[Harvard_MGH_variants_w_results$ENCODE_liver_cCRE=="TRUE"]<-"ENCODE Hepatocyte cCRE"

Harvard_MGH_variants_w_results$ENCODE_liver_cCRE[Harvard_MGH_variants_w_results$ENCODE_liver_cCRE=="FALSE" & Harvard_MGH_variants_w_results$ENCODE_element!="Not In ENCODE cCRE"]<-"ENCODE Non-Tissue-Specific cCRE"

Harvard_MGH_variants_w_results$ENCODE_liver_cCRE[Harvard_MGH_variants_w_results$ENCODE_liver_cCRE=="FALSE"]<-"Not In ENCODE cCRE"

Harvard_MGH_variants_w_results$ENCODE_liver_cCRE<-factor(Harvard_MGH_variants_w_results$ENCODE_liver_cCRE,levels=c("Not In ENCODE cCRE","ENCODE Non-Tissue-Specific cCRE","ENCODE Hepatocyte cCRE"))

Harvard_MGH_variants_w_results$aPC_Conservation_g.t._20_fact<-factor(Harvard_MGH_variants_w_results$aPC_Conservation_g.t._20)
sig_pred_fun_var<-Harvard_MGH_variants_w_results%>%filter((LDL_efflux_p_FDR_lt_.05|LDL_uptake_p_FDR_lt_.05)&Predicted_Functional&variant_category=="noncoding")


#### Miami plots w/ results
var_list<-list(c("neg_log10_LDL_efflux_p_FDR","neg_log10_LDL_uptake_p_FDR"))
for(j in 1:length(var_list)){
  vars<-unlist(var_list[j])
  top_var<-vars[1]
  bottom_var<-vars[2]
  
  if(top_var=="neg_log10_LDL_efflux_p_FDR"){
    top_lab="-log10(Efflux p-value)"
  }
  if(bottom_var=="neg_log10_LDL_uptake_p_FDR"){
    bottom_lab="-log10(Uptake p-value)"
  }
  
  if(sum(grepl("cV2F",vars)>0) & sum(grepl("aPC",vars)>0)
     | ("DNase"%in%vars) | "H3K27ac"%in%vars| "abs_chromBPnet_logfc"%in%vars){
    diff_y_lim=TRUE
  }else{
    diff_y_lim=FALSE
  }
  
  if(top_var%in%c("DNase","H3K27ac","H3K4me3","abs_chromBPnet_logfc")){
    suggest_line_upper=NA
  }else{
    if(grepl("cV2F",top_var)
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
    ){
      suggest_line_upper=-10*log10(1-.75+.0001)
    }else{
      # For now, number for -log10(.05)
      suggest_line_upper=-1*log10(.05)
    }
    
  }
  if(bottom_var%in%c("DNase","H3K27ac","H3K4me3")){
    suggest_line_lower=NA
  }else{
    if(bottom_var%in%c("cV2F","Liver.cV2F")
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
    ){
      suggest_line_lower=-10*log10(1-.75+.0001)
    }else{
      suggest_line_lower=-1*log10(.05)
    }
    
  }
  
  color_var<-"Predicted_Functional";color_var_levels<-rev(c("TRUE","FALSE"))
  #shape_var<-"ENCODE_element";shape_var_levels<-c("Not In ENCODE cCRE","ENCODE Enhancer","ENCODE Promoter")
  shape_var<-"ENCODE_liver_cCRE";shape_var_levels<-c("Not In ENCODE cCRE","ENCODE Non-Tissue-Specific cCRE","ENCODE Hepatocyte cCRE")
  #shape_var<-"aPC_Conservation_g.t._20";shape_var_levels<-c("FALSE","TRUE")
  
  Harvard_MGH_subset<-Harvard_MGH_variants_w_results%>%dplyr::select(variant,variant_category,!!sym(color_var),!!sym(shape_var),Genecode_Comprehensive_Exonic_Category,nearest_gene_gencode,CAGE_noncoding_category,genes_within_500000,rsID,ENCODE_element,DHS_category,chr,Position,!!sym(top_var),!!sym(bottom_var),TOPMed_Bravo_AF,Predicted_Functional)
  Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9,-10)]<-apply(Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9,-10)],2,FUN=as.numeric)
  
  #Need to drop indels, which will have missing funcitonal predictions
  Harvard_MGH_subset<-Harvard_MGH_subset%>%filter(!is.na(Predicted_Functional))
  
  Harvard_MGH_subset_long<-as.data.frame(pivot_longer(Harvard_MGH_subset,cols=c(all_of(top_var),all_of(bottom_var))
                                                      ,names_to="score"))
  Harvard_MGH_subset_long$color_var<-color_var
  Harvard_MGH_subset_long[,color_var]<-factor(Harvard_MGH_subset_long[,color_var],levels=color_var_levels)
  
  Harvard_MGH_subset_long$shape_var<-shape_var
  Harvard_MGH_subset_long[,shape_var]<-factor(Harvard_MGH_subset_long[,shape_var],levels=shape_var_levels)
  
  Harvard_MGH_subset_long$variant_category2<-factor(Harvard_MGH_subset_long$variant_category,levels=c("noncoding","coding"))
  
  miami_plot_anno_shape(Harvard_MGH_subset_long,only_upper_plot = FALSE,filtered_chrs = NULL
                  ,color_var="Predicted_Functional"
                  ,plot.title="Harvard_MGH Variants",upper_ylab=top_lab,lower_ylab=bottom_lab,shape_var=shape_var
                  ,file_name=paste0(plot_wd,"Harvard_MGH_miami_",paste(vars,collapse="_"),"_v1.pdf")
                  ,split_var="score",split_at_var = top_var,data_col="value",suggestive_line_upper = suggest_line_upper
                  ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
  miami_plot_anno_shape(Harvard_MGH_subset_long,only_upper_plot = FALSE,filtered_chrs = NULL
                        ,color_var="Predicted_Functional"
                        ,plot.title="Harvard_MGH Variants",upper_ylab=top_lab,lower_ylab=bottom_lab,shape_var="variant_category2",shape_var_levels=c("noncoding","coding")
                        ,file_name=paste0(plot_wd,"Harvard_MGH_miami2_",paste(vars,collapse="_"),"_v1.pdf")
                        ,split_var="score",split_at_var = top_var,data_col="value",suggestive_line_upper = suggest_line_upper
                        ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
  chr_list<-c(1,1,1,10,11,11,12,12,12,14,17,19,19,2,2,20,5,6,7,9,9,9,1,1)
  gene_list<-c("LDLRAP1","ARID1A","DR1","NHLRC2","TMEM258","FADS2","MVK"
               ,"HNF1A","SCARB1","RAD51B","ZNF652","SMARCA4","LDLR","GPN1"
               ,"MRPL33","HNF4A","HMGCR","MYLIP","NPC1L1","DENND4C"
               ,"RPS6","ABCA1","GPN1","NUDC")
  for(jj in 1:length(gene_list)){
    g=gene_list[jj]
    chr_val<-chr_list[jj]
    Harvard_MGH_subset_long_g<-Harvard_MGH_subset_long%>%filter(grepl(g,genes_within_500000)&chr==chr_val)
    if(nrow(Harvard_MGH_subset_long_g)>0){
      miami_plot_anno_shape(Harvard_MGH_subset_long_g,only_upper_plot = FALSE,only_filtered_plots = TRUE,filtered_chrs = chr_val,shape_var=shape_var
                      ,color_var="Predicted_Functional"
                      ,plot.title=paste0("Variants w/in 500kb of ",g),upper_ylab=top_lab,lower_ylab=bottom_lab
                      ,file_name=paste0(plot_wd,"Harvard_MGH_miami_",paste(vars,collapse="_"),"_v1_",g,".pdf")
                      ,split_var="score",split_at_var = top_var,data_col="value"
                      ,including_coding_variants = TRUE,label_genes_manual=TRUE,genes_to_label_manual=g
                      ,suggestive_line_upper = suggest_line_upper
                      ,suggestive_line_lower = suggest_line_lower
                      ,diff_y_limits = diff_y_lim)
    }
  }
}
for(j in 1:length(var_list)){
  vars<-unlist(var_list[j])
  top_var<-vars[1]
  bottom_var<-vars[2]
  
  if(top_var=="neg_log10_LDL_efflux_p_FDR"){
    top_lab="-log10(Efflux p-value)"
  }
  if(bottom_var=="neg_log10_LDL_uptake_p_FDR"){
    bottom_lab="-log10(Uptake p-value)"
  }
  
  if(sum(grepl("cV2F",vars)>0) & sum(grepl("aPC",vars)>0)
     | ("DNase"%in%vars) | "H3K27ac"%in%vars| "abs_chromBPnet_logfc"%in%vars){
    diff_y_lim=TRUE
  }else{
    diff_y_lim=FALSE
  }
  
  if(top_var%in%c("DNase","H3K27ac","H3K4me3","abs_chromBPnet_logfc")){
    suggest_line_upper=NA
  }else{
    if(grepl("cV2F",top_var)
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
    ){
      suggest_line_upper=-10*log10(1-.75+.0001)
    }else{
      # For now, number for -log10(.05)
      suggest_line_upper=-1*log10(.05)
    }
    
  }
  if(bottom_var%in%c("DNase","H3K27ac","H3K4me3")){
    suggest_line_lower=NA
  }else{
    if(bottom_var%in%c("cV2F","Liver.cV2F")
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
    ){
      suggest_line_lower=-10*log10(1-.75+.0001)
    }else{
      suggest_line_lower=-1*log10(.05)
    }
    
  }
  
  color_var<-"Predicted_Functional";color_var_levels<-rev(c("TRUE","FALSE"))
  #shape_var<-"ENCODE_element";shape_var_levels<-c("Not In ENCODE cCRE","ENCODE Enhancer","ENCODE Promoter")
  shape_var<-"ENCODE_liver_cCRE";shape_var_levels<-c("Not In ENCODE cCRE","ENCODE Non-Tissue-Specific cCRE","ENCODE Hepatocyte cCRE")
  #shape_var<-"aPC_Conservation_g.t._20";shape_var_levels<-c("FALSE","TRUE")
  
  Harvard_MGH_subset<-sig_pred_fun_var%>%dplyr::select(variant,variant_category,!!sym(color_var),!!sym(shape_var),Genecode_Comprehensive_Exonic_Category,nearest_gene_gencode,CAGE_noncoding_category,genes_within_500000,rsID,ENCODE_element,DHS_category,chr,Position,!!sym(top_var),!!sym(bottom_var),TOPMed_Bravo_AF,Predicted_Functional)
  Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9,-10)]<-apply(Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9,-10)],2,FUN=as.numeric)
  
  #Need to drop indels, which will have missing funcitonal predictions
  Harvard_MGH_subset<-Harvard_MGH_subset%>%filter(!is.na(Predicted_Functional))
  
  Harvard_MGH_subset_long<-as.data.frame(pivot_longer(Harvard_MGH_subset,cols=c(all_of(top_var),all_of(bottom_var))
                                                      ,names_to="score"))
  Harvard_MGH_subset_long$color_var<-color_var
  Harvard_MGH_subset_long[,color_var]<-factor(Harvard_MGH_subset_long[,color_var],levels=color_var_levels)
  
  Harvard_MGH_subset_long$shape_var<-shape_var
  Harvard_MGH_subset_long[,shape_var]<-factor(Harvard_MGH_subset_long[,shape_var],levels=shape_var_levels)
  
  Harvard_MGH_subset_long$variant_category2<-factor(Harvard_MGH_subset_long$variant_category,levels=c("noncoding","coding"))
  
  miami_plot_anno_shape(Harvard_MGH_subset_long,only_upper_plot = FALSE,filtered_chrs = NULL
                        ,color_var="Predicted_Functional"
                        ,plot.title="Harvard_MGH Variants",upper_ylab=top_lab,lower_ylab=bottom_lab,shape_var=shape_var
                        ,file_name=paste0(plot_wd,"Harvard_MGH_miami_sig_pred_func",paste(vars,collapse="_"),"_v1.pdf")
                        ,split_var="score",split_at_var = top_var,data_col="value",suggestive_line_upper = suggest_line_upper
                        ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
  miami_plot_anno_shape(Harvard_MGH_subset_long,only_upper_plot = FALSE,filtered_chrs = NULL
                        ,color_var="Predicted_Functional"
                        ,plot.title="Harvard_MGH Variants",upper_ylab=top_lab,lower_ylab=bottom_lab,shape_var="variant_category2",shape_var_levels=c("noncoding","coding")
                        ,file_name=paste0(plot_wd,"Harvard_MGH_miami2_sig_pred_func",paste(vars,collapse="_"),"_v1.pdf")
                        ,split_var="score",split_at_var = top_var,data_col="value",suggestive_line_upper = suggest_line_upper
                        ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
  chr_list<-c(1,1,1,10,11,11,12,12,12,14,17,19,19,2,2,20,5,6,7,9,9,9,1,1)
  gene_list<-c("LDLRAP1","ARID1A","DR1","NHLRC2","TMEM258","FADS2","MVK"
               ,"HNF1A","SCARB1","RAD51B","ZNF652","SMARCA4","LDLR","GPN1"
               ,"MRPL33","HNF4A","HMGCR","MYLIP","NPC1L1","DENND4C"
               ,"RPS6","ABCA1","GPN1","NUDC")
  for(jj in 1:length(gene_list)){
    g=gene_list[jj]
    chr_val<-chr_list[jj]
    Harvard_MGH_subset_long_g<-Harvard_MGH_subset_long%>%filter(grepl(g,genes_within_500000)&chr==chr_val)
    if(nrow(Harvard_MGH_subset_long_g)>0){
      miami_plot_anno_shape(Harvard_MGH_subset_long_g,only_upper_plot = FALSE,only_filtered_plots = TRUE,filtered_chrs = chr_val,shape_var=shape_var
                            ,color_var="Predicted_Functional"
                            ,plot.title=paste0("Variants w/in 500kb of ",g),upper_ylab=top_lab,lower_ylab=bottom_lab
                            ,file_name=paste0(plot_wd,"Harvard_MGH_miami_sig_pred_func",paste(vars,collapse="_"),"_v1_",g,".pdf")
                            ,split_var="score",split_at_var = top_var,data_col="value"
                            ,including_coding_variants = TRUE,label_genes_manual=TRUE,genes_to_label_manual=g
                            ,suggestive_line_upper = suggest_line_upper
                            ,suggestive_line_lower = suggest_line_lower
                            ,diff_y_limits = diff_y_lim)
    }
  }
}

###End results plots

ggplot(data=Harvard_MGH_variants_w_results,aes(x=-1*log10(LDL_efflux_p_FDR),y=-1*log10(LDL_uptake_p_FDR)))+geom_point()

var_list<-list(
  c(#"abs_chromBPnet_logfc",
    "aPC_Transcription_Factor","aPC_Transcription_Factor","aPC_Epigenetics_Transcription")
  ,c("aPC_Transcription_Factor","aPC_Conservation")
  ,c("aPC_Protein_Function","aPC_Conservation")
  ,c("aPC_Epigenetics_Active","aPC_Conservation")
  ,c("aPC_Epigenetics_Active","aPC_Epigenetics_Repressed")
  ,c("MACIE_protein","MACIE_conserved")
  ,c("MACIE_conserved","MACIE_regulatory")
  ,c("MACIE_anyclass","aPC_Transcription_Factor")
  #,c("cV2F","Liver.cV2F")
  ,c("DNase","H3K27ac")
  ,c("DNase","H3K4me3"))
#var_list<-list(c("MACIE_anyclass","MACIE_regulatory"))
var_list<-var_list[1]

for(j in 1:length(var_list)){
  vars<-unlist(var_list[j])
  top_var<-vars[1]
  bottom_var<-vars[2]
  
  if(sum(grepl("cV2F",vars)>0) & sum(grepl("aPC",vars)>0)
     | ("DNase"%in%vars) | "H3K27ac"%in%vars| "abs_chromBPnet_logfc"%in%vars){
    diff_y_lim=TRUE
  }else{
    diff_y_lim=FALSE
  }
  
  if(top_var%in%c("DNase","H3K27ac","H3K4me3","abs_chromBPnet_logfc")){
    suggest_line_upper=NA
  }else{
    if(grepl("cV2F",top_var)
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
       ){
      suggest_line_upper=-10*log10(1-.75+.0001)
    }else{
      suggest_line_upper=20
    }
    
  }
  if(bottom_var%in%c("DNase","H3K27ac","H3K4me3")){
    suggest_line_lower=NA
  }else{
    if(bottom_var%in%c("cV2F","Liver.cV2F")
       #c("MACIE_regulatory","MACIE_protein","MACIE_conserved","MACIE_anyclass")
       ){
      suggest_line_lower=-10*log10(1-.75+.0001)
    }else{
      suggest_line_lower=20
    }
    
  }
  
  color_var<-"variant_category";color_var_levels<-rev(c("coding","noncoding"))
  
  Harvard_MGH_subset<-Harvard_MGH_variants%>%dplyr::select(variant,!!sym(color_var),Genecode_Comprehensive_Exonic_Category,nearest_gene_gencode,CAGE_noncoding_category,genes_within_500000,rsID,ENCODE_element,DHS_category,chr,Position,!!sym(top_var),!!sym(bottom_var),TOPMed_Bravo_AF)
  Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9)]<-apply(Harvard_MGH_subset[,c(-1,-2,-3,-4,-5,-6,-7,-8,-9)],2,FUN=as.numeric)
  Harvard_MGH_subset_long<-as.data.frame(pivot_longer(Harvard_MGH_subset,cols=c(all_of(top_var),all_of(bottom_var))
                                              ,names_to="score"))
  Harvard_MGH_subset_long$color_var<-color_var
  Harvard_MGH_subset_long[,color_var]<-factor(Harvard_MGH_subset_long[,color_var],levels=color_var_levels)
  
  miami_plot_anno(Harvard_MGH_subset_long,only_upper_plot = FALSE,filtered_chrs = c(1,11,19)
                  ,color_var="variant_category"
                  ,plot.title="Harvard_MGH Variants",upper_ylab=top_var,lower_ylab=bottom_var
                  ,file_name=paste0(plot_wd,"Harvard_MGH_miami_",paste(vars,collapse="_"),"_v1.pdf")
                  ,split_var="score",split_at_var = top_var,data_col="value",suggestive_line_upper = suggest_line_upper
                  ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
  chr_list<-c(1,1,1,10,11,11,12,12,12,14,17,19,19,2,2,20,5,6,7,9,9,9,1,1)
  gene_list<-c("LDLRAP1","ARID1A","DR1","NHLRC2","TMEM258","FADS2","MVK"
               ,"HNF1A","SCARB1","RAD51B","ZNF652","SMARCA4","LDLR","GPN1"
               ,"MRPL33","HNF4A","HMGCR","MYLIP","NPC1L1","DENND4C"
               ,"RPS6","ABCA1","GPN1","NUDC")
  for(jj in 1:length(gene_list)){
    g=gene_list[jj]
    chr_val<-chr_list[jj]
    Harvard_MGH_subset_long_g<-Harvard_MGH_subset_long%>%filter(grepl(g,genes_within_500000)&chr==chr_val)
    if(nrow(Harvard_MGH_subset_long_g)>0){
      miami_plot_anno(Harvard_MGH_subset_long_g,only_upper_plot = FALSE,only_filtered_plots = TRUE,filtered_chrs = chr_val
                      ,color_var="variant_category"
                      ,plot.title=paste0("Harvard_MGH Variants w/in 500kb of ",g),upper_ylab=top_var,lower_ylab=bottom_var
                      ,file_name=paste0(plot_wd,"Harvard_MGH_miami_",paste(vars,collapse="_"),"_v1_",g,".pdf")
                      ,split_var="score",split_at_var = top_var,data_col="value"
                      ,including_coding_variants = TRUE,label_genes_manual=TRUE,genes_to_label_manual=g
                      ,suggestive_line_upper = suggest_line_upper
                      ,suggestive_line_lower = suggest_line_lower,diff_y_limits = diff_y_lim)
    }
  }
}

# Make some basic tablesin
#category vars at 1 should match mean_grouping_var?
create_summary_tables<-function(obj=Harvard_MGH_variants
                                ,mean_vars=NULL
                                ,mean_grouping_var=NULL
                                ,sum_vars
                                ,category_vars=c("TOPMed_AF_category","aPC_Protein_Function_category"),grouping_vars,first_col_name="Variant_Category",table_fn="mean"
                                ,table_title="Values by Variant Category"){
  #browser()
  full_table<-tibble()
for(grouping_var in grouping_vars){
    #browser()
    lol<-obj%>%filter(!is.na(!!sym(grouping_var)))%>%group_by(!!sym(grouping_var))%>%
      #filter(!is.na(!!sym(grouping_var)))%>%
      mutate("N"=n())
    if(nrow(lol)>0){
      if(!is.null(mean_vars)){
        lol<-lol%>%mutate(across(all_of(mean_vars),.fns=list(mean=mean)
                                 ,na.rm=TRUE,.names="{fn}_{col}"))
        if(!is.null(mean_grouping_var)){
          lol_temp<-lol%>%group_by(!!sym(grouping_var),!!sym(mean_grouping_var))%>%
            summarise(across(all_of(mean_vars),.fns=list(mean=mean)
                          ,na.rm=TRUE,.names="{fn}_{col}"))%>%ungroup()
            
            lol_temp2<-lol_temp%>%pivot_wider(values_from=all_of(paste0("mean_",mean_vars)),names_from=!!sym(mean_grouping_var))
            colnames(lol_temp2)<-gsub("_TRUE",paste0("_with_",mean_grouping_var,"_TRUE"),colnames(lol_temp2))
            colnames(lol_temp2)<-gsub("_FALSE",paste0("_with_",mean_grouping_var,"_FALSE"),colnames(lol_temp2))
            lol<-left_join(lol,lol_temp2,by=grouping_var)
        }
        
      }
      if(!is.null(sum_vars)){
        lol<-lol%>%mutate(across(all_of(sum_vars),.fns=list(sum=sum)
                                 ,na.rm=TRUE,.names="{fn}_{col}"))
      }
      lol<-lol%>%ungroup()
      #browser()
      #lol2<-lol
      #lol<-lol2
      
      lol<-dummy_cols(lol,select_columns = category_vars)%>%
        dplyr::select(!!sym(grouping_var),N,contains("mean")
                      ,contains("sum"),contains(eval(category_vars)))%>%
        dplyr::select(-all_of(category_vars))%>%group_by(!!sym(grouping_var))%>%
        mutate(across(all_of(contains(category_vars)),.fns=list(sum=sum)
                      ,na.rm=TRUE,.names="# with_{col}"))%>%ungroup()
      
      lol<-lol%>%dplyr::select(-contains(paste0("# with_mean_",mean_vars,"_with_",category_vars[1])))
      
      table_cols<-c(colnames(lol)[grepl("mean|#|sum",colnames(lol))])
      lol<-lol%>%arrange(desc(N))%>%group_by(!!sym(grouping_var),)%>%
        mutate(across(all_of(table_cols),round,2))%>%
        dplyr::select(!!sym(grouping_var),N,table_cols)%>%distinct()%>%ungroup()
      #browser()
      colnames(lol)<-gsub("# with_in","#_in",colnames(lol))
      colnames(lol)<-gsub("# with_Predicted","#_Predicted",colnames(lol))
      #
      #Order Columns
      #lol2<-lol
      
      
      #lol<-lol2
      lol<-lol%>%dplyr::select(all_of(grouping_var),N,
                               contains("#"),any_of(paste0("mean_",mean_vars))
                               #,contains(mean_vars)
                               ,contains(category_vars),everything())
      
      
      # lol<-lol%>%dplyr::select(all_of(grouping_var),N
      #                          #,contains(mean_vars)
      #                          ,contains(category_vars),everything())
      colnames(lol)<-gsub("_"," ",colnames(lol))
      colnames(lol)<-gsub("category ","",colnames(lol))
      colnames(lol)<-gsub("mean ","Mean\n",colnames(lol))
      colnames(lol)<-gsub(" ","\n",colnames(lol))
      
      #only first column name will differ, needs
      # to be the same to bind data frames together
      colnames(lol)[1]<-first_col_name
      
      full_table<-bind_rows(full_table,lol)
    }
  }
  
  #colnames(lol)[1]<-"Exonic\nCategory"
  
  #colnames(lol)<-gsub("TOPMed_","",colnames(lol))
  
  #grep("([^ ]* ){2}(.*)|.*", "\\2", colnames(lol))
  #sub("([^ ]* ){2}(.*)|.*", "\n\\2", colnames(lol))
  full_table<-full_table%>%dplyr::select(-all_of(contains("FALSE")))
  Table <- ggtexttable(full_table, rows = NULL) %>% tab_add_title(text = table_title, face = "bold")
  #browser()
  print(Table)
}
include_coding_variants=TRUE
#stop("1775")
for(a in 1){
  pdf(file=paste0(data_wd,"/variants/Harvard_MGH_upset_results.pdf"),width=10,height=6)  
  zzz<-0
  for(cat in c("noncoding")){
    zzz<-zzz+1
    if(cat%in%c("noncoding","coding")){
      Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category==cat)
    }
    if(cat=="all"){
      Harvard_MGH_var<-Harvard_MGH_variants_w_results
    }
    
    colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
    df<-Harvard_MGH_var%>%dplyr::select(variant,Predicted_Functional_aPC,Predicted_Functional_MACIE,Predicted_Functional_ClinVar,Predicted_Functional_Any,Predicted_Functional_liver_ASE,Predicted_Functional_liver_TLand,Predicted_Functional_chromBPnet,liver_caQTL,Predicted_Functional_liver_cV2F,Predicted_Functional_cV2F,`in_HepG2_ATAC-seq_peak`,`in_CATlas_peak_Hepatocyte`,liver_caQTL,LDL_uptake_p_FDR_lt_.05,LDL_efflux_p_FDR_lt_.05,in_ENCODE_cts_adipose_cCRE,in_ENCODE_cts_liver_cCRE,in_ENCODE_cts_pancreas_cCRE,in_ENCODE_cts_small_intestine_cCRE,liver_TLand_top5pct,liver_caQTL,sig_pvalue_LDL_TOPMed_F8,sig_pvalue_LDL_UKB_200K)%>%dplyr::rename(Predicted_Functional_Any=Predicted_Functional_Any)%>%pivot_longer(-variant)
    
    if(cat%in%c("noncoding")){plot_title<-"Harvard_MGH Variants: Non-Coding Only"}
    if(cat=="coding"){plot_title<-"Harvard_MGH Variants: Coding Only"}
    if(cat=="all"){plot_title<-"Harvard_MGH Variants: All Variants"}
    
    df1<-df%>%group_by(variant)%>%summarize(output=paste(name[value],collapse=","))%>%ungroup()%>%filter(!output=="")
    
    all_cat<-as.list(strsplit(df1$output,","))
    df1$all_cats<-all_cat
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL",
            "Predicted_Functional_Any"
            ,"Predicted_Functional_liver_ASE"
            ,"Predicted_Functional_aPC"
            ,"Predicted_Functional_liver_cV2F"
            ,"Predicted_Functional_MACIE"
            ,"Predicted_Functional_chromBPnet"
            ,"Predicted_Functional_liver_TLand")}else{
              sets<-c("Predicted_Functional_Any"
                      ,"Predicted_Functional_liver_ASE"
                      ,"Predicted_Functional_aPC"
                      ,"Predicted_Functional_liver_cV2F"
                      ,"Predicted_Functional_MACIE"
                      ,"Predicted_Functional_chromBPnet"
                      ,"Predicted_Functional_liver_TLand")
            }
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL","Predicted_Functional_Any",
                           "Predicted_Functional_liver_ASE"
                           ,"Predicted_Functional_aPC"
                           ,"Predicted_Functional_liver_cV2F"
                           ,"Predicted_Functional_MACIE"
                           ,"Predicted_Functional_chromBPnet"
                           ,"in_HepG2_ATAC-seq_peak"
                           ,"in_CATlas_peak_Hepatocyte"
                           ,"Predicted_Functional_liver_TLand"
    )}else{
      sets<-c("Predicted_Functional_Any"
        ,"Predicted_Functional_liver_ASE"
        ,"Predicted_Functional_aPC"
        ,"Predicted_Functional_liver_cV2F"
        ,"Predicted_Functional_MACIE"
        ,"Predicted_Functional_chromBPnet"
        ,"in_HepG2_ATAC-seq_peak"
        ,"in_CATlas_peak_Hepatocyte"
        ,"Predicted_Functional_liver_TLand"
      )}
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL",
      "in_HepG2_ATAC-seq_peak"
      ,"Predicted_Functional_Any"
      ,"LDL_efflux_p_FDR_lt_.05"
      ,"LDL_uptake_p_FDR_lt_.05"
      ,"sig_pvalue_LDL_UKB_200K"
      ,"sig_pvalue_LDL_TOPMed_F8"
    )
    }else{
      sets<-c("in_HepG2_ATAC-seq_peak"
              ,"Predicted_Functional_Any"
              ,"LDL_efflux_p_FDR_lt_.05"
              ,"LDL_uptake_p_FDR_lt_.05"
              ,"sig_pvalue_LDL_UKB_200K")
    }
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    
   # grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL",
      "in_HepG2_ATAC-seq_peak"
      ,"Predicted_Functional_Any"
      ,"LDL_efflux_p_FDR_lt_.05"
      ,"LDL_uptake_p_FDR_lt_.05"
      ,"sig_pvalue_LDL_UKB_200K"
      ,"sig_pvalue_LDL_TOPMed_F8"
    )
    }else{
      sets<-c("in_HepG2_ATAC-seq_peak"
              ,"Predicted_Functional_Any"
              ,"GLGC_pvalue_lt_5e8"
              ,"sig_pvalue_LDL_UKB_200K")
    }
    
    print(df1 %>%filter(grepl("p_FDR_lt_.05",output))%>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=c(#"Predicted_Functional_Any",
                            "Predicted_Functional_liver_ASE"
                            ,"Predicted_Functional_aPC"
                            ,"Predicted_Functional_MACIE"
                            ,"Predicted_Functional_chromBPnet"
                            ,"Predicted_Functional_cV2F"
                            ,"Predicted_Functional_liver_TLand"
                            ,"Predicted_Functional_liver_cV2F"
                            ,"LDL_efflux_p_FDR_lt_.05"
                            ,"LDL_uptake_p_FDR_lt_.05")
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(paste0("Experimental Sig Variants:",plot_title),x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c("Predicted_Functional_Any"
                                 ,"Predicted_Functional_liver_ASE"
                                 ,"Predicted_Functional_aPC"
                                 ,"Predicted_Functional_MACIE"
                                 ,"Predicted_Functional_chromBPnet"
                                 ,"Predicted_Functional_liver_TLand"
                                 ,"Predicted_Functional_liver_cV2F"
                                 ,"sig_pvalue_LDL_TOPMed_F8"
                                 ,"sig_pvalue_LDL_UKB_200K")
    }else{
      sets<-c("Predicted_Functional_Any"
              ,"Predicted_Functional_liver_ASE"
              ,"Predicted_Functional_aPC"
              ,"Predicted_Functional_MACIE"
              ,"Predicted_Functional_chromBPnet"
              ,"Predicted_Functional_liver_TLand"
              ,"Predicted_Functional_liver_cV2F"
              ,"sig_pvalue_LDL_UKB_200K")
        }
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
   
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=c("Predicted_Functional_Any"
                                  ,"LDL_uptake_p_FDR_lt_.05"
                                  ,"LDL_efflux_p_FDR_lt_.05")
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=1, gp=gpar(fontsize=20))
    
    print(df1 %>%filter((grepl("LDL_efflux",output)|grepl("LDL_uptake",output))&grepl("Pred",output))
          %>%distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=c("Predicted_Functional_liver_ASE"
                                  ,"Predicted_Functional_aPC"
                                  ,"Predicted_Functional_MACIE"
                                  ,"Predicted_Functional_chromBPnet"
                                  ,"Predicted_Functional_liver_TLand"
                                  ,"Predicted_Functional_liver_cV2F"
                                  ,"LDL_uptake_p_FDR_lt_.05"
                                  ,"LDL_efflux_p_FDR_lt_.05")
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL",
                                            "in_HepG2_ATAC-seq_peak"
                                 ,"in_CATlas_peak_Hepatocyte"
                                 ,"in_ENCODE_cts_liver_cCRE"
                                 ,"Predicted_Functional_Any"
                                 #,"abs_chromBPnet_logfc_top10pct"
                                 ,"LDL_uptake_p_FDR_lt_.05"
                                 ,"LDL_efflux_p_FDR_lt_.05")
    }else{
      sets<-c("in_HepG2_ATAC-seq_peak"
              ,"in_CATlas_peak_Hepatocyte"
              ,"in_ENCODE_cts_liver_cCRE"
              ,"Predicted_Functional_Any"
              #,"abs_chromBPnet_logfc_top10pct"
              ,"LDL_uptake_p_FDR_lt_.05"
              ,"LDL_efflux_p_FDR_lt_.05")
    }
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
    if(cat%in%c("noncoding","all")){sets<-c(#"liver_caQTL",
                                            "in_HepG2_ATAC-seq_peak"
                                 ,"Predicted_Functional_Any"
                                 ,"LDL_uptake_p_FDR_lt_.05"
                                 ,"LDL_efflux_p_FDR_lt_.05"
                                 ,"sig_pvalue_LDL_UKB_200K"
                                 ,"sig_pvalue_LDL_TOPMed_F8"
                                 )
    }else{
      sets<-c("in_HepG2_ATAC-seq_peak"
              ,"Predicted_Functional_Any"
              ,"LDL_uptake_p_FDR_lt_.05"
              ,"LDL_efflux_p_FDR_lt_.05"
              ,"sig_pvalue_LDL_UKB_200K")
    }
    
    print(df1 %>%
            distinct(variant,.keep_all = TRUE) %>%
            unnest(cols = all_cats) %>% distinct()%>%
            mutate(ctmember=1) %>%
            pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
            as.data.frame() %>%
            UpSetR::upset(nsets=999
                          ,sets=sets
                          ,keep.order = FALSE
                          ,order.by='freq'
                          ,nintersects=NA))
    #grid.text(plot_title,x = .6, y=0.95, gp=gpar(fontsize=20))
    
  }
  
  
  # for(cts_var in c("max_aPC","aPC_Conservation","aPC_Epigenetics_Active","aPC_Epigenetics_Repressed","aPC_Transcription_Factor","MACIE_conserved","MACIE_anyclass","AS_PHRED_max","DNase_PHRED","abs_chromBPnet_logfc","chromBPnet_logfc_mu","LDL_efflux_p_FDR","LDL_uptake_p_FDR")){
  #   lol<-roc(response=Harvard_MGH_var$LDL_efflux_p_FDR_lt_.05,
  #            predictor=Harvard_MGH_var%>%pull(cts_var)
  #            ,direction="<")
  #   print(ggroc(lol)+ theme_minimal() + ggtitle(paste0("LDL Efflux: ",cts_var)) +
  #           geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed"))
  #   
  #   lol<-roc(response=Harvard_MGH_var$LDL_uptake_p_FDR_lt_.05,
  #            predictor=Harvard_MGH_var%>%pull(cts_var)
  #            ,direction="<")
  #   print(ggroc(lol)+ theme_minimal() + ggtitle(paste0("LDL Uptake: ",cts_var)) +
  #           geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed"))
  # }
  
  # if(include_coding_variants==TRUE){
  #   Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="coding")
  #   colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  #   df<-Harvard_MGH_var%>%dplyr::select(variant,Predicted_Functional_aPC,Predicted_Functional_MACIE,Predicted_Functional_ClinVar,Predicted_Functional_Any,Predicted_Functional_ASE,Predicted_Functional_liver_TLand,Predicted_Functional_chromBPnet,Predicted_Functional_liver_cV2F,Predicted_Functional_cV2F,abs_chromBPnet_logfc_top10pct,`in_HepG2_ATAC-seq_peak`,`in_CATlas_peak_Hepatocyte`,liver_caQTL,LDL_uptake_p_FDR_lt_.05,LDL_efflux_p_FDR_lt_.05,in_ENCODE_cts_adipose_cCRE,in_ENCODE_cts_liver_cCRE,in_ENCODE_cts_pancreas_cCRE,in_ENCODE_cts_small_intestine_cCRE,liver_TLand_top5pct,sig_pvalue_LDL_TOPMed_F8,sig_pvalue_LDL_UKB_200K)%>%dplyr::rename(Predicted_Functional_Any=Predicted_Functional_Any)%>%pivot_longer(-variant)
  #   df1<-df%>%group_by(variant)%>%summarize(output=paste(name[value],collapse=","))%>%ungroup()%>%filter(!output=="")
  #   
  #   all_cat<-as.list(strsplit(df1$output,","))
  #   df1$all_cats<-all_cat
  #   print(df1 %>%
  #           distinct(variant,.keep_all = TRUE) %>%
  #           unnest(cols = all_cats) %>% 
  #           mutate(ctmember=1) %>%
  #           pivot_wider(names_from = all_cats, values_from = ctmember, values_fill = list(ctmember = 0)) %>%
  #           as.data.frame() %>%
  #           UpSetR::upset(nsets=999,
  #                         keep.order = FALSE
  #                         ,order.by='freq'
  #                         ,nintersects=NA
  #                         ,sets=c("Predicted_Functional_Any"
  #                                ,"Predicted_Functional_ASE"
  #                                ,"Predicted_Functional_aPC"
  #                                ,"Predicted_Functional_MACIE"
  #                                ,"Predicted_Functional_chromBPnet"
  #                                ,"Predicted_Functional_liver_TLand"
  #                                ,"Predicted_Functional_liver_cV2F"
  #                                ,"sig_pvalue_LDL_TOPMed_F8"
  #                                ,"sig_pvalue_LDL_UKB_200K")))
  #   grid.text("Harvard_MGH Variants: Coding Only",x = .6, y=0.95, gp=gpar(fontsize=20))
  # }
  
  dev.off()
  
  pdf(file=paste0(data_wd,"/variants/Harvard_MGH_tables_v6.pdf"),width=16,height=16)
  
  Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="noncoding")%>%dplyr::select(-sum_pred_func)
  
  colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  
  Harvard_MGH_var$Significant_pvalue<-"Non Signif. p-value"
  Harvard_MGH_var$Significant_pvalue[(Harvard_MGH_var$LDL_efflux_p_FDR_lt_.05|Harvard_MGH_var$LDL_uptake_p_FDR_lt_.05)]<-"Signif. p-value"
  
  # create_summary_tables(obj=Harvard_MGH_var,
  #                       category_vars=c("Predicted_Functional_Any","in_HepG2_ATAC-seq_peak"
  #                                       ,"MACIE_anyclass_g.t._20","aPC_Transcription_Factor_g.t._20","Allele_Specific_Evidence"
  #                                       ,"abs_chromBPnet_logfc_top10pct","abs_chromBPnet_logfc_top1pct")
  #                       #,mean_vars=c("MACIE_anyclass_PHRED")
  #                       ,mean_vars=NULL
  #                       ,sum_vars=NULL
  #                       ,grouping_vars=c("variant_category","ENCODE_element"
  #                                        ,"CAGE_noncoding_category2",
  #                                        "DHS_category"
  #                                        ,"Drug_Response","Pathogenic","functional_category")
  #                       ,first_col_name = "Variant Category"
  #                       ,table_title = "Harvard_MGH Variants: Non-Coding Only")
  
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("Predicted_Functional_Any","Predicted_Functional_aPC","Predicted_Functional_MACIE","Predicted_Functional_ClinVar"
,"Predicted_Functional_ASE"                                        ,"Predicted_Functional_liver_ASE","Predicted_Functional_chromBPnet","Predicted_Functional_liver_TLand","Predicted_Functional_cV2F","Predicted_Functional_liver_cV2F"
                                        ,"in_HepG2_ATAC-seq_peak","liver_caQTL","caQTL_eQTL_coloc")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("functional_category"
                                         #"variant_category"
                                         ,"ENCODE_element","Significant_pvalue"
                                         #,"CAGE_noncoding_category2",
                                         #"DHS_category",
                                         #,"Drug_Response"
                                         #,"Pathogenic"
                        )
                        ,first_col_name = "Variant Category"
                        ,table_title = "Base Editing Variants: Non-Coding Only")
  
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("Predicted_Functional_Any")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_grouping_var="Predicted_Functional_Any"
                        ,mean_vars=c("DNase","H3K27ac","H3K4me3")
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","ENCODE_element"
                                         ,"CAGE_noncoding_category2",
                                         "DHS_category"
                                         ,"Drug_Response","Pathogenic")
                        ,first_col_name = "Variant Category"
                        ,table_title = "Harvard_MGH Variants: Non-Coding Only")
  
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("TOPMed_AF_category","Predicted_Functional_Any"),
                        mean_vars="TOPMed_Bravo_AF",sum_vars=NULL
                        ,table_title="Harvard_MGH Variants: Non-Coding Only"
                        ,grouping_vars=c("variant_category","ENCODE_element"
                                         ,"CAGE_noncoding_category2"
                                         ,"DHS_category"
                                         ,"Drug_Response","Pathogenic"))
  #https://docs.google.com/spreadsheets/d/1zcnkqh1n8dLQXvZjp6cW20vppB1IWAXGQpK6jKcIQRw/edit#gid=1629529355
  if(include_coding_variants==TRUE){
    Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="coding")%>%dplyr::select(-sum_pred_func)
    colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
    create_summary_tables(obj=Harvard_MGH_var,
                          category_vars=c("Predicted_Functional_Any"
                                          ,"in_HepG2_ATAC-seq_peak","MACIE_anyclass_g.t._20")
                          #,mean_vars=c("MACIE_anyclass_PHRED")
                          ,mean_vars=NULL
                          ,sum_vars=NULL
                          ,grouping_vars=c("variant_category"
                                           ,"Genecode_Comprehensive_Exonic_Category"
                                           ,"disruptive_missense"
                                           ,"Drug_Response","Pathogenic","functional_category")
                          ,first_col_name = "Variant Category"
                          ,table_title = "Harvard_MGH Variants: Coding Only")
  }
  chr_list<-c(1,1,1,10,11,11,12,12,12,14,17,19,19,2,2,20,5,6,7,9,9,9,1,1)
  gene_list<-c("LDLRAP1","ARID1A","DR1","NHLRC2","TMEM258","FADS2","MVK"
               ,"HNF1A","SCARB1","RAD51B","ZNF652","SMARCA4","LDLR","GPN1"
               ,"MRPL33","HNF4A","HMGCR","MYLIP","NPC1L1","DENND4C"
               ,"RPS6","ABCA1","GPN1","NUDC")
  #chr_list<-numeric()
  for(jj in 1:length(chr_list)){
    g<-gene_list[jj]
    chr_val<-chr_list[jj]
    Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter((chr==chr_val&grepl(g,genes_within_500000)
                                                    #chr==1&grepl("PCSK9",genes_within_500000)
    ),variant_category=="noncoding")
    colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
    if(nrow(Harvard_MGH_var)>0){
      create_summary_tables(obj=Harvard_MGH_var,
                            category_vars=c("Predicted_Functional_Any"
                                            ,"in_HepG2_ATAC-seq_peak","MACIE_anyclass_g.t._20")
                            #,mean_vars=c("MACIE_anyclass_PHRED")
                            ,mean_vars=NULL
                            ,sum_vars=NULL
                            ,grouping_vars=c("variant_category","ENCODE_element"
                                             ,"CAGE_noncoding_category2"
                                             ,"DHS_category"
                                             ,"Drug_Response","Pathogenic","functional_category")
                            ,first_col_name = "Variant Category"
                            ,table_title = paste0("Harvard_MGH Variants: Non-Coding Within 100kb of ",g))
      create_summary_tables(obj=Harvard_MGH_var,
                            category_vars=c("Predicted_Functional_Any")
                            #,mean_vars=c("MACIE_anyclass_PHRED")
                            ,mean_vars=c("DNase","H3K27ac","H3K4me3")
                            ,mean_grouping_var="Predicted_Functional_Any"
                            ,sum_vars=NULL
                            ,grouping_vars=c("variant_category","ENCODE_element"
                                             ,"CAGE_noncoding_category2",
                                             "DHS_category"
                                             ,"Drug_Response","Pathogenic")
                            ,first_col_name = "Variant Category"
                            ,table_title =  paste0("Harvard_MGH Variants: Non-Coding Within 100kb of ",g))
      if(include_coding_variants==TRUE){
        Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter((chr==chr_val&grepl(g,genes_within_500000)
                                                        #chr==1&grepl("PCSK9",genes_within_500000)
        ),
        variant_category=="coding")%>%dplyr::select(-sum_pred_func)
        colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
        if(nrow(Harvard_MGH_var)>0){
          create_summary_tables(obj=Harvard_MGH_var,
                                category_vars=c("Predicted_Functional_Any"
                                                ,"in_HepG2_ATAC-seq_peak","MACIE_anyclass_g.t._20")
                                #,mean_vars=c("MACIE_anyclass_PHRED")
                                ,mean_vars=NULL
                                ,sum_vars=NULL
                                ,grouping_vars=c("variant_category"
                                                 ,"Genecode_Comprehensive_Exonic_Category"
                                                 ,"disruptive_missense"
                                                 ,"Drug_Response","Pathogenic","functional_category")
                                ,first_col_name = "Variant Category"
                                ,table_title = paste0("Harvard_MGH Variants: Coding Within 100kb of ",g))
          
        }
      }
    }
    
  }
  
  Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="noncoding")%>%dplyr::select(-sum_pred_func)
  colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("MACIE_anyclass_category","Predicted_Functional_Any")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","ENCODE_element","CAGE_noncoding_category2")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Non-Coding Only")
  
  Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="noncoding")%>%dplyr::select(-sum_pred_func)
  colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("DNase_PHRED_category","H3K27ac_PHRED_category","H3K4me3_PHRED_category"),
                        mean_vars=c("DNase_PHRED","H3K27ac_PHRED","H3K4me3"),sum_vars=NULL
                        ,grouping_vars=c("variant_category","ENCODE_element","CAGE_noncoding_category2")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Non-Coding Only")
  
  Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="coding")%>%dplyr::select(-sum_pred_func)
  colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("MACIE_anyclass_category","Predicted_Functional_Any")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","Genecode_Comprehensive_Exonic_Category"
                                         ,"disruptive_missense"
                                         ,"Drug_Response","Pathogenic")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Coding Only")
  Harvard_MGH_var<-Harvard_MGH_variants_w_results%>%filter(variant_category=="coding")%>%dplyr::select(-sum_pred_func)
  colnames(Harvard_MGH_var)[colnames(Harvard_MGH_var)=="Predicted_Functional"]<-"Predicted_Functional_Any"
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("MACIE_protein_category","MACIE_conserved_category","MACIE_anyclass_category","Predicted_Functional_Any")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","Genecode_Comprehensive_Exonic_Category"
                                         ,"disruptive_missense"
                                         ,"Drug_Response","Pathogenic")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Coding Only")
  
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("Predicted_Functional_Any","MACIE_protein_g.t._20","MACIE_conserved_g.t._20","MACIE_anyclass_g.t._20")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","Genecode_Comprehensive_Exonic_Category"
                                         ,"disruptive_missense"
                                         ,"Drug_Response","Pathogenic")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Coding Only")
  
  create_summary_tables(obj=Harvard_MGH_var,
                        category_vars=c("Predicted_Functional_Any","MACIE_protein_g.t._20","MACIE_conserved_g.t._20","MACIE_anyclass_g.t._20")
                        #,mean_vars=c("MACIE_anyclass_PHRED")
                        ,mean_vars=NULL
                        ,sum_vars=NULL
                        ,grouping_vars=c("variant_category","Genecode_Comprehensive_Exonic_Category"
                                         ,"disruptive_missense"
                                         ,"Drug_Response","Pathogenic")
                        ,first_col_name = "Nearest Gene"
                        ,table_title = "Harvard_MGH Variants: Coding Only")
  
  
  # create_summary_tables(obj=Harvard_MGH_variants,
  #                       category_vars=c("PolyPhenCat"),
  #                       mean_vars=c("aPC_Conservation"),sum_vars=NULL
  #                       ,grouping_vars=c("Genecode_Comprehensive_Exonic_Category"
  #                                        ,"disruptive_missense"
  #                                        ,"Drug_Response","Pathogenic")
  #                       ,table_title = "Harvard_MGH Variants")
  dev.off()
}
stop("2620")

Harvard_MGH_variants_w_results$abs_LDL_efflux_z<-abs(Harvard_MGH_variants_w_results$LDL_efflux_mu_z_adj)

Harvard_MGH_variants_w_results$abs_LDL_uptake_z<-abs(Harvard_MGH_variants_w_results$LDL_uptake_mu_z_adj)

Harvard_MGH_variants_w_results_sigonly<-Harvard_MGH_variants_w_results%>%filter(LDL_efflux_p_FDR<.05 | LDL_uptake_p_FDR<.05)

Harvard_MGH_variants_w_results$LDL_efflux_z_sq<-Harvard_MGH_variants_w_results$LDL_efflux_mu_z_adj^2
lm1<-lm(LDL_efflux_z_sq~aPC_Conservation+aPC_Epigenetics_Active+abs_chromBPnet_logfc,data=Harvard_MGH_variants_w_results)
summary(lm1)

Harvard_MGH_variants_w_results$LDL_efflux_mu_z_adj_gt2<-abs(Harvard_MGH_variants_w_results$LDL_efflux_mu_z_adj)>1.96

Harvard_MGH_variants_w_results$LDL_uptake_mu_z_adj_gt2<-abs(Harvard_MGH_variants_w_results$LDL_uptake_mu_z_adj)>1.96

pdf(paste0(data_wd,"/variants/results_ind_plots.pdf"),width=10,height=6)
for(yvar in c("aPC_Conservation","aPC_Epigenetics_Active")){
  print(ggplot(data=Harvard_MGH_variants_w_results
         ,aes(x=LDL_efflux_mu_z_adj,y=!!sym(yvar)
              ,color=LDL_efflux_mu_z_adj_gt2
              ,shape=!!sym(paste0(yvar,"_g.t._20"))))+
          scale_color_viridis_d()+
    geom_point()+geom_hline(aes(yintercept=20,color="red"))+
    geom_vline(aes(xintercept=-1.96,color="red"))+
    geom_vline(aes(xintercept=1.96,color="red"))+
    xlab("LDL Efflux Z-score")+guides(color="none",shape="none")+theme_cowplot())
  
  print(ggplot(data=Harvard_MGH_variants_w_results
               ,aes(x=LDL_uptake_mu_z_adj,y=!!sym(yvar)
                    ,color=LDL_uptake_mu_z_adj_gt2
                    ,shape=!!sym(paste0(yvar,"_g.t._20"))))+
          scale_color_viridis_d()+
          geom_point()+geom_hline(aes(yintercept=20,color="red"))+
          geom_vline(aes(xintercept=-1.96,color="red"))+
          geom_vline(aes(xintercept=1.96,color="red"))+
          xlab("LDL Uptake Z-score")+guides(color="none",shape="none")+theme_cowplot())
}
dev.off()

plot(Harvard_MGH_variants_w_results$LDL_efflux_mu_z_adj,Harvard_MGH_variants_w_results$aPC_Conservation)

library(caret)
library(randomForest)
library(caTools)
split <- sample(1:nrow(Harvard_MGH_variants_w_results), size = ceiling(.5*nrow(Harvard_MGH_variants_w_results)),replace = FALSE)

data_train <- Harvard_MGH_variants_w_results[split,]
data_test <- Harvard_MGH_variants_w_results[setdiff(c(1:nrow(Harvard_MGH_variants_w_results)),split),]

rf_e <- randomForest(LDL_eff~aPC_Conservation+aPC_Epigenetics_Active+abs_chromBPnet_logfc+liver_TLand, data=data_train%>%filter(!is.na(aPC_Conservation),!is.na(LDL_efflux_p_FDR_lt_.05)), proximity=TRUE,importance=TRUE)
print(rf_e)
plot(rf_e)
varImpPlot(rf_e,type=1)

rf_e_class <- randomForest(as.factor(LDL_efflux_p_FDR_lt_.05)~aPC_Conservation+aPC_Epigenetics_Active+abs_chromBPnet_logfc+liver_TLand, data=data_train%>%filter(!is.na(aPC_Conservation),!is.na(LDL_efflux_p_FDR_lt_.05)), proximity=TRUE,importance=TRUE)
print(rf_e_class)
plot(rf_e_class)
varImpPlot(rf_e_class,type=1)
confusionMatrix(rf_e_class)

rf_u <- randomForest(abs_LDL_uptake_z~aPC_Conservation+aPC_Epigenetics_Active+abs_chromBPnet_logfc+liver_TLand, data=data_train%>%filter(!is.na(aPC_Conservation),!is.na(LDL_uptake_p_FDR_lt_.05)), proximity=TRUE,importance=TRUE)
print(rf_u)
plot(rf_u)
varImpPlot(rf_u,type=1)

rf_u_class <- randomForest(as.factor(LDL_uptake_p_FDR_lt_.05)~aPC_Conservation+aPC_Epigenetics_Active+abs_chromBPnet_logfc+liver_TLand, data=data_train%>%filter(!is.na(aPC_Conservation),!is.na(LDL_uptake_p_FDR_lt_.05)), proximity=TRUE,importance=TRUE)
print(rf_u_class)
plot(rf_u_class)
varImpPlot(rf_u_class,type=1)


#####LDLRAP1 LD analysis
############
LDLRAP1_LD_out<-Harvard_MGH_variants_w_results%>%filter(grepl("LDLRAP1",genes_within_500000))%>%dplyr::select(rsID)

write_tsv(LDLRAP1_LD_out,file=paste0(data_wd,"/variants/LDLRAP1_LD_out.tsv"))

for(id in c("rs35345536","rs34594460")){
  TOPLD_LDLRAP1_results_rs<-read_csv(paste0(data_wd,"/variants/TOPLD_results_LDLRAP1_locus_",id,".csv"))
  
  TOPLD_LDLRAP1_results_rs$rsID<-TOPLD_LDLRAP1_results_rs$rsID1
  TOPLD_LDLRAP1_results_rs$rsID[TOPLD_LDLRAP1_results_rs$rsID==id]<-TOPLD_LDLRAP1_results_rs$rsID2[TOPLD_LDLRAP1_results_rs$rsID==id]
  
  print(id)
  print(data.frame(Harvard_MGH_variants_w_results%>%filter(rsID%in%TOPLD_LDLRAP1_results_rs$rsID)%>%dplyr::select(rsID,Predicted_Functional,LDL_uptake_p_FDR,LDL_efflux_p_FDR)%>%arrange(LDL_uptake_p_FDR)))
}

Harvard_MGH_variants_w_results%>%filter(rsID%in%TOPLD_LDLRAP1_results_rs34594460$rsID)%>%dplyr::select(rsID,Predicted_Functional,LDL_uptake_p_FDR,LDL_efflux_p_FDR)


LDLRAP1_LD<-Harvard_MGH_variants_w_results%>%filter(grepl("LDLRAP1",genes_within_500000))



############ end LDLRAP1 
Harvard_MGH_variants_w_results<-Harvard_MGH_variants_w_results%>%dplyr::select(-ref.y,-alt.y)%>%rename(ref=ref.x,alt=alt.x)

SPDI_ref_file<-read_tsv("~/OneDrive - Harvard University/PostDoc/Data/variants/GRCh38p14_chrom_map.tsv")%>%dplyr::select(chrom,refseq)
SPDI_ref_file$chrom=as.numeric(SPDI_ref_file$chrom)
variants<-Harvard_MGH_variants_w_results%>%dplyr::select(variant,chr,Position,ref,alt)
variants<-left_join(variants,SPDI_ref_file,by=c("chr"="chrom"))
variants$SPDI_ID<-paste0(variants$refseq,":",variants$Position-1,":",variants$ref,":",variants$alt)
#Remove indels, we dont make predictions for those anyways
variants<-variants%>%filter(!is.na(ref))
SPDI_variants<-variants%>%dplyr::select(variant,SPDI_ID,refseq)
Harvard_MGH_variants_w_results_SPDI<-left_join(Harvard_MGH_variants_w_results,SPDI_variants,by="variant")%>%rename(chrRefSeqID=refseq,position=Position,ReferenceAllele=ref,AlternativeAllele=alt,SPDI=SPDI_ID)%>%dplyr::select(chrRefSeqID,chr,position,ReferenceAllele,AlternativeAllele,SPDI,rsID,everything())

write_tsv(Harvard_MGH_variants_w_results_SPDI,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_variants_w_results_SPDI.tsv")

write_tsv(Harvard_MGH_variants_w_results,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_variant_predictions_w_experimental_results.tsv")

lol<-Harvard_MGH_variants_w_results%>%filter(Predicted_Functional)%>%dplyr::select(chr,Position,ref.x,alt.x)

colnames(lol)<-c("chr","position","REF","ALT")

rofl<-Harvard_MGH_variants_w_results%>%filter(Predicted_Functional)%>%dplyr::select(chr,Position,ref.x,alt.x)
colnames(rofl)<-c("chr","position","REF","ALT")
write_csv(rofl,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Predicted_Functional_Harvard_MGH_variants.csv")

table(rofl$in_Hep_ATAC_peak,rofl$sum_pred_func)

subset2<-Harvard_MGH_variants%>%dplyr::select(variant,contains("GNOMAD"))%>%
  pivot_longer(cols=-c("variant"),names_to="class",values_to="GNOMAD_AF")

subset2$class<-as.factor(subset2$class)
subset2$temp_class<-gsub("_GNOMAD_AF","",subset2$class)
subset2$sex<-NA
subset2$sex[grepl("Male",subset2$temp_class)]<-"Male"
subset2$sex[grepl("Female",subset2$temp_class)]<-"Female"
subset2<-subset2%>%filter(!is.na(subset2$sex),!temp_class%in%c("Male","Female"))

subset2$population<-gsub("_.*","",subset2$temp_class)

cols<-viridis_pal()(length(unique(subset2$class)))




subset<-Harvard_MGH_variants%>%dplyr::select(variant,AFR_GNOMAD_AF,AMR_GNOMAD_AF,
                              EAS_GNOMAD_AF,NFE_GNOMAD_AF,SAS_GNOMAD_AF,FIN_GNOMAD_AF,
                              AMI_GNOMAD_AF,ASJ_GNOMAD_AF,OTH_GNOMAD_AF,Total_GNOMAD_AF)%>%
  pivot_longer(cols=2:10,names_to="class",values_to="GNOMAD_AF")
subset$population<-gsub("_GNOMAD_AF","",subset$class)
subset$population<-as.factor(subset$population)

cols<-viridis_pal()(length(unique(subset$population)))

ggplot(subset,aes(x=population,y=GNOMAD_AF,color=class,fill=class))+geom_boxplot()+scale_color_manual(name="Category",values=cols)+ggplot2::theme(legend.position = "none")

pdf("~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_boxplots.pdf",width=8,height=6)
ggplot(subset,aes(x=population,y=GNOMAD_AF,fill=class))+geom_boxplot()+ggplot2::theme(legend.position = "none")
ggplot(subset2,aes(x=population,y=GNOMAD_AF,fill=sex))+geom_boxplot()
dev.off()

#######################################
# 
#######################################



break
Harvard_MGH_var_out<-Harvard_MGH_variants%>%dplyr::select(variant,rsID)%>%separate(variant,c("chr","start","ref","alt"))
Harvard_MGH_var_out[,1]<-as.numeric(Harvard_MGH_var_out%>%pull(1))
Harvard_MGH_var_out[,2]<-as.numeric(Harvard_MGH_var_out%>%pull(2))
Harvard_MGH_var_out$end<-Harvard_MGH_var_out$start
Harvard_MGH_var_out<-Harvard_MGH_var_out%>%dplyr::select(chr,start,end,ref,alt,rsID)
Harvard_MGH_var_out$chr<-paste0("chr",Harvard_MGH_var_out%>%pull(chr))

write_tsv(Harvard_MGH_var_out,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_variants_for_Christopher_hg38.tsv")

######################################################################
# liftOver
######################################################################
setwd(data_wd)
ch=import.chain("Hg38Tohg19.over.chain")

temp_df<-Harvard_MGH_variants%>%mutate(hg38_variant=variant)%>%dplyr::select(hg38_variant,variant)%>%separate(variant,c("chr","start","ref","alt"))
temp_df$start<-as.numeric(temp_df$start)
temp_df$end<-temp_df$start
temp_df$chr<-paste0("chr",temp_df$chr)
hg38_GRanges<-makeGRangesFromDataFrame(temp_df,start.field="start",end.field="end",keep.extra.columns = TRUE)

Harvard_MGH_hg19<-annoGR2DF(unlist(liftOver(hg38_GRanges,ch)))
Harvard_MGH_hg19_out<-Harvard_MGH_hg19%>%dplyr::select(chr,start,ref,alt)
Harvard_MGH_hg19_out$chr<-gsub("chr","",Harvard_MGH_hg19_out$chr)
colnames(Harvard_MGH_hg19_out)[colnames(Harvard_MGH_hg19_out)=="start"]<-"pos"

write_tsv(Harvard_MGH_hg19_out,file="~/OneDrive - Harvard University/PostDoc/Data/variants/Harvard_MGH_hg19_out.tsv")

######################################################################
# end
######################################################################
