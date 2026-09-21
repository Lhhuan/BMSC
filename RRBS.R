library(methylKit)
library(dplyr)

setwd("/home/data/sdzl24/Project/aging/output/")

file_list1=list("./J31/03_alignment/J31_sorted.bam","./J32/03_alignment/J32_sorted.bam","./J33/03_alignment/J33_sorted.bam","./J34/03_alignment/J34_sorted.bam","./J35/03_alignment/J35_sorted.bam","./J36/03_alignment/J36_sorted.bam","./J37/03_alignment/J37_sorted.bam","./J38/03_alignment/J38_sorted.bam","./J39/03_alignment/J39_sorted.bam","./J40/03_alignment/J40_sorted.bam","./J41/03_alignment/J41_sorted.bam","./J42/03_alignment/J42_sorted.bam")
obj1 = processBismarkAln(location=file_list1, sample.id=list("J31","J32","J33","J34","J35","J36","J37","J38","J39","J40","J41","J42"),assembly="hg38", save.folder=NULL,
  save.context = c("CpG"), read.context = "CpG", nolap = FALSE,
  mincov= 10, minqual = 20, phred64 = FALSE, treatment=c(rep(0,6),rep(1,6)),
  save.db = FALSE) 
save(obj1,file="./Diff/A_B/05_A_B_obj.Rdata")

sample=c("J31","J32","J33","J34","J35","J36","J37","J38","J39","J40","J41","J42")
lapply(1:12,function(i){
    pdf(paste0("./",sample[i],"/05_CoverageStats.pdf"),width=8,height=8)
        getCoverageStats(obj1[[i]],plot=TRUE,both.strands=FALSE)
    dev.off()
    print(i)
})

lapply(1:12,function(i){
    pdf(paste0("./",sample[i],"/05_MethylationStats.pdf"),width=8,height=8)
        getMethylationStats(obj1[[i]],plot=TRUE,both.strands=FALSE)
    dev.off()
    print(i)
})

filtered_obj<- filterByCoverage(obj1,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
norm_obj <- normalizeCoverage(filtered_obj)

#####################################################################################
#                                   DMR
#####################################################################################
regions <- tileMethylCounts(norm_obj,win.size=1000,step.size=1000,cov.bases = 0)
meth =unite(regions)
#######
meth_levels=percMethylation(meth,rowids=TRUE)%>%data.frame
meth_levels1 <- bind_cols(getData(meth)[,1:3],meth_levels)
write.table(meth_levels1,"./Diff/A_B/05_methylation_level.txt", sep="\t", quote=F, row.names=F, col.names=T)
meth_levels1$start <- meth_levels1$start - 1
write.table(meth_levels1,"./Diff/A_B/05_methylation_level.bed.txt", sep="\t", quote=F, row.names=F, col.names=T)
##################################
pdf("./Diff/A_B/05_MethylationCorrelation.pdf",width=8,height=8)
getCorrelation(meth,plot=TRUE)
dev.off()

pdf("./Diff/A_B/05_cluster_sample.pdf",width=8,height=8)
clusterSamples(meth,plot=TRUE)
dev.off()

pdf("./Diff/A_B/05_PCA_screeplot.pdf",width=8,height=8)
PCASamples(meth,screeplot=TRUE)
dev.off()

pdf("./Diff/A_B/05_PCASamples.pdf",width=8,height=8)
PCASamples(meth)
dev.off()
####################################
myDiff = calculateDiffMeth(meth,mc.cores=2)
#######
library(ggplot2)
library(cowplot)
library(tibble)
library(ggrepel)

p_theme<-theme(panel.grid =element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

exclude_chr <- getData(myDiff) %>% filter(!(chr %in%paste0("chr",1:22)))

pdf("./Diff/A_B/05_diffMethPerChr.pdf",height = 6,width = 6)
diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude = unique(exclude_chr$chr))
dev.off()

myDiff25p_hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p_hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff_tiles25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

myDiffall=getMethylDiff(myDiff,difference=1,qvalue=1)
df=myDiffall
df$type=ifelse(df$qvalue>0.05,'none',
            ifelse( df$meth.diff >25,'hyperDMR', 
                    ifelse( df$meth.diff< -25,'hypoDMR','none') )
)
head(df)
df1 <-getData(df)
df1$class <- "A_B"
save(df1,file="./Diff/A_B/05_diff_region.Rdata")
position_math <- unite(norm_obj)

pdf("./Diff/A_B/05_position_MethylationCorrelation.pdf",width=8,height=8)
getCorrelation(position_math,plot=TRUE)
dev.off()

pdf("./Diff/A_B/05_position_cluster_sample.pdf",width=8,height=8)
clusterSamples(position_math,plot=TRUE)
dev.off()

pdf("./Diff/A_B/05_position_PCA_screeplot.pdf",width=8,height=8)
PCASamples(position_math,screeplot=TRUE)
dev.off()

pdf("./Diff/A_B/05_position_PCASamples.pdf",width=8,height=8)
PCASamples(position_math)
dev.off()

myDiff_DMP = calculateDiffMeth(position_math,mc.cores=2)
exclude_chr <- getData(myDiff_DMP) %>% filter(!(chr %in%paste0("chr",1:22)))
pdf("./Diff/A_B/05_diffMeth_DMP_PerChr.pdf",height = 6,width = 6)
diffMethPerChr(myDiff_DMP,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude = unique(exclude_chr$chr))
dev.off()


myDiffall_DMP=getMethylDiff(myDiff_DMP,difference=1,qvalue=1)
df_DMP=myDiffall_DMP
df_DMP$type=ifelse(df_DMP$qvalue>0.05,'none',
            ifelse( df_DMP$meth.diff >25,'hyperDMP', 
                    ifelse( df_DMP$meth.diff< -25,'hypoDMP','none') )
)
head(df_DMP)
df_DMP1 <-getData(df_DMP)
df_DMP1$class <- "A_B"
save(df_DMP1,file="./Diff/A_B/05_diff_DMP.Rdata")
################################################################# B_C ####################################################
file_list2=list("./J37/03_alignment/J37_sorted.bam","./J38/03_alignment/J38_sorted.bam","./J39/03_alignment/J39_sorted.bam","./J40/03_alignment/J40_sorted.bam","./J41/03_alignment/J41_sorted.bam","./J42/03_alignment/J42_sorted.bam","./J43/03_alignment/J43_sorted.bam","./J44/03_alignment/J44_sorted.bam","./J45/03_alignment/J45_sorted.bam","./J47/03_alignment/J47_sorted.bam","./J48/03_alignment/J48_sorted.bam","./J49/03_alignment/J49_sorted.bam")

obj2 = processBismarkAln(location=file_list2, sample.id=list("J37","J38","J39","J40","J41","J42","J43","J44","J45","J47","J48","J49"),assembly="hg38", save.folder=NULL,
  save.context = c("CpG"), read.context = "CpG", nolap = FALSE,
  mincov= 10, minqual = 20, phred64 = FALSE, treatment=c(rep(0,6),rep(1,6)),
  save.db = FALSE)

save(obj2,file="./Diff/B_C/05_B_C_obj.Rdata")

sample=c("J37","J38","J39","J40","J41","J42","J43","J44","J45","J47","J48","J49")
lapply(1:12,function(i){
    pdf(paste0("./",sample[i],"/05_CoverageStats.pdf"),width=8,height=8)
        getCoverageStats(obj2[[i]],plot=TRUE,both.strands=FALSE)
    dev.off()
    print(i)
})

lapply(1:12,function(i){
    pdf(paste0("./",sample[i],"/05_MethylationStats.pdf"),width=8,height=8)
        getMethylationStats(obj2[[i]],plot=TRUE,both.strands=FALSE)
    dev.off()
    print(i)
})

# filter 
filtered_obj<- filterByCoverage(obj2,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
norm_obj <- normalizeCoverage(filtered_obj)

#####################################################################################
#                                   DMR
#####################################################################################

regions <- tileMethylCounts(norm_obj,win.size=1000,step.size=1000,cov.bases = 0)
meth.min =unite(regions,min.per.group=1L)
meth =unite(regions)
meth_levels=percMethylation(meth,rowids=TRUE)%>%data.frame
meth_levels1 <- bind_cols(getData(meth)[,1:3],meth_levels)
write.table(meth_levels1,"./Diff/B_C/05_methylation_level.txt", sep="\t", quote=F, row.names=F, col.names=T)
meth_levels1$start <- meth_levels1$start - 1
write.table(meth_levels1,"./Diff/B_C/05_methylation_level.bed.txt", sep="\t", quote=F, row.names=F, col.names=T)
pdf("./Diff/B_C/05_MethylationCorrelation.pdf",width=8,height=8)
getCorrelation(meth,plot=TRUE)
dev.off()


pdf("./Diff/B_C/05_cluster_sample.pdf",width=8,height=8)
clusterSamples(meth,plot=TRUE)
dev.off()

pdf("./Diff/B_C/05_PCA_screeplot.pdf",width=8,height=8)
PCASamples(meth,screeplot=TRUE)
dev.off()

pdf("./Diff/B_C/05_PCASamples.pdf",width=8,height=8)
PCASamples(meth)
dev.off()

myDiff = calculateDiffMeth(meth,mc.cores=2)

pdf("./Diff/B_C/05_diffMethPerChr.pdf",height = 6,width = 6)
diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude = unique(exclude_chr$chr))
dev.off()

myDiff25p_hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
myDiff25p_hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
myDiff_tiles25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

myDiffall=getMethylDiff(myDiff,difference=1,qvalue=1)
df=myDiffall
df$type=ifelse(df$qvalue>0.05,'none',
            ifelse( df$meth.diff >25,'hyperDMR', 
                    ifelse( df$meth.diff< -25,'hypoDMR','none') )
)
head(df)
df1 <-getData(df)
df1$class <- "B_C"
save(df1,file="./Diff/B_C/05_diff_region.Rdata")
position_math <- unite(norm_obj)

pdf("./Diff/B_C/05_position_MethylationCorrelation.pdf",width=8,height=8)
getCorrelation(position_math,plot=TRUE)
dev.off()

pdf("./Diff/B_C/05_position_cluster_sample.pdf",width=8,height=8)
clusterSamples(position_math,plot=TRUE)
dev.off()

pdf("./Diff/B_C/05_position_PCA_screeplot.pdf",width=8,height=8)
PCASamples(position_math,screeplot=TRUE)
dev.off()

pdf("./Diff/B_C/05_position_PCASamples.pdf",width=8,height=8)
PCASamples(position_math)
dev.off()

myDiff_DMP = calculateDiffMeth(position_math,mc.cores=2)
exclude_chr <- getData(myDiff_DMP) %>% filter(!(chr %in%paste0("chr",1:22)))
pdf("./Diff/B_C/05_diffMeth_DMP_PerChr.pdf",height = 6,width = 6)
diffMethPerChr(myDiff_DMP,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude = unique(exclude_chr$chr))
dev.off()

myDiffall_DMP=getMethylDiff(myDiff_DMP,difference=1,qvalue=1)
df_DMP=myDiffall_DMP
df_DMP$type=ifelse(df_DMP$qvalue>0.05,'none',
            ifelse( df_DMP$meth.diff >25,'hyperDMP', 
                    ifelse( df_DMP$meth.diff< -25,'hypoDMP','none') )
)
head(df_DMP)
df_DMP1 <-getData(df_DMP)
df_DMP1$class <- "B_C"
save(df_DMP1,file="./Diff/B_C/05_diff_DMP.Rdata")

#################################################################################

setwd("/home/data/sdzl24/Project/aging/output/Diff/")
load("./A_B/05_diff_region.Rdata")
A_B <- df1
load("./B_C/05_diff_region.Rdata")
B_C <- df1
rm(df1)
df <- bind_rows(A_B,B_C,A_C)
df <- filter(df,type!="none")
df <-filter(df,chr%in%paste0("chr",1:22))
df$start <- df$start -1 
df$start <-as.integer(df$start)
save(df,file="06_3_class_diff_region.Rdata")
write.table(df,"06_3_class_diff_region.bed.txt",row.names = F, col.names = T,quote =F,sep="\t")


library(ChIPseeker)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(clusterProfiler)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


########### gene position
lapply(unique(df$class),function(i){
    tmp_df <- filter(df,class==i)
    write.table(tmp_df,paste0("./",i,"/06_diff_region.bed.txt"),row.names = F, col.names = T,quote =F,sep="\t")
    write.table(tmp_df[,1:3],paste0("./",i,"/06_diff_region.bed"),row.names = F, col.names = F,quote =F,sep="\t")
    peak=readPeakFile(peakfile=paste0("./",i,"/06_diff_region.bed.txt"), header=T, as = 'GRanges')
    peakAnno=annotatePeak(peak = peak, tssRegion=c(-3000,3000), TxDb = txdb, annoDb="org.Hs.eg.db")
    peakAnnodf <- data.frame(peakAnno)
    peakAnnodf$start <- peakAnnodf$start -1
    peakAnnodf <- peakAnnodf[,-c(6)]
    peakAnnodf$anno2 <- peakAnnodf$annotation
    peakAnnodf$anno2 <- gsub("Promoter .*", "Promoter", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Intron .*", "Intron", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Exon .*", "Exon", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("3' UTR", "3'UTR", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("5' UTR", "5'UTR", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Distal Intergenic .*", "Intergenic", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Distal Intergenic", "Intergenic", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Downstream .*", "Intergenic", peakAnnodf$anno2)
    save(peakAnnodf,file=paste0("./",i,"/06_diff_region_chipseeker_anno.bed.Rdata"))
    write.table(peakAnnodf,paste0("./",i,"/06_diff_region_chipseeker_anno.bed.txt"),row.names = F, col.names = T,quote =F,sep="\t")
     df <- unique(peakAnnodf[,c(1:3,9,23)])%>%group_by(type,anno2)%>%summarise(number=n())%>%data.frame
    df_totals <- df %>%group_by(type) %>%summarise(total = sum(number))%>%data.frame
    df <-left_join(df,df_totals,by="type")
    df$percentage <- df$number/df$total*100
    df$type <- gsub("hyperDMR","Hyper-DMR",df$type)
    df$type <- gsub("hypoDMR","Hypo-DMR",df$type)
    save(df,file=paste0("./",i,"/06_DMR_chipseeker_anno_statistic.Rdata"))
    write.table(df,paste0("./",i,"/06_DMR_chipseeker_anno_statistic.txt"),row.names = F, col.names = T,quote =F,sep="\t")
})

############   CGI
library(annotatr)
lapply(unique(df$class),function(i){
    print(i)
    tmp_df <- filter(df,class==i)
    dm_file=paste0("./",i,"/06_diff_region.bed")
    dm_regions = read_regions(con = dm_file, genome = 'hg38')

    annots = c('hg38_cpgs')
    annotations = build_annotations(genome = 'hg38', annotations = annots)  #1to5kb upstream of TSS
    dm_annotated = annotate_regions(
        regions = dm_regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)
    df_dm_annotated = data.frame(dm_annotated)
    df_dm_annotated$start <- df_dm_annotated$start -1
    df_dm_annotated$annot.start <- df_dm_annotated$annot.start -1
    GCI_anno <- unique(df_dm_annotated[,c(6:8,11,15)])
    write.table(GCI_anno,paste0("./",i,"/06_for_DMR_CPG_anno.bed"),row.names = F, col.names = F,quote =F,sep="\t")
    system(paste0("less ./",i,"/06_for_DMR_CPG_anno.bed |sort -k1,1 -k2,2n >./",i,"/06_for_DMR_CPG_anno_sorted.bed"))
    system(paste0("bedtools intersect -a  ./",i,"/06_diff_region.bed -b ./",i,"/06_for_DMR_CPG_anno_sorted.bed -wo > ./",i,"/06_DMR_CPG_anno_from_bedtools.bed" ))
    df_dm_annotated <- read.csv(paste0("./",i,"/06_DMR_CPG_anno_from_bedtools.bed"),header =F,sep = "\t")%>% data.frame
    colnames(df_dm_annotated)<-c("chr","start","end","annot.chr","annot.start","annot.end","annot.id","annot.type","overlap_width")
    tmp_anno <- df_dm_annotated%>%group_by(chr,start,end)%>%summarise(overlap_width=max(overlap_width))%>% data.frame
    df_dm_annotated1 <- left_join(tmp_anno,df_dm_annotated,by=c("chr","start","end","overlap_width"))
    all_duplicated_rows <- duplicated(df_dm_annotated1[,1:4]) | duplicated(df_dm_annotated1[,1:4], fromLast = TRUE)
    df_dm_annotated1[all_duplicated_rows,]
    df_dm_annotated1 <- df_dm_annotated1[-46,]

    colnames(df_dm_annotated)[1] <-"chr"
    fdat <- left_join(tmp_df,unique(df_dm_annotated1[,c(1:3,8:9)]),by=c("chr","start","end"))
    fdat$annot.type<- gsub("hg38_cpg_shores","CPG shore",fdat$annot.type)
    fdat$annot.type <- gsub("hg38_cpg_islands","CPG island",fdat$annot.type)
    fdat$annot.type <- gsub("hg38_cpg_inter","Open sea",fdat$annot.type)
    fdat$annot.type <- gsub("hg38_cpg_shelves","CPG shelf",fdat$annot.type)
    save(fdat,file=paste0("./",i,"/06_DMR_CPG_anno.bed.Rdata"))
    write.table(fdat,paste0("./",i,"/06_DMR_CPG_anno.bed.txt"),row.names = F, col.names = T,quote =F,sep="\t")

    df <- unique(fdat[,c(1:3,8,11)])%>%group_by(type,annot.type)%>%summarise(number=n())%>%data.frame
    df_totals <- df %>%group_by(type) %>%summarise(total = sum(number))%>%data.frame
    df <-left_join(df,df_totals,by="type")
    df$percentage <- df$number/df$total*100
    df$type <- gsub("hyperDMR","Hyper-DMR",df$type)
    df$type <- gsub("hypoDMR","Hypo-DMR",df$type)
    save(df,file=paste0("./",i,"/06_DMR_CGI_anno_statistic.Rdata"))
    write.table(df,paste0("./",i,"/06_DMR_CGI_anno_statistic.txt"),row.names = F, col.names = T,quote =F,sep="\t")
    df_totals <- unique(df[,c(1,4)])
    df<- filter(df,annot.type!="NA")
})



##################enrichment

library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

lapply(c("hyperDMR","hypoDMR"),function(i){
    lapply(c("A_B","B_C"),function(j){
        sub_dir <- paste0("/home/data/sdzl24/Project/aging/output/Diff/",j,"/",i)
            if(dir.exists(sub_dir)){
                print(c(sub_dir, "is exists\n"))
            }else{
                dir.create(sub_dir,recursive=TRUE)
            }
        setwd(sub_dir)
        load("../06_diff_region_chipseeker_anno.bed.Rdata") #peakAnnodf
        fdat <- filter(peakAnnodf,!(is.na(geneId)))
        promoter <- filter(fdat,anno2=="Promoter"&type==i) 
        write.table(promoter,"06_diff_region_promoter.bed.txt",row.names = F, col.names = T,quote =F,sep="\t")
        k <- "promoter"
        dd <- read.csv(paste0("06_diff_region_",k,".bed.txt"),header =T,sep = "\t")%>% data.frame
        genes <- unique(dd[,c("geneId","SYMBOL")])
        colnames(genes) <- c("ENTREZID","SYMBOL")
        
        save(genes,file=paste0("06_",k,"_regulated_gene.Rdata"))
        write.table(genes,paste0("06_",k,"_regulated_gene.txt"),row.names = F, col.names = T,quote =F,sep="\t")
        kegg <- enrichKEGG(genes$ENTREZID, 
            organism = 'hsa', 
            keyType = 'kegg', 
            pvalueCutoff = 0.05,
            pAdjustMethod = 'fdr', 
            minGSSize = 10,
            maxGSSize = 500,
            qvalueCutoff = 0.05,
            use_internal_data = FALSE)
            kegg1 <-data.frame(kegg)
        if(nrow(kegg1)>0){
            kegg1$Category <-"KEGG"
            kegg1$Generatio1 <-sapply(1:nrow(kegg1),function(x){
                as.numeric(strsplit(kegg1$GeneRatio[x],split = '/')[[1]][1])/as.numeric(strsplit(kegg1$GeneRatio[x],split = '/')[[1]][2])
            })
            kegg1 <-kegg1[order(kegg1$Generatio1),]
            kegg1$Description<-factor(kegg1$Description,levels=kegg1$Description)
            write.table(kegg1,paste0(k,"_KEGG.txt"),col.names=T,row.names =F,quote=F,sep="\t")
            write.xlsx(kegg1s,paste0(k,"_KEGG.xlsx"), rowNames = FALSE)
        }

        types <-c("BP","CC","MF")
        #================================================
        tmp <-lapply(types,function(type){
            go <- enrichGO(genes$SYMBOL, 
                OrgDb = 'org.Hs.eg.db',
                keyType = 'SYMBOL',
                ont=type,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'fdr', 
                minGSSize = 10,
                maxGSSize = 500,
                qvalueCutoff = 0.05,   
                readable = FALSE,
                pool = FALSE)
            go1 <-data.frame(go)
            if(nrow(go1)>1){
                go1$Category <-type
                go1$Generatio1 <-sapply(1:nrow(go1),function(x){
                    as.numeric(strsplit(go1$GeneRatio[x],split = '/')[[1]][1])/as.numeric(strsplit(go1$GeneRatio[x],split = '/')[[1]][2])
                })
                go1 <-go1[order(go1$Generatio1),]
                go1$Description<-factor(go1$Description,levels=go1$Description)
                write.table(go1,paste0(k,"_GO_",type,".txt"),col.names=T,row.names =F,quote=F,sep="\t")
                write.xlsx(go1,paste0(k,"_GO_",type,".xlsx"), rowNames = FALSE)
            #############
                return(go1)
            }else{return(NULL)}
        })
        #============================================
        ego1 <- do.call(rbind,tmp)
        if(is.null(ego1)){
            print(paste0("go no reslut"))
        }else{
            write.table(ego1,paste0(k,"_go.txt"),col.names=T,row.names =F,quote=F,sep="\t")
            write.xlsx(ego1,paste0(k,"_go.xlsx"), rowNames = FALSE)
        }
    })
})
