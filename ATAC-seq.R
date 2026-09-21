library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
setwd("/home/data/sdzl24/Project/aging/output_for_ATAC/")

dat1 <- read.csv("sample_ID_aligment.txt",header =T,sep = "\t")%>% data.frame
dat1 <- (dat1[,c(2,1)])%>%unique
colnames(dat1) <- c("SampleID","Tissue")
dat1$Factor <- dat1$Tissue
dat1$Replicate <- c(1:6,1:6,1:5)
dat1$bamReads <-NA 
dat1$Peaks <-NA
for(i in 1:nrow(dat1)){
    t_ID <- dat1[i,"SampleID"]
    dat1[i,"bamReads"] <-paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/",t_ID,"/alignment/rmdup/",t_ID,".bam")
    dat1[i,"Peaks"] <-paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/",t_ID,"/peak/filter_q0.05/",t_ID,"_peaks.narrowPeak")
    print(i)
}
dat1$PeakCaller <- "MACS2"
dat1$PeakFormat <- "narrow"

write.csv(dat1,"06_for_Diffbind.csv",row.names=F)
dat_AB <- filter(dat1,Tissue%in%c("A","B"))
dat_BC <- filter(dat1,Tissue%in%c("B","C"))
write.csv(dat_AB,"06_for_Diffbind_A_B.csv",row.names=F)
write.csv(dat_BC,"06_for_Diffbind_B_C.csv",row.names=F)

library(DiffBind)
lapply(c("A_B","B_C"),function(i){
    print(c(i,"start"))
    dbObj <- dba(sampleSheet=paste0("06_for_Diffbind_",i,".csv"))
    dbObj1 <- dba.blacklist(dbObj, blacklist=DBA_BLACKLIST_HG38,greylist=FALSE)
    peakdata.BL <- dba.show(dbObj1)$Intervals
    #########################################
    dbObj2 <- dba.count(dbObj1, bUseSummarizeOverlaps=TRUE)
    pdf(paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/figure/06_pca_heatmap_after_count.pdf"))
    dba.plotPCA(dbObj2, attributes=DBA_TISSUE, label=DBA_ID)
    plot(dbObj2)
    dev.off()
    dbObj <- dba.contrast(dbObj2, categories=DBA_TISSUE,minMembers = 2)
    dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
    dba.show(dbObj, bContrasts=T)
    pdf(paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/figure/06_VENN.pdf"))
    dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
    dev.off()
    comp1.edgeR <- DiffBind::dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)%>%data.frame
    comp1.deseq <- DiffBind::dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)%>%data.frame
    #to bed
    comp1.edgeR$start <- comp1.edgeR$start -1
    comp1.deseq$start <- comp1.deseq$start -1 
    write.table(comp1.edgeR,paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_edgeR.txt"), sep="\t", quote=F, col.names = NA)
    edge.bed <- filter(comp1.edgeR,FDR < 0.05)
    write.table(edge.bed, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_edgeR_sig.bed"), sep="t", quote=F, row.names=F, col.names=F)
    write.table(comp1.deseq,paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2.txt"), sep="\t", quote=F, col.names = NA)
    deseq.bed <-filter(comp1.deseq,FDR < 0.05)
    write.table(deseq.bed, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
    deseq.bed$Method <- "DESeq2"
    edge.bed$Method <- "edgeR"
    Fdat <- bind_rows(deseq.bed,edge.bed)
    Fdat$type=ifelse( Fdat$Fold >0,'hyper', "hypo")
    save(Fdat,file=paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig.Rdata"))
    write.table(Fdat, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
    system(paste0("less /home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig.bed |cut -f1-3,12-13|sort -k1,1 -k2,2n |gzip >/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig_sort.bed.gz"))
    print(c(i,"finish"))
})


library(ChIPseeker)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(clusterProfiler)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
head(seqlevels(txdb))
library(ggsci)

library(Hmisc)

lapply(c("A_B","A_C"),function(i){
    print(c(i,"start"))
    peak=readPeakFile(peakfile=paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig.bed"), header=F, as = 'GRanges')
    peakAnno=annotatePeak(peak = peak, tssRegion=c(-3000,3000), TxDb = txdb, annoDb="org.Hs.eg.db")
    peakAnnodf <- data.frame(peakAnno)
    peakAnnodf$start <- peakAnnodf$start -1
    peakAnnodf <- peakAnnodf[,-c(4:5)]
    if(i=="A_B"){
        colnames(peakAnnodf)[4:13] <-c("width","strand","Conc","Conc_B","Conc_A","Fold","p.value","FDR","Method","type")
    }else if(i=="A_C"){
        colnames(peakAnnodf)[4:13] <-c("width","strand","Conc","Conc_C","Conc_A","Fold","p.value","FDR","Method","type")
    }

    peakAnnodf$anno2 <- peakAnnodf$annotation
    peakAnnodf$anno2 <- gsub("Promoter .*", "Promoter", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Intron .*", "Intron", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Exon .*", "Exon", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("3' UTR", "3'UTR", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("5' UTR", "5'UTR", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Distal Intergenic .*", "Intergenic", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Distal Intergenic", "Intergenic", peakAnnodf$anno2)
    peakAnnodf$anno2 <- gsub("Downstream .*", "Intergenic", peakAnnodf$anno2)
    write.table(peakAnnodf, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig_anno.bed.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    system(paste0("sed -n '2,$p' /home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig_anno.bed.txt |cut -f1-3,12-13,20-21,23-24,26 |sort -k1,1 -k2,2n| gzip> /home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_diff_deseq2_edgeR_sig_anno_sort.bed.gz"))
   
    hyper <- filter(peakAnnodf,type=="hyper")
    hypo <- filter(peakAnnodf,type=="hypo")
    write.table(hyper, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_hyper_deseq2_edgeR_sig_anno.bed.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    write.table(hypo, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_hypo_deseq2_edgeR_sig_anno.bed.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    lapply(c("hyper","hypo"),function(j){
        print(c(j,"start"))
        df_tmp <- read.csv(paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_",j,"_deseq2_edgeR_sig_anno.bed.txt"),header = T,sep = "\t") %>% data.frame
        df <-table(df_tmp$anno2)%>%data.frame
        df$per <- df$Freq/sum(df$Freq)*100
        df$per<- round(df$per,2)
        df$per1 <- paste0(as.character(df$per),"%")
        # mycolor = pal_d3("category20")(6)
        mycolor = pal_jama("default")(6)
        write.table(df, paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_",j,"_deseq2_edgeR_anno_statistic.txt"), sep="\t", quote=F, row.names=F, col.names=T)
        pdf(paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/figure/06_pie_ATAC_",j,"_deseq2_edgeR_sig_anno.pdf"),width=6,height=5)
            # pie(df$Freq, cex=1.5,col = mycolor,labels = df$per1, radius = 1,main=paste0(i,": ",j))
            pie(df$Freq, cex=1.5,col = mycolor,labels = df$per1, radius = 1,main=paste0(j))
            legend("topright", c("3'UTR","5'UTR","Exon","Intergenic","Intron","Promoter"), cex = 1, fill = mycolor)
        dev.off()
    })
    print(c(i,"finish"))

    df <- unique(peakAnnodf[,c(1:3,13,26)])%>%group_by(type,anno2)%>%summarise(number=n())%>%data.frame
    df_totals <- df %>%group_by(type) %>%summarise(total = sum(number))%>%data.frame
    df <-left_join(df,df_totals,by="type")
    df$percentage <- df$number/df$total*100
    df$type <- gsub("hyper","Hyper",df$type)
    df$type <- gsub("hypo","Hypo",df$type)
    save(df,file=paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_all_deseq2_edgeR_anno_statistic.Rdata"))
    write.table(df,paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",i,"/06_ATAC_all_deseq2_edgeR_anno_statistic.txt"),row.names = F, col.names = T,quote =F,sep="\t")
    
})


library(clusterProfiler)
library(org.Hs.eg.db)

lapply(c("hyper","hypo"),function(i){
    lapply(c("A_B","A_C"),function(j){
        # sub_dir <- paste0("/home/data/sdzl24/Project/aging/output/Diff/",j,"/",i)
        sub_dir <- paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",j,"/",i)
            if(dir.exists(sub_dir)){
                print(c(sub_dir, "is exists\n"))
            }else{
                dir.create(sub_dir,recursive=TRUE)
            }
        setwd(sub_dir)
        fdat <- read.csv(paste0("/home/data/sdzl24/Project/aging/output_for_ATAC/ATAC/merge/all/Diff/",j,"/06_ATAC_",i,"_deseq2_edgeR_sig_anno.bed.txt"),header = T,sep = "\t") %>% data.frame
        fdat <- filter(fdat,!(is.na(geneId)))
        promoter <- filter(fdat,anno2=="Promoter")
        write.table(promoter,"08_diff_ATAC_promoter.bed.txt",row.names = F, col.names = T,quote =F,sep="\t")
        k="promoter"
        genes <- unique(promoter[,c("geneId","SYMBOL")])
        colnames(genes) <- c("ENTREZID","SYMBOL")
        save(genes,file=paste0("08_",k,"_regulated_gene.Rdata"))
        write.table(genes,paste0("08_",k,"_regulated_gene.txt"),row.names = F, col.names = T,quote =F,sep="\t")
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
                return(go1)
            }else{return(NULL)}
        })
        #============================================
        ego1 <- do.call(rbind,tmp)
        if(is.null(ego1)){
            print(paste0("go no reslut"))
        }else{
            write.table(ego1,paste0(k,"_go.txt"),col.names=T,row.names =F,quote=F,sep="\t")
        }
        print(c(i,j,k))
    })
})
