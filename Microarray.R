library(dplyr)
library(limma)

lapply(c("A_B","B_C"),function(i){
  setwd("/home/data/sdzl24/Project/aging/data/Microarray/Rawdata/")
  print(c(i,"start"))
  targets <- readTargets(paste0("/home/data/sdzl24/Project/aging/script_for_microarray/",i,"_targets.txt"))
  rawData <- read.maimages(targets, source="agilent", green.only=TRUE) 
  normData <- backgroundCorrect(rawData, method="normexp")
  normData <- normalizeBetweenArrays(normData, method="none") 

  design <- model.matrix(~ 0 + factor(c(rep(0,3),rep(1,3)))) 
  colnames(design) <- c("Control", "Treatment")
  fit <- lmFit(normData, design) 
  contrast <- makeContrasts(Treatment - Control, levels=design) 
  fit2 <- contrasts.fit(fit, contrast) 
  fit2 <- eBayes(fit2) 
  results <- topTable(fit2, adjust="BH",number=nrow(fit2))
  colnames(results)[11] <-"log2FC"
  #####################################
  rs1 <- results %>%filter(!startsWith(SystematicName, "ENST"))
  rs2 <-  results[grepl("^ENST", results$SystematicName), ]

  library(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  rs_ENST <- getBM(
    attributes = c("refseq_mrna", "ensembl_transcript_id"),
    filters = "refseq_mrna",
    values = rs1$SystematicName,
    mart = ensembl
  )%>%unique
  colnames(rs_ENST) <- c("SystematicName","ensembl_transcript_id")
  rs1<-left_join(rs1,rs_ENST,by="SystematicName")
  rs2$ensembl_transcript_id <- rs2$SystematicName
  Fresults <- bind_rows(rs1,rs2)
  rownames(Fresults)<-1:nrow(Fresults)
  ###############################################################
  gene1 <- Fresults[grepl("^ENST", Fresults$GeneName), ]
  gene2 <- filter(Fresults,!(GeneName%in%gene1$GeneName))
  gene_ensg_symbol <- getBM(
    attributes = c("hgnc_symbol", "ensembl_transcript_id"),
    # attributes = c("hgnc_symbol", "ensembl_transcript_id",'ensembl_gene_id'),
    filters =  "ensembl_transcript_id",
    values = gene1$GeneName,
    mart = ensembl
  )%>%unique
  colnames(gene_ensg_symbol)<-c("genename2","GeneName")
  gene1 <- left_join(gene1,gene_ensg_symbol,by="GeneName")
  gene2$genename2 <- gene2$GeneName
  Frs <- bind_rows(gene1,gene2)
  ###############################

  write.table(Frs, paste0(i,"/01_microarry_differential_gene.txt"), sep="\t", quote=F, row.names=F, col.names=T)
  tmp_df <- Frs%>%group_by(GeneName,SystematicName,Description,genename2)%>%summarise(adj.P.Val=min(adj.P.Val),P.Value=min(P.Value))%>%data.frame
  fdat <- inner_join(Frs,tmp_df,by=c("GeneName","SystematicName","Description","genename2","adj.P.Val","P.Value"))
  save(fdat,file=paste0(i,"/01_microarry_differential_gene_refine.Rdata"))
  save(normData,file=paste0(i,"/01_microarry_normdata.Rdata"))
  write.table(fdat, paste0(i,"/01_microarry_differential_gene_refine.txt"), sep="\t", quote=F, row.names=F, col.names=T)
  print(c(i,"finish"))
})
#######################################################  functional enrichment
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)

lapply(c("A_B","B_C"),function(i){
    j <- 0.005
    print(c(i,j,"start"))
    sub_dir <- paste0("/home/data/sdzl24/Project/aging/output_for_microarray/",i,"/Pvalue/",j)
    if(dir.exists(sub_dir)){
        print(c(sub_dir, "is exists\n"))
    }else{
        dir.create(sub_dir,recursive=TRUE)
    }
    setwd(sub_dir)
    load(paste0("/home/data/sdzl24/Project/aging/output_for_microarray/",i,"/01_microarry_differential_gene_refine.Rdata")) #fdat
    df <- fdat[,c("GeneName","SystematicName","log2FC","P.Value","adj.P.Val","genename2")]%>%unique
    df_tmp2 <-df%>%group_by(GeneName,SystematicName,genename2)%>%summarise(adj.P.Val=min(adj.P.Val),P.Value=min(P.Value))%>%data.frame
    df2 <- inner_join(df,df_tmp2,by=c("GeneName","SystematicName","genename2","adj.P.Val","P.Value"))
    df2$Change <- ifelse(df2$P.Value < j & abs(df2$log2FC) >= 0,ifelse(df2$log2FC > 0, 'Up', 'Down'), 'Stable')
    df2$PValue_cutoff <- j
    save(df2,file="01_microarry_gene_direction.Rdata")
    write.table(df2, "01_microarry_gene_direction.txt", sep="\t", quote=F, row.names=F, col.names=T)
})

library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

lapply(c("A_B","B_C"),function(i){
  j <-0.005
    print(c(i,j,"start"))
    sub_dir <- paste0("/home/data/sdzl24/Project/aging/output_for_microarray/",i,"/Pvalue/",j)
    setwd(sub_dir)
    load("01_microarry_gene_direction.Rdata") #df2
    down<- filter(df2,Change=="Down")
    up<- filter(df2,Change=="Up")
    save(down,file="01_microarry_gene_down.Rdata")
    save(up,file="01_microarry_gene_up.Rdata")
    write.table(down, "01_microarry_gene_down.txt", sep="\t", quote=F, row.names=F, col.names=T)
    write.table(up, "01_microarry_gene_up.txt", sep="\t", quote=F, row.names=F, col.names=T)
    lapply(c("down","up"),function(k){
      tmp_df <- read.csv(paste0("01_microarry_gene_",k,".txt"),header =T,sep = "\t")%>% data.frame
      colnames(tmp_df)[6] <-"SYMBOL"
      data = tmp_df[,"SYMBOL"]%>%unique()%>%as.vector()
      genes <- select(org.Hs.eg.db, keys=data, 
                    columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
      save(genes,file=paste0("01_",k,"_regulated_gene.Rdata"))
      write.table(genes,paste0("01_",k,"_regulated_gene.txt"),row.names = F, col.names = T,quote =F,sep="\t")
      write.xlsx(genes,paste0("01_",k,"_regulated_gene.xlsx"), rowNames = FALSE)
      ###############################################################
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
            write.xlsx(kegg1,paste0(k,"_KEGG.xlsx"), rowNames = FALSE)
            ###############################################################
            p1 <-ggplot(kegg1,aes( x = reorder(Description,Generatio1),y=Generatio1))+
                geom_segment(aes(yend=0, xend = Description)) +
                coord_flip()+
                geom_point(aes(color = Category, size = p.adjust)) + scale_color_manual(values=c("#8f85c3"))+
                labs(y="GeneRatio",x=NULL)+scale_size_continuous(range=c(10, 2))+
                guides(colour = guide_legend(override.aes = list(size=8), keyheight = 3 ),size = guide_legend(keyheight =2)) +
                P_theme+
                theme_cowplot(8)
            ggsave(paste0(k,"_KEGG.pdf"),p1, height=6,width=7)
        }

        types <-c("BP","CC","MF")
        #####################################################
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
        ######################################################
        ego1 <- do.call(rbind,tmp)
        if(is.null(ego1)){
            print(paste0("go no reslut"))
        }else{
          write.table(ego1,paste0(k,"_go.txt"),col.names=T,row.names =F,quote=F,sep="\t")
          write.xlsx(ego1,paste0(k,"_go.xlsx"), rowNames = FALSE)
      }
      print(c(i,j,k))
    })
  })
})
