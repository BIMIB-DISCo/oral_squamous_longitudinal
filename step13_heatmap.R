# set working directory
baseDir <- "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# load the libraries
library("gplots")
library("RColorBrewer")

# load data
SRAfinalTable <- readRDS("data_info/SRAfinalTable_full.rds")
load(file="final_matters_arising/total_allele_reads_matrix.RData")
load(file="final_matters_arising/alt_allele_reads_matrix.RData")

# clusters labels
clusters_labels <- array(NA,c(dim(alt_allele_reads_matrix)[1],1))
rownames(clusters_labels) <- rownames(total_allele_reads_matrix)
colnames(clusters_labels) <- "cluster"
for(j in c("HN120P","HN120PCR","HN120PCRDH","HN120M","HN120MCR","HN120MCRDH","HN137P","HN137PCR","HN137PCRDH","HN137M","HN137MCR")) {
    if(j!="HN137P") {
        curr <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j),"Run"])
        clusters_labels[curr,"cluster"] <- j
    }
    else {
        curr_run_single_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="SINGLE"),"Run"])
        curr_run_single_end <- sort(c(rownames(total_allele_reads_matrix)[1:80],curr_run_single_end[which(curr_run_single_end%in%rownames(total_allele_reads_matrix))]))
        curr_run_paired_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="PAIRED"),"Run"])
        clusters_labels[curr_run_single_end,"cluster"] <- "HN137P (single-end)"
        clusters_labels[curr_run_paired_end,"cluster"] <- "HN137P (paired-end)"
    }
}

# make data object
invalid_variants <- c("3_49359890","15_44711516","17_81510318","18_3277865")
valid_variants <- 
data <- alt_allele_reads_matrix[rownames(clusters_labels),]
data <- data[,colnames(alt_allele_reads_matrix)[which(!colnames(alt_allele_reads_matrix)%in%invalid_variants)]]
data[which(data<3)] <- 0
data[which(data>0)] <- 1

# order data based on clustering labels
ordered_labels <- c("HN120P","HN120PCR","HN120PCRDH","HN120M","HN120MCR","HN120MCRDH","HN137P (single-end)","HN137P (paired-end)","HN137PCR","HN137PCRDH","HN137M","HN137MCR")
ordered_idx <- NULL
ordered_idx_labels <- NULL
cont <- 0
for(i in ordered_labels) {
    ordered_idx <- c(ordered_idx,as.numeric(which(clusters_labels[,"cluster"]==i)))
    cont <- cont + 1
    ordered_idx_labels <- c(ordered_idx_labels,rep(cont,length(as.numeric(which(clusters_labels[,"cluster"]==i)))))
    
}
data <- data[ordered_idx,]

# make heatmap
my_palette = colorRampPalette(c("blue","white","red"))(n = 299)
mat_data = data.matrix(t(data))
heatmap.2(mat_data,main="Heatmap",notecol="black",density.info="none",trace="none",margins=c(12,13),col=my_palette,dendrogram="row",Rowv=TRUE,Colv=FALSE,ColSideColors=as.character(ordered_idx_labels))
legend("topright",legend=ordered_labels,col=1:12,lty=1,lwd=5,cex=.4)
