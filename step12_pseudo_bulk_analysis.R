# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# setting
cov_thr <- 100
vaf_thr_min <- 0.25
vaf_thr_common <- 0.00

# load data
cellInfo <- readRDS("data_info/SRAfinalTable_full.rds")
scMutInfo <- readRDS(file=paste0("results/scMutInfo_reduced.rds"))
load("results/depth_matrix.RData")

# HN120P
HN120_tot <- NULL
HN120P_tot <- depth_matrix[as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN120P"),"Run"]),]
HN120_tot <- rbind(HN120_tot,colSums(HN120P_tot))
HN120PCR_tot <- depth_matrix[as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN120PCR"),"Run"]),]
HN120_tot <- rbind(HN120_tot,colSums(HN120PCR_tot))
HN120PCRDH_tot <- depth_matrix[as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN120PCRDH"),"Run"]),]
HN120_tot <- rbind(HN120_tot,colSums(HN120PCRDH_tot))

# HN137P
HN137_tot <- NULL
HN137P_tot <- depth_matrix[rownames(depth_matrix)[which(rownames(depth_matrix)%in%as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN137P"),"Run"]))],]
HN137_tot <- rbind(HN137_tot,colSums(HN137P_tot))
HN137PCR_tot <- depth_matrix[as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN137PCR"),"Run"]),]
HN137_tot <- rbind(HN137_tot,colSums(HN137PCR_tot))
HN137PCRDH_tot <- depth_matrix[as.character(cellInfo[which(as.character(cellInfo[,"source_name"])=="HN137PCRDH"),"Run"]),]
HN137_tot <- rbind(HN137_tot,colSums(HN137PCRDH_tot))

# select well covered regions
HN120_137 <- rbind(HN120_tot,HN137_tot)
min_coverage <- NULL
for(i in 1:ncol(HN120_tot)) {
    min_coverage <- c(min_coverage,min(HN120_137[,i]))
}
HN120_137 <- HN120_137[,which(min_coverage>=cov_thr)]
rownames(HN120_137) <- c("HN120P","HN120PCR","HN120PCRDH","HN137P","HN137PCR","HN137PCRDH")
HN120_137_tot_counts <- t(HN120_137)

# alt counts matrix
HN120_137_alt_counts <- array(0,dim(HN120_137_tot_counts))
rownames(HN120_137_alt_counts) <- rownames(HN120_137_tot_counts)
colnames(HN120_137_alt_counts) <- colnames(HN120_137_tot_counts)
variants <- gsub(" ","",paste0(as.character(scMutInfo[,"Chr"]),"_",as.character(scMutInfo[,"PosStart"])))
cont <- 0
for(i in rownames(HN120_137_alt_counts)) {
    curr <- scMutInfo[which(as.character(variants)==i),]
    for(j in as.character(unique(curr[,"TimePoint"]))) {
        HN120_137_alt_counts[i,j] <- sum(as.numeric(as.character(curr[which(as.character(curr[,"TimePoint"])==j),"depth_alt"])))
    }
    cont <- cont + 1
    cat(cont/nrow(HN120_137_alt_counts),"\n")
}

# pseudo bulk results
HN120_137_freq <- HN120_137_alt_counts / HN120_137_tot_counts

# remove variants below vaf_thr_min both in HN120P and HN137P
invalid <- which(as.numeric(HN120_137_freq[,"HN120P"])<vaf_thr_min&as.numeric(HN120_137_freq[,"HN137P"])<vaf_thr_min)
HN120_137_tot_counts <- HN120_137_tot_counts[-invalid,]
HN120_137_alt_counts <- HN120_137_alt_counts[-invalid,]
HN120_137_freq <- HN120_137_freq[-invalid,]

# remove variants above vaf_thr_common both in HN120P and HN137P
invalid <- which(as.numeric(HN120_137_freq[,"HN120P"])>vaf_thr_common&as.numeric(HN120_137_freq[,"HN137P"])>vaf_thr_common)
HN120_137_tot_counts <- HN120_137_tot_counts[-invalid,]
HN120_137_alt_counts <- HN120_137_alt_counts[-invalid,]
HN120_137_freq <- HN120_137_freq[-invalid,]

# order based on support
new_order <- order(HN120_137_freq[,"HN120P"],1-HN120_137_freq[,"HN137P"])
HN120_137_tot_counts <- HN120_137_tot_counts[new_order,]
HN120_137_alt_counts <- HN120_137_alt_counts[new_order,]
HN120_137_freq <- HN120_137_freq[new_order,]

# save the results
save(HN120_137_tot_counts,file="pseudo_bulk/HN120_137_tot_counts.RData")
save(HN120_137_alt_counts,file="pseudo_bulk/HN120_137_alt_counts.RData")
save(HN120_137_freq,file="pseudo_bulk/HN120_137_freq.RData")
