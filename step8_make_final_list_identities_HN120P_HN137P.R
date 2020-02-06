# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# load data
scMutInfo <- readRDS(file=paste0("results/scMutInfo_reduced.rds"))
load(file="results/total_allele_reads_matrix.RData")
load(file="results/alternative_allele_reads_matrix.RData")

# HN120P
time_point <- "HN120P"
curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
HN120P_identity_variants <- NULL
for(i in 1:ncol(alternative_allele_reads_matrix)) {
    curr <- length(which(alternative_allele_reads_matrix[curr_cells,i]>=3))
    HN120P_identity_variants <- c(HN120P_identity_variants,curr)
}

# HN137P
time_point <- "HN137P"
curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
HN137P_identity_variants <- NULL
for(i in 1:ncol(alternative_allele_reads_matrix)) {
    curr <- length(which(alternative_allele_reads_matrix[curr_cells,i]>=3))
    HN137P_identity_variants <- c(HN137P_identity_variants,curr)
}

# save results
cell_lines_identities <- cbind(HN120P_identity_variants,HN137P_identity_variants)
colnames(cell_lines_identities) <- c("HN120P","HN137P")
rownames(cell_lines_identities) <- colnames(alternative_allele_reads_matrix)
cell_lines_identities <- cell_lines_identities[-which(as.numeric(cell_lines_identities[,"HN120P"])>0&as.numeric(cell_lines_identities[,"HN137P"])>0),]
save(cell_lines_identities,file="results/cell_lines_identities.RData")
