# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# load data
scMutInfo <- readRDS(file=paste0("results/scMutInfo_reduced.rds"))
load("results/depth_matrix.RData")

# remove positions with total read counts <5 in >50% of the cells in each time point
valid_HN120P <- NULL
valid_HN120PCR <- NULL
valid_HN120PCRDH <- NULL
valid_HN137P <- NULL
valid_HN137PCR <- NULL
valid_HN137PCRDH <- NULL
for(i in 1:ncol(depth_matrix)) {
    # HN120P
    time_point <- "HN120P"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN120P <- c(valid_HN120P,is.valid)
    # HN120PCR
    time_point <- "HN120PCR"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN120PCR <- c(valid_HN120PCR,is.valid)
    # HN120PCRDH
    time_point <- "HN120PCRDH"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN120PCRDH <- c(valid_HN120PCRDH,is.valid)
    # HN137P
    time_point <- "HN137P"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN137P <- c(valid_HN137P,is.valid)
    # HN137PCR
    time_point <- "HN137PCR"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN137PCR <- c(valid_HN137PCR,is.valid)
    # HN137PCRDH
    time_point <- "HN137PCRDH"
    curr_data <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==time_point),]
    curr_cells <- sort(unique(as.character(curr_data[,"Run"])))
    num_cells <- length(curr_cells)
    low_coverage_cells <- length(which(as.numeric(depth_matrix[curr_cells,i])<5))
    is.valid <- TRUE
    if((low_coverage_cells/num_cells)>0.50) {
        is.valid <- FALSE
    }
    valid_HN137PCRDH <- c(valid_HN137PCRDH,is.valid)
    cat(i/ncol(depth_matrix),"\n")
}
valid_positions <- which(valid_HN120P&valid_HN120PCR&valid_HN120PCRDH&valid_HN137P&valid_HN137PCR&valid_HN137PCRDH)
depth_matrix <- depth_matrix[,valid_positions]

# make variants matrix
variants <- array(0,dim(depth_matrix))
rownames(variants) <- rownames(depth_matrix)
colnames(variants) <- colnames(depth_matrix)
cont <- 0
for(i in rownames(variants)) {
    curr <- scMutInfo[which(as.character(scMutInfo[,"Run"])==i),]
    listVariants <- strsplit(as.character(curr[,"UIDsnp"]), "_")
    listVariants <- unlist(lapply(listVariants, function(m) paste0(m[2], "_", m[3])))
    listVariantsFormatted <- unlist(lapply(strsplit(listVariants, "_"), function(m) paste0(m[1], ":", m[2], "-", m[2])))
    listVariantsFormatted <- unique(listVariantsFormatted)
    names_row <- NULL
    for(j in listVariantsFormatted) {
        names_row <- c(names_row,gsub(":","_",strsplit(j,"-")[[1]][1]))
    }
    curr_valid <- names_row[which(names_row%in%colnames(variants))]
    for(j in curr_valid) {
        variants[i,j] <- mean(as.numeric(as.character(curr[,"depth_alt"]))[which(as.character(listVariants)==j)])
    }
    cont <- cont + 1
    cat(cont/nrow(variants),"\n")
}

# save final data
total_allele_reads_matrix <- depth_matrix
total_allele_reads_matrix <- total_allele_reads_matrix[sort(rownames(total_allele_reads_matrix)),]
total_allele_reads_matrix <- total_allele_reads_matrix[,sort(colnames(total_allele_reads_matrix))]
alternative_allele_reads_matrix <- variants
alternative_allele_reads_matrix <- alternative_allele_reads_matrix[sort(rownames(alternative_allele_reads_matrix)),]
alternative_allele_reads_matrix <- alternative_allele_reads_matrix[,sort(colnames(alternative_allele_reads_matrix))]
save(total_allele_reads_matrix,file="results/total_allele_reads_matrix.RData")
save(alternative_allele_reads_matrix,file="results/alternative_allele_reads_matrix.RData")
