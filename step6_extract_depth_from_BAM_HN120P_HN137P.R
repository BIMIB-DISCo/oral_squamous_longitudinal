##################################################################
# GET DEPTH FOR EACH SELECTED MUTATION FROM BAMs
##################################################################

# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/"
setwd(baseDir)

# load required libraries
library("foreach")
library("doParallel")

# load data
cellInfo <- readRDS("scripts/data_info/SRAfinalTable.rds")
scMutInfo <- readRDS(file=paste0("scripts/results/scMutInfo_reduced.rds"))
scMutInfo <- as.data.frame(scMutInfo,stringsAsFactors=FALSE)

# get list of variants occurring in at least 1 cell across all time points
listMut <- strsplit(scMutInfo$UIDsnp, "_")
scMutInfo$Mut_pos <- unlist(lapply(listMut, function(m) paste0(m[2], "_", m[3])))
listMutFormatted <- unique(scMutInfo$Mut_pos)
listMutFormatted <- unlist(lapply(strsplit(listMutFormatted, "_"), function(m) paste0(m[1], ":", m[2], "-", m[2])))

BAMlist1 <- paste(dir(path = paste0(baseDir, "HN120Primary/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ")

all_bam <- dir(path = paste0(baseDir, "HN137Primary/picard"), pattern = "s4.bam$", full.names = TRUE)
valid_bam <- NULL
for(i in all_bam) {
    curr <- strsplit(i,"/")[[1]]
    curr <- gsub(".bbq_s4.bam","",curr[length(curr)])
    if(length(which(as.character(cellInfo$Run)==curr))>0) {
        valid_bam <- c(valid_bam,i)
    }
}
BAMlist2 <- paste(valid_bam, collapse=" ")

BAMlist <- paste0(BAMlist1," ",BAMlist2)

numCores <- 40
registerDoParallel(numCores)

res <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

# save results
extracted_coverage <- res
save(extracted_coverage,file="scripts/results/extracted_coverage.RData")

##################################################################
# MAKE DEPTH MATRIX
##################################################################

names_row <- NULL
for(i in listMutFormatted) {
    names_row <- c(names_row,gsub(":","_",strsplit(i,"-")[[1]][1]))
}
names_col <- gsub(".bbq_s4.bam","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN120Primary/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Primary/picard/","",BAMlist)))
names_col <- strsplit(names_col," ")[[1]]

depth_tbl <- array(NA,c(length(names_row),length(names_col)))
rownames(depth_tbl) <- names_row
colnames(depth_tbl) <- names_col

for(i in 1:nrow(depth_tbl)) {
    curr <- strsplit(res[i],"\t")[[1]]
    curr <- as.numeric(curr[3:length(curr)])
    for(j in 1:ncol(depth_tbl)) {
        depth_tbl[i,j] <- curr[j]
    }
}
depth_tbl <- t(depth_tbl)
depth_tbl <- depth_tbl[sort(rownames(depth_tbl)),]
depth_tbl <- depth_tbl[,sort(colnames(depth_tbl))]

depth_matrix <- depth_tbl
save(depth_matrix,file="scripts/results/depth_matrix.RData")
