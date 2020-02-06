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
load("scripts/final_matters_arising/final_list_genomic_regions.RData")

# get list of final selected variants
listMut <- final_list_genomic_regions
listMutFormatted <- unlist(lapply(strsplit(listMut, "_"), function(m) paste0(m[1], ":", m[2], "-", m[2])))

BAMlist1 <- paste(dir(path = paste0(baseDir, "HN120Primary/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN120Primary

all_bam <- dir(path = paste0(baseDir, "HN137Primary/picard"), pattern = "s4.bam$", full.names = TRUE)
valid_bam <- NULL
for(i in all_bam) {
    curr <- strsplit(i,"/")[[1]]
    curr <- gsub(".bbq_s4.bam","",curr[length(curr)])
    if(length(which(as.character(cellInfo$Run)==curr))>0) {
        valid_bam <- c(valid_bam,i)
    }
}
BAMlist2 <- paste(valid_bam, collapse=" ") # HN137Primary

BAMlist3 <- paste(dir(path = paste0(baseDir, "HN120Metast/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN120Metast

BAMlist4 <- paste(dir(path = paste0(baseDir, "HN137MCR/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN137MCR

BAMlist5 <- paste(dir(path = paste0(baseDir, "HN137Metast_PAIRED/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN137Metast_PAIRED

BAMlist6 <- paste(dir(path = paste0(baseDir, "HN137Primary_DUPLICATED/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN137Primary_DUPLICATED

BAMlist7 <- paste(dir(path = paste0(baseDir, "HN137Primary_PAIRED/picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ") # HN137Primary_PAIRED

BAMlist <- paste0(BAMlist1," ",BAMlist2," ",BAMlist3," ",BAMlist4," ",BAMlist5," ",BAMlist6," ",BAMlist7)

numCores <- 40
registerDoParallel(numCores)

res1 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist1), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res2 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist2), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res3 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist3), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res4 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist4), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res5 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist5), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res6 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist6), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res7 <- foreach(i = listMutFormatted, .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",i," ", BAMlist7), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

res <- list()
res[["bam1"]] <- res1
res[["bam2"]] <- res2
res[["bam3"]] <- res3
res[["bam4"]] <- res4
res[["bam5"]] <- res5
res[["bam6"]] <- res6
res[["bam7"]] <- res7

# save results
extracted_coverage <- res
save(extracted_coverage,file="scripts/final_matters_arising/extracted_coverage.RData")

##################################################################
# MAKE DEPTH MATRIX
##################################################################

names_row <- NULL
for(i in listMutFormatted) {
    names_row <- c(names_row,gsub(":","_",strsplit(i,"-")[[1]][1]))
}
names_col <- gsub(".bbq_s4.bam","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN120Primary/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Primary/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN120Metast/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137MCR/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Metast_PAIRED/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Primary_DUPLICATED/picard/","",gsub("/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Primary_PAIRED/picard/","",BAMlist))))))))
names_col <- strsplit(names_col," ")[[1]]

depth_tbl <- array(NA,c(length(names_row),length(names_col)))
rownames(depth_tbl) <- names_row
colnames(depth_tbl) <- names_col

for(i in 1:nrow(depth_tbl)) {
    curr <- strsplit(res[["bam1"]],"\t")[[1]]
    curr1 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam2"]],"\t")[[1]]
    curr2 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam3"]],"\t")[[1]]
    curr3 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam4"]],"\t")[[1]]
    curr4 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam5"]],"\t")[[1]]
    curr5 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam6"]],"\t")[[1]]
    curr6 <- as.numeric(curr[3:length(curr)])
    curr <- strsplit(res[["bam7"]],"\t")[[1]]
    curr7 <- as.numeric(curr[3:length(curr)])
    curr <- as.numeric(c(curr1,curr2,curr3,curr4,curr5,curr6,curr7))
    for(j in 1:ncol(depth_tbl)) {
        depth_tbl[i,j] <- curr[j]
    }
}
depth_tbl <- t(depth_tbl)

total_allele_reads_matrix <- depth_tbl
total_allele_reads_matrix <- total_allele_reads_matrix[sort(unique(rownames(total_allele_reads_matrix))),]
save(total_allele_reads_matrix,file="scripts/final_matters_arising/total_allele_reads_matrix.RData")

##################################################################
# MAKE ALT MATRIX
##################################################################

# process mutations data
variants_part1 <- readRDS("scripts/results/cells_aggregate_info.rds")
variants_part2 <- readRDS("scripts/results/cells_aggregate_info_part2.rds")
variants <- rbind(variants_part1,variants_part2)
variants_id <- gsub(" ","",paste0(as.character(variants[,"Chr"]),"_",as.character(variants[,"PosStart"])))
variants <- variants[which(variants_id%in%colnames(total_allele_reads_matrix)),]
variants <- variants[,c("Run","TimePoint","Chr","PosStart","PosEnd","REF","ALT","depth_alt")]
variants <- variants[which(as.character(variants[,"REF"])%in%c("A","C","G","T")),]
variants <- variants[which(as.character(variants[,"ALT"])%in%c("A","C","G","T")),]
invalid_entries <- c("1_173865494_173865494_G_T","2_73958749_73958749_A_G","2_237098197_237098197_A_T","3_40457354_40457354_C_A","12_124913493_124913493_G_A")
invalid_entries <- which(gsub(" ","",paste0(as.character(variants[,"Chr"]),"_",as.character(variants[,"PosStart"]),"_",as.character(variants[,"PosEnd"]),"_",as.character(variants[,"REF"]),"_",as.character(variants[,"ALT"])))%in%invalid_entries)
variants <- variants[-invalid_entries,]

# make final matrix
alt_allele_reads_matrix <- array(0,dim(total_allele_reads_matrix))
rownames(alt_allele_reads_matrix) <- rownames(total_allele_reads_matrix)
colnames(alt_allele_reads_matrix) <- colnames(total_allele_reads_matrix)
for(i in 1:nrow(variants)) {
    curr <- grep(as.character(variants[i,"Run"]),rownames(alt_allele_reads_matrix))
    alt_allele_reads_matrix[curr,gsub(" ","",paste0(as.character(variants[i,"Chr"]),"_",as.character(variants[i,"PosStart"])))] <- as.numeric(as.character(variants[i,"depth_alt"]))
}
save(alt_allele_reads_matrix,file="scripts/final_matters_arising/alt_allele_reads_matrix.RData")
