# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# load required libraries
library("biomaRt")

# ensembl to get known rsIDs
ensembl_snp <- useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp")

# load data
scMutInfo <- readRDS(file=paste0("results/scMutInfo_reduced.rds"))
positions <- gsub(" ","",paste0(as.character(scMutInfo[,"Chr"]),"_",as.character(scMutInfo[,"PosStart"])))
load(file="results/total_allele_reads_matrix.RData")
load(file="results/alternative_allele_reads_matrix.RData")
load(file="results/cell_lines_identities.RData")

# select counts only for variants of interest
alternative_allele_reads_matrix <- alternative_allele_reads_matrix[,sort(unique(rownames(cell_lines_identities)))]
total_allele_reads_matrix <- total_allele_reads_matrix[,sort(unique(rownames(cell_lines_identities)))]

# HN120P identity
HN120P_identity <- NULL
valid_variants_info <- unique(scMutInfo[which(positions%in%sort(unique(names(which(cell_lines_identities[,"HN120P"]>0))))),c("Chr","PosStart","PosEnd","REF","ALT")])
for(i in sort(unique(names(which(cell_lines_identities[,"HN120P"]>0))))) {
    for(j in c("HN120P","HN137PCR","HN137PCRDH","HN137P","HN120PCR","HN120PCRDH")) {
        curr <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==j),] # consider current cell line
        curr_cells <- sort(unique(as.character(curr[,"Run"])))
        mutation <- i # consider current variant
        chr <- strsplit(mutation,"_")[[1]][1]
        posStart <- strsplit(mutation,"_")[[1]][2]
        posEnd <- posStart
        REF <- as.character(valid_variants_info[which(as.character(valid_variants_info[,"Chr"])==chr&gsub(" ","",as.character(valid_variants_info[,"PosStart"]))==posStart),"REF"])
        ALT <- as.character(valid_variants_info[which(as.character(valid_variants_info[,"Chr"])==chr&gsub(" ","",as.character(valid_variants_info[,"PosStart"]))==posStart),"ALT"])
        cell_line <- j
        total_cells <- length(curr_cells) # number of cells in current cell line
        valid_cells <- length(which(total_allele_reads_matrix[curr_cells,mutation]>=5)) # number of cells covering the considered variant
        mutated_cells <- length(which(alternative_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]>=3)) # number of cells with the variant among the well covered ones
        mutated_frequency <- mutated_cells / valid_cells
        average_alt_depth <- mean(alternative_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]) # mean alt depth among the valid cells
        average_total_depth <- mean(total_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]) # mean total depth among the valid cells
        curr_res <- c(mutation,chr,posStart,posEnd,REF,ALT,NA,NA,NA,cell_line,total_cells,valid_cells,mutated_cells,mutated_frequency,average_alt_depth,average_total_depth)
        HN120P_identity <- rbind(HN120P_identity,curr_res)
    }
}
rownames(HN120P_identity) <- 1:nrow(HN120P_identity)
colnames(HN120P_identity) <- c("mutation","chr","posStart","posEnd","REF","ALT","rsID","variant","MAF","cell_line","total_cells","valid_cells","mutated_cells","mutated_frequency","average_alt_depth","average_total_depth")

# HN137P identity
HN137P_identity <- NULL
valid_variants_info <- unique(scMutInfo[which(positions%in%sort(unique(names(which(cell_lines_identities[,"HN137P"]>0))))),c("Chr","PosStart","PosEnd","REF","ALT")])
for(i in sort(unique(names(which(cell_lines_identities[,"HN137P"]>0))))) {
    for(j in c("HN120P","HN137PCR","HN137PCRDH","HN137P","HN120PCR","HN120PCRDH")) {
        curr <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])==j),] # consider current cell line
        curr_cells <- sort(unique(as.character(curr[,"Run"])))
        mutation <- i # consider current variant
        chr <- strsplit(mutation,"_")[[1]][1]
        posStart <- strsplit(mutation,"_")[[1]][2]
        posEnd <- posStart
        REF <- as.character(valid_variants_info[which(as.character(valid_variants_info[,"Chr"])==chr&gsub(" ","",as.character(valid_variants_info[,"PosStart"]))==posStart),"REF"])
        ALT <- as.character(valid_variants_info[which(as.character(valid_variants_info[,"Chr"])==chr&gsub(" ","",as.character(valid_variants_info[,"PosStart"]))==posStart),"ALT"])
        cell_line <- j
        total_cells <- length(curr_cells) # number of cells in current cell line
        valid_cells <- length(which(total_allele_reads_matrix[curr_cells,mutation]>=5)) # number of cells covering the considered variant
        mutated_cells <- length(which(alternative_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]>=3)) # number of cells with the variant among the well covered ones
        mutated_frequency <- mutated_cells / valid_cells
        average_alt_depth <- mean(alternative_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]) # mean alt depth among the valid cells
        average_total_depth <- mean(total_allele_reads_matrix[names(which(total_allele_reads_matrix[curr_cells,mutation]>=5)),mutation]) # mean total depth among the valid cells
        curr_res <- c(mutation,chr,posStart,posEnd,REF,ALT,NA,NA,NA,cell_line,total_cells,valid_cells,mutated_cells,mutated_frequency,average_alt_depth,average_total_depth)
        HN137P_identity <- rbind(HN137P_identity,curr_res)
    }
}
rownames(HN137P_identity) <- 1:nrow(HN137P_identity)
colnames(HN137P_identity) <- c("mutation","chr","posStart","posEnd","REF","ALT","rsID","variant","MAF","cell_line","total_cells","valid_cells","mutated_cells","mutated_frequency","average_alt_depth","average_total_depth")

# rsIDs for HN120P_identity
variants <- unique(HN120P_identity[,c("chr","posStart","posEnd","REF","ALT")])
rsIDs_HN120P_identity <- list()
for(i in 1:nrow(variants)) {
    temp <- NULL
    temp <- getBM(attributes=c("refsnp_id","allele","minor_allele","minor_allele_freq"),filters=c("chr_name","start","end"),values=list(as.character(variants[i,"chr"]),as.character(variants[i,"posStart"]),as.character(variants[i,"posEnd"])),mart=ensembl_snp)
    rsIDs_HN120P_identity[[i]] <- temp
    cat(i/nrow(variants),"\n")
}
rsIDs_HN120P_variants <- variants

# rsIDs for HN137P_identity
variants <- unique(HN137P_identity[,c("chr","posStart","posEnd","REF","ALT")])
rsIDs_HN137P_identity <- list()
for(i in 1:nrow(variants)) {
    temp <- NULL
    temp <- getBM(attributes=c("refsnp_id","allele","minor_allele","minor_allele_freq"),filters=c("chr_name","start","end"),values=list(as.character(variants[i,"chr"]),as.character(variants[i,"posStart"]),as.character(variants[i,"posEnd"])),mart=ensembl_snp)
    rsIDs_HN137P_identity[[i]] <- temp
    cat(i/nrow(variants),"\n")
}
rsIDs_HN137P_variants <- variants

# manually verify rsIDs

# rsIDs for HN120P_identity
rsIDs_HN120P_identity[[1]] <- rsIDs_HN120P_identity[[1]][8,]
rsIDs_HN120P_identity[[2]] <- rsIDs_HN120P_identity[[2]][2,]
rsIDs_HN120P_identity[[4]] <- rsIDs_HN120P_identity[[4]][2,]
rsIDs_HN120P_identity[[5]] <- rsIDs_HN120P_identity[[5]][1,]
rsIDs_HN120P_identity[[8]] <- rsIDs_HN120P_identity[[8]][2,]
rsIDs_HN120P_identity[[11]] <- rsIDs_HN120P_identity[[11]][3,]
rsIDs_HN120P_identity[[18]] <- rsIDs_HN120P_identity[[18]][1,]
rsIDs_HN120P_identity[[23]] <- rsIDs_HN120P_identity[[23]][2,]
rsIDs_HN120P_identity[[26]] <- rsIDs_HN120P_identity[[26]][1,]
rsIDs_HN120P_identity[[27]] <- rsIDs_HN120P_identity[[27]][1,]
rsIDs_HN120P_identity[[30]] <- NA
rsIDs_HN120P_identity[[37]] <- rsIDs_HN120P_identity[[37]][2,]
rsIDs_HN120P_identity[[38]] <- rsIDs_HN120P_identity[[38]][2,]
rsIDs_HN120P_identity[[40]] <- rsIDs_HN120P_identity[[40]][2,]
rsIDs_HN120P_identity[[47]] <- rsIDs_HN120P_identity[[47]][2,]
rsIDs_HN120P_identity[[54]] <- rsIDs_HN120P_identity[[54]][2,]
rsIDs_HN120P_identity[[55]] <- rsIDs_HN120P_identity[[55]][2,]
rsIDs_HN120P_identity[[58]] <- rsIDs_HN120P_identity[[58]][2,]
rsIDs_HN120P_identity[[67]] <- rsIDs_HN120P_identity[[67]][1,]

# rsIDs for HN137P_identity
rsIDs_HN137P_identity[[1]] <- rsIDs_HN137P_identity[[1]][1,]
rsIDs_HN137P_identity[[14]] <- rsIDs_HN137P_identity[[14]][1,]
rsIDs_HN137P_identity[[20]] <- rsIDs_HN137P_identity[[20]][1,]
rsIDs_HN137P_identity[[21]] <- rsIDs_HN137P_identity[[21]][2,]
rsIDs_HN137P_identity[[24]] <- rsIDs_HN137P_identity[[24]][2,]
rsIDs_HN137P_identity[[25]] <- rsIDs_HN137P_identity[[25]][2,]
rsIDs_HN137P_identity[[36]] <- rsIDs_HN137P_identity[[36]][2,]
rsIDs_HN137P_identity[[37]] <- rsIDs_HN137P_identity[[37]][3,]
rsIDs_HN137P_identity[[38]] <- NA
rsIDs_HN137P_identity[[41]] <- rsIDs_HN137P_identity[[41]][2,]
rsIDs_HN137P_identity[[51]] <- rsIDs_HN137P_identity[[51]][2,]
rsIDs_HN137P_identity[[65]] <- rsIDs_HN137P_identity[[65]][2,]
rsIDs_HN137P_identity[[68]] <- NA
rsIDs_HN137P_identity[[74]] <- rsIDs_HN137P_identity[[74]][3,]
rsIDs_HN137P_identity[[77]] <- NA
rsIDs_HN137P_identity[[79]] <- rsIDs_HN137P_identity[[79]][2,]
rsIDs_HN137P_identity[[82]] <- rsIDs_HN137P_identity[[82]][2,]
rsIDs_HN137P_identity[[93]] <- rsIDs_HN137P_identity[[93]][1,]
rsIDs_HN137P_identity[[94]] <- rsIDs_HN137P_identity[[94]][2,]
rsIDs_HN137P_identity[[101]] <- rsIDs_HN137P_identity[[101]][16,]
rsIDs_HN137P_identity[[106]] <- rsIDs_HN137P_identity[[106]][1,]
rsIDs_HN137P_identity[[107]] <- rsIDs_HN137P_identity[[107]][2,]
rsIDs_HN137P_identity[[111]] <- rsIDs_HN137P_identity[[111]][2,]

# make final results

# HN120P_identity
for(i in 1:length(rsIDs_HN120P_identity)) {
    if(!is.na(rsIDs_HN120P_identity[[i]])) {
        curr_pos <- as.numeric(which(HN120P_identity[,"chr"]==rsIDs_HN120P_variants[i,"chr"]&HN120P_identity[,"posStart"]==rsIDs_HN120P_variants[i,"posStart"]))
        HN120P_identity[curr_pos,"rsID"] <- as.character(rsIDs_HN120P_identity[[i]][,"refsnp_id"])
        HN120P_identity[curr_pos,"variant"] <- as.character(rsIDs_HN120P_identity[[i]][,"allele"])
        HN120P_identity[curr_pos,"MAF"] <- as.numeric(rsIDs_HN120P_identity[[i]][,"minor_allele_freq"])
    }
}

# HN137P_identity
for(i in 1:length(rsIDs_HN137P_identity)) {
    if(!is.na(rsIDs_HN137P_identity[[i]])) {
        curr_pos <- as.numeric(which(HN137P_identity[,"chr"]==rsIDs_HN137P_variants[i,"chr"]&HN137P_identity[,"posStart"]==rsIDs_HN137P_variants[i,"posStart"]))
        HN137P_identity[curr_pos,"rsID"] <- as.character(rsIDs_HN137P_identity[[i]][,"refsnp_id"])
        HN137P_identity[curr_pos,"variant"] <- as.character(rsIDs_HN137P_identity[[i]][,"allele"])
        HN137P_identity[curr_pos,"MAF"] <- as.numeric(rsIDs_HN137P_identity[[i]][,"minor_allele_freq"])
    }
}

# save the results
save(alternative_allele_reads_matrix,file="final_outputs/alternative_allele_reads_matrix.RData")
save(total_allele_reads_matrix,file="final_outputs/total_allele_reads_matrix.RData")
save(HN120P_identity,file="final_outputs/HN120P_identity.RData")
save(HN137P_identity,file="final_outputs/HN137P_identity.RData")
