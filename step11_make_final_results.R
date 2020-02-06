# set working directory
baseDir = "/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/"
setwd(baseDir)

# load data
SRAfinalTable = readRDS("data_info/SRAfinalTable_full.rds")
load(file="final_outputs/HN120P_identity.RData")
load(file="final_outputs/HN137P_identity.RData")
load(file="final_matters_arising/total_allele_reads_matrix.RData")
load(file="final_matters_arising/alt_allele_reads_matrix.RData")
load(file="final_matters_arising/final_list_genomic_regions.RData")

# make final data
HN120P_137P_identities <- array(NA,c(2196,17))
rownames(HN120P_137P_identities) <- 1:nrow(HN120P_137P_identities)
colnames(HN120P_137P_identities) <- c("mutation_id","chr","posStart","posEnd","REF","ALT","rsID","variant","MAF","cell_identity","cell_line","total_cells","valid_cells","mutated_cells","mutated_frequency","average_alt_depth","average_total_depth")

# HN120P identity
muts <- final_list_genomic_regions[which(final_list_genomic_regions%in%unique(HN120P_identity[,c("mutation")]))]
cont <- 0
for(i in muts) {
      for(j in c("HN120P","HN120PCR","HN120PCRDH","HN120M","HN120MCR","HN120MCRDH","HN137P","HN137PCR","HN137PCRDH","HN137M","HN137MCR")) {
            if(j!="HN137P") {
                  curr <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j),"Run"])
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN120P_identity[which(HN120P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN120P"
                  HN120P_137P_identities[cont,"cell_line"] <- j
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
            }
            else {
                  curr_run_single_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="SINGLE"),"Run"])
                  curr_run_single_end <- sort(c(rownames(total_allele_reads_matrix)[1:80],curr_run_single_end[which(curr_run_single_end%in%rownames(total_allele_reads_matrix))]))
                  curr_run_paired_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="PAIRED"),"Run"])
                  # single-end
                  curr <- curr_run_single_end
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN120P_identity[which(HN120P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN120P"
                  HN120P_137P_identities[cont,"cell_line"] <- "HN137P_single_end"
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  # paired-end
                  curr <- curr_run_paired_end
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN120P_identity[which(HN120P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN120P"
                  HN120P_137P_identities[cont,"cell_line"] <- "HN137P_paired_end"
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
            }
      }
}

# HN137P identity
muts <- final_list_genomic_regions[which(final_list_genomic_regions%in%unique(HN137P_identity[,c("mutation")]))]
for(i in muts) {
      for(j in c("HN120P","HN120PCR","HN120PCRDH","HN120M","HN120MCR","HN120MCRDH","HN137P","HN137PCR","HN137PCRDH","HN137M","HN137MCR")) {
            if(j!="HN137P") {
                  curr <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j),"Run"])
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN137P_identity[which(HN137P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN137P"
                  HN120P_137P_identities[cont,"cell_line"] <- j
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
            }
            else {
                  curr_run_single_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="SINGLE"),"Run"])
                  curr_run_single_end <- sort(c(rownames(total_allele_reads_matrix)[1:80],curr_run_single_end[which(curr_run_single_end%in%rownames(total_allele_reads_matrix))]))
                  curr_run_paired_end <- as.character(SRAfinalTable[which(as.character(SRAfinalTable[,"source_name"])==j&as.character(SRAfinalTable[,"LibraryLayout"])=="PAIRED"),"Run"])
                  # single-end
                  curr <- curr_run_single_end
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN137P_identity[which(HN137P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN137P"
                  HN120P_137P_identities[cont,"cell_line"] <- "HN137P_single_end"
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  # paired-end
                  curr <- curr_run_paired_end
                  curr_total_allele_reads_matrix <- total_allele_reads_matrix[curr,i]
                  curr_alt_allele_reads_matrix <- alt_allele_reads_matrix[curr,i]
                  cont <- cont + 1
                  HN120P_137P_identities[cont,"mutation_id"] <- i
                  HN120P_137P_identities[cont,c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")] <- as.character(HN137P_identity[which(HN137P_identity[,"mutation"]==i)[1],c("chr","posStart","posEnd","REF","ALT","rsID","variant","MAF")])
                  HN120P_137P_identities[cont,"cell_identity"] <- "HN137P"
                  HN120P_137P_identities[cont,"cell_line"] <- "HN137P_paired_end"
                  HN120P_137P_identities[cont,"total_cells"] <- length(curr_total_allele_reads_matrix)
                  HN120P_137P_identities[cont,"valid_cells"] <- length(which(curr_total_allele_reads_matrix>=5))
                  HN120P_137P_identities[cont,"mutated_cells"] <- length(which(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))]>=3))
                  HN120P_137P_identities[cont,"mutated_frequency"] <- as.numeric(HN120P_137P_identities[cont,"mutated_cells"]) / as.numeric(HN120P_137P_identities[cont,"total_cells"])
                  HN120P_137P_identities[cont,"average_alt_depth"] <- mean(curr_alt_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
                  HN120P_137P_identities[cont,"average_total_depth"] <- mean(curr_total_allele_reads_matrix[names(which(curr_total_allele_reads_matrix>=5))])
            }
      }
}

# save results
save(HN120P_137P_identities,file="final_matters_arising/HN120P_137P_identities.RData")
