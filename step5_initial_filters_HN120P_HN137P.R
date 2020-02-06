# set working directory
setwd("/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/")

# read data
scMutInfo <- readRDS(file=paste0("results/cells_aggregate_info.rds"))
scMutInfo <- as.matrix(scMutInfo)

# select cells of time points of interest
line120p <- c("HN120P","HN120PCR","HN120PCRDH")
line137p <- c("HN137P","HN137PCR","HN137PCRDH")
scMutInfo_reduced <- scMutInfo[which(as.character(scMutInfo[,"TimePoint"])%in%c(line120p,line137p)),]

# remove indels and other structural variants
scMutInfo_reduced <- scMutInfo_reduced[which(as.character(scMutInfo_reduced[,"REF"])%in%c("A","C","G","T")),]
scMutInfo_reduced <- scMutInfo_reduced[which(as.character(scMutInfo_reduced[,"ALT"])%in%c("A","C","G","T")),]

# remove variants mapped on mitochondrial genes
scMutInfo_reduced <- scMutInfo_reduced[which(as.character(scMutInfo_reduced[,"Chr"])%in%c(as.character(1:22),"X","Y")),]

# consider only variants observed at a frequency of at least 20% (at least 18 cells) in either HN120P or HN137P
# (this is done to remove variants likely due to artifacts or shared across the two cell lines)

# HN120P
HN120P_variants <- scMutInfo_reduced[which(as.character(scMutInfo_reduced[,"TimePoint"])=="HN120P"),]
HN120P_variants <- gsub(" ","",paste0(as.character(HN120P_variants[,"Chr"]),"_",as.character(HN120P_variants[,"PosStart"]),"_",as.character(HN120P_variants[,"REF"]),"_",as.character(HN120P_variants[,"ALT"])))
HN120P_variants <- sort(table(HN120P_variants),decreasing=TRUE)
HN120P_variants <- sort(unique(names(which(HN120P_variants>=18)))) # variants representative for HN120P identity

# HN137P
HN137P_variants <- scMutInfo_reduced[which(as.character(scMutInfo_reduced[,"TimePoint"])=="HN137P"),]
HN137P_variants <- gsub(" ","",paste0(as.character(HN137P_variants[,"Chr"]),"_",as.character(HN137P_variants[,"PosStart"]),"_",as.character(HN137P_variants[,"REF"]),"_",as.character(HN137P_variants[,"ALT"])))
HN137P_variants <- sort(table(HN137P_variants),decreasing=TRUE)
HN137P_variants <- sort(unique(names(which(HN137P_variants>=18)))) # variants representative for HN137P identity

# get final list of valid variants
all_variants <- sort(unique(c(HN120P_variants,HN137P_variants)))
invalid_variants <- intersect(HN120P_variants,HN137P_variants)
valid_variants <- sort(unique(all_variants[which(!all_variants%in%invalid_variants)]))
all_variants <- gsub(" ","",paste0(as.character(scMutInfo_reduced[,"Chr"]),"_",as.character(scMutInfo_reduced[,"PosStart"]),"_",as.character(scMutInfo_reduced[,"REF"]),"_",as.character(scMutInfo_reduced[,"ALT"])))
scMutInfo_reduced <- scMutInfo_reduced[which(all_variants%in%valid_variants),]

# save the results
saveRDS(scMutInfo_reduced,file=paste0("results/scMutInfo_reduced.rds"))
