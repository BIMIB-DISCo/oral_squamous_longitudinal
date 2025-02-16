# set working directory
setwd("/data_directory/BimiB/share/oral_squamous_longitudinal/scripts/")

# get list of files to be considered and cells information
cellInfo <- readRDS("data_info/SRAfinalTable.rds")
listSCHN120Primary <- dir(path="/data_directory/BimiB/share/oral_squamous_longitudinal/HN120Primary/annotated/",pattern=".variant_function",full.names=TRUE)
listSCHN120Primary <- listSCHN120Primary[-grep("exonic_variant_function",listSCHN120Primary)]
listSCHN137Primary <- dir(path="/data_directory/BimiB/share/oral_squamous_longitudinal/HN137Primary/annotated/",pattern=".variant_function",full.names=TRUE)
listSCHN137Primary <- listSCHN137Primary[-grep("exonic_variant_function",listSCHN137Primary)]
listSCHN120Metast <- dir(path="/data_directory/BimiB/share/oral_squamous_longitudinal/HN120Metast/annotated/",pattern=".variant_function",full.names=TRUE)
listSCHN120Metast <- listSCHN120Metast[-grep("exonic_variant_function",listSCHN120Metast)]
listSC <- c(listSCHN120Primary,listSCHN137Primary,listSCHN120Metast)

# create data structure to save the results
tbFinal <- setNames(data.frame(matrix(ncol=13,nrow=0)), c("scID","TimePoint","Run","GSM","Gene","Chr","PosStart","PosEnd","REF","ALT","rs_ID","depth_tot","depth_alt"))
cn <- colnames(tbFinal)

# consider each cell
cont <- 0
for(sc in listSC) {
    
    # get current cell name
    run <- strsplit(x = sc, split = "\\.")[[1]][1]
    run <- strsplit(x = run, split = "SRR")[[1]][2]
    run <- paste0("SRR",run)

    if(length(which(cellInfo$Run==run))==1) {
    
        # read vcf for current cell
        tmpTbl <- read.table(file = sc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        tmpTbl$Run <- run
        
        # get information about source_name of experiment
        tmpTimePoint <- as.character(cellInfo$source_name[cellInfo$Run==run])
        tmpTbl$TimePoint <- tmpTimePoint
        tmpRHH <- as.character(cellInfo$RHH[cellInfo$Run==run])
        tmpTbl$scID <- tmpRHH
        tmpGSM <- as.character(cellInfo$Sample.Name[cellInfo$Run==run])
        tmpTbl$GSM <- tmpGSM
        
        # get genes names
        tmpTbl$V2 <- sapply(as.character(tmpTbl$V2), function(x){strsplit(x, "\\(")[[1]][1]})
        
        # get the remaining information
        OtherInfo <- tmpTbl[, "V17"]
        tmpTbl <- tmpTbl[, c("scID", "TimePoint", "Run", "GSM", "V2", "V3", "V4", "V5", "V6", "V7", "V10")]
        colnames(tmpTbl) <- c("scID", "TimePoint", "Run", "GSM", "Gene", "Chr", "PStart", "PEnd", "Ref", "Alt", "rs_ID")
        
        # extract quality and coverage information
        for(m in 1:length(OtherInfo)) {
            
            info <- unlist(strsplit(OtherInfo[m],":"))
            
            if(length(info) < 5) {
                depth_tot <- NA
                depth_alt <- NA
            }
            else {
                RefAlt <- unlist(strsplit(info[2], ","))
                depth_tot <- sum(as.numeric(RefAlt))
                depth_alt <- as.numeric(RefAlt[2])
            }
            tmpTbl$depth_tot[m] <- depth_tot
            tmpTbl$depth_alt[m] <- depth_alt
            
        }
        
        # save results for current cell
        colnames(tmpTbl) <- cn
        tbFinal <- rbind(tbFinal,tmpTbl)
    }
    
    cont = cont + 1
    cat(paste0(cont, " done | ", round(cont/length(listSC),digits = 2),"%\r"))
    
}

tbFinal$UIDsnp <- paste0(tbFinal$Gene,"_",tbFinal$Chr,"_",tbFinal$PosStart,"_",tbFinal$ALT)
saveRDS(tbFinal,file=paste0("results/cells_aggregate_info.rds"))
