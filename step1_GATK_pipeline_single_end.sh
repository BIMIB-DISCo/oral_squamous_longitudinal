# PIPELINE FROM FASTQ to ANNOVAR
# 
# REQUIRED EXTERNAL TOOLS:
#    STAR 2.7.3a
#    trimmomatic 0.39
#    samtools 1.9
#    gatk3 3.8-1-0-gf15c1c3ef
#    picard 2.21.3
#    annovar

# SETTING: base directory
basedir="/data_directory/BimiB/share/oral_squamous_longitudinal/"

# SETTING ANNOVAR: to be downloaded from the tool Website
annovar_c2a_cmd="perl /data_directory/BimiB/share/tools/annovar/convert2annovar.pl"
annovar_av_cmd="perl /data_directory/BimiB/share/tools/annovar/annotate_variation.pl"

# SETTING DIRECTORIES
fastqDir=${basedir}fastq/
fastqTrimDir=${basedir}trimmeds/
genomeDir_first_pass=${basedir}reference/genome_index/
genomeDir=${basedir}reference/second_hand_genome_index/
runDir_first_pass=${basedir}STAR_pass1/
runDir=${basedir}STAR_pass2/
singleSamDir=${basedir}singleSAM/
singleBamDir=${basedir}singleBAM/
picardDir=${basedir}picard/
vcfDir=${basedir}vcf/
annovarDir=${basedir}annotated/

# SETTING REFERENCE FILES
dbSNP_vcf=${basedir}reference/1000GENOMES-phase_3.vcf
genomeFa=${basedir}reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
annovar_db=${basedir}reference/annovar/RefGene_hg38/

# IMPORT FASTQs
fastqlist="$(find "$fastqDir" -type f  -exec basename {} \; | sort | cut -f 1 -d '.' | tr '\n' ' ')"

jobs=20

#===B0===

echo "trimming fastq"

mkdir -p $fastqTrimDir
if [ ! -d "$fastqTrimDir" ]; then
    echo "Error mkdir"
    exit 1
fi

# Default parameters. No illumina adapters
parallel -j $jobs trimmomatic SE -phred33 -threads 3 -quiet ${fastqDir}{}.fastq.gz ${fastqTrimDir}{}.trim.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -summary ${fastqTrimDir}{}.trim.summary ::: $fastqlist

#===E0===



#===B1===

echo "STAR -- generate genome indexing file"

mkdir -p $genomeDir_first_pass
if [ ! -d "$genomeDir_first_pass" ]; then
    echo "Error mkdir"
    exit 1
fi

STAR --runMode genomeGenerate --genomeDir $genomeDir_first_pass --genomeFastaFiles $genomeFa --runThreadN 50

#===E1===



#===B2===

echo "STAR -- first pass"

filelist="$(find "$fastqTrimDir" -type f -name "*.fastq.gz"| sort | tr '\n' ','| sed 's/,$//g')"


mkdir -p $runDir_first_pass
if [ ! -d "$runDir_first_pass" ]; then
    echo "Error mkdir"
    exit 1
fi

cd $runDir_first_pass

STAR --runMode alignReads --genomeDir $genomeDir_first_pass --readFilesCommand zcat --readFilesIn $filelist --runThreadN 50

cd $basedir

#===E2===



#===B3===

echo "STAR -- regenerate genome indexing file"

mkdir -p $genomeDir

if [ ! -d "$genomeDir" ]; then
    echo "Error mkdir"
    exit 1
fi

STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFa --sjdbFileChrStartEnd "${runDir_first_pass}/SJ.out.tab" --sjdbOverhang 75 --runThreadN 50

#===E3===



#===B4===

echo "STAR -- second pass"

mkdir -p $runDir
if [ ! -d "$runDir" ]; then
    echo "Error mkdir"
    exit 1
fi

cd $runDir

ulimit -n 2048

for file in $fastqlist; do
echo $(basename $file)
STAR --runMode alignReads --genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 --genomeDir $genomeDir --readFilesCommand zcat --readFilesIn "${fastqTrimDir}${file}.trim.fastq" --runThreadN 25 --outFileNamePrefix $file --outSAMtype BAM SortedByCoordinate

done

cd $basedir

#===E4===


# REMOVE GENOME FROM SHARED MEMORY
STAR --runMode alignReads --genomeLoad Remove --genomeDir $genomeDir


#===B=samtools===

echo "indexing fasta"

#indexing and dictionary creation of the fasta file inside 'reference' folder
samtools faidx $genomeFa
picard CreateSequenceDictionary R=$genomeFa

echo "The fasta file is indexed."

#===E=samtools===


# @2. Add read groups, sort, mark duplicates, and create index
# The above step produces a SAM file, which we then put through the usual Picard processing steps: adding read group information, sorting, marking duplicates and indexing. 


#===B=Picard===

echo "add or replace read groups"

mkdir -p $picardDir
if [ ! -d "$picardDir" ]; then
    echo "Error mkdir"
    exit 1
fi

# GATK BEST PRACTICE
parallel -j 40 picard AddOrReplaceReadGroups I="${runDir}{}Aligned.sortedByCoord.out.bam" O="${picardDir}{}.arrg_s1.bam" SO=coordinate RGID=0 RGPL=ILLUMINA RGLB=lib1 RGPU=group RGSM=cell ::: $fastqlist
parallel -j 20 picard MarkDuplicates I="${picardDir}{}.arrg_s1.bam" O="${picardDir}{}.md_s2.bam" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M="${picardDir}output_metrics_{}.txt" ::: $fastqlist
parallel -j 50 samtools index -b "${picardDir}{}.md_s2.bam" ::: $fastqlist
parallel -j 20 gatk3 -Xmx32G -T SplitNCigarReads -R ${genomeFa} -I "${picardDir}{}.md_s2.bam" -o "${picardDir}{}.snc_s3.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS ::: $fastqlist
parallel -j 20 gatk3 -Xmx32G -T BaseRecalibrator -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -o "${picardDir}{}.table" -knownSites ${dbSNP_vcf} ::: $fastqlist
parallel -j 40 gatk3 -Xmx3G -T PrintReads -kpr -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -BQSR "${picardDir}{}.table" -o "${picardDir}{}.bbq_s4.bam" ::: $fastqlist

echo "PICARD DONE."

#===E===



#===B=HaplotypeCaller===

echo "HaplotypeCaller and Filtration"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
    echo "Error mkdir"
    exit 1
fi

parallel -j 40 gatk3 -T HaplotypeCaller -R ${genomeFa} -I "${picardDir}{}.bbq_s4.bam" -o "${vcfDir}{}.raw.snps.indels.vcf" -dontUseSoftClippedBases -stand_call_conf 20 --dbsnp ${dbSNP_vcf} ::: $fastqlist
parallel -j 40 gatk3 -T VariantFiltration -R ${genomeFa} -V "${vcfDir}{}.raw.snps.indels.vcf" -window 35 -cluster 3 -filterName "FS" -filter \"FS '>' 30.0\" -filterName "QD" -filter \"QD '<' 2.0\" -o "${vcfDir}{}.filtered.vcf" ::: $fastqlist

#===E=HaplotypeCaller===



#===B=Annovar===

echo "Annotation with annovar"

mkdir -p $annovarDir
if [ ! -d "$annovarDir" ]; then
    echo "Error mkdir"
    exit 1
fi

parallel -j 40 "${annovar_c2a_cmd}" -format vcf4 "${vcfDir}{}.filtered.vcf" --outfile "${annovarDir}{}.anninput" --includeinfo ::: $fastqlist
parallel -j 40 "${annovar_av_cmd}" -out ${annovarDir}{} -exonicsplicing -build hg38 "${annovarDir}{}.anninput" "${annovar_db}" ::: $fastqlist

#===E=Annovar===
