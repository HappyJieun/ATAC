#################################
## 2023-12-10                                 ##
## ATAC_seq pipline code                   ##
## ATAC_pip.R                                 ##
## 2019270926                                ##    
## JIEUN LEE                                   ##  
#################################

## set Repositories
setRepositories(ind = 1:7)

## Insert project name
PROJECT <- "PRJNA756705"
species <- "gallus_gallus"

## sample download
# setwd(paste0(getwd(), "/", sample_dir))
# system(paste(paste0(getwd(), "/", "download_via_parallel.sh"), paste0(getwd(), "/", "accession_list.txt")))
# setwd(base_dir)

## set. work dir
WORK_DIR <- paste0("/disk7/bilje/ATAC_seq/", PROJECT)
setwd(WORK_DIR)

## check sample
Sample_DIR <- paste0(WORK_DIR, "/0.sample/")

sample <- list.files(path=Sample_DIR, pattern = ".fastq.gz")
# length(sample)

## sample name
sample_names <- c()
for (sample in sample) {
  name <- strsplit(sample, "_")[[1]][1]  
  sample_names <- c(sample_names, name)  
}
sample_names <- unique(sample_names)
length(sample_names)


## 1. trimmomatic - cut adapt  
# make dir 
Trim_DIR <- paste0(WORK_DIR, "/1.trim/")
system(paste0("mkdir -p ", Trim_DIR))
system(paste0("mkdir -p ", Trim_DIR, "log/"))

trim <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("cutadapt",
                paste0(Sample_DIR, sample_name, "_1.fastq.gz"),
                paste0(Sample_DIR, sample_name, "_2.fastq.gz"),
                paste0("-o ", Trim_DIR, sample_name, "_1.fastq.gz"),
                paste0("-p ", Trim_DIR, sample_name, "_2.fastq.gz"),
                "-l 50 --cores=16",
                paste0("> ", Trim_DIR, "log/", sample_name, ".log"),
                sep = " ")
  trim <- c(trim, temp)
}

for(i in 1:length(trim)){
  system(trim[i])
}

## 2. reference crawling and download
library(RCurl)
library(XML)

# make dir 
Ref_DIR <- paste0(WORK_DIR, "/2.reference/")
system(paste0("mkdir -p ", Ref_DIR))

# download
ref_url <- "https://ftp.ensembl.org/pub/release-106/fasta/"

ref <- getHTMLLinks(getURL(paste0(ref_url, species, "/", "dna/"), ftp.use.epsv = T, dirlistonly = T))

if(length((grep("dna.primary_assembly", ref)))==1){
  filenames <- ref[grep("dna.primary_assembly.fa.gz", ref)]
}else{
  filenames <- ref[grep("dna.toplevel.fa.gz", ref)]
}

system(paste("wget", paste0(ref_url, species, "/", "dna/", filenames), "-P", Ref_DIR, "-o", paste0(Ref_DIR, "reference_download_log.txt")))

## 2. indexing - bwa
# unzip
system(paste0("gzip -d ", Ref_DIR, "*.gz"))

# reference indexing 
system(paste("bwa index",
             "-p indexed_genome",
             paste0(Ref_DIR, "*.fa"),
             paste0("1> ", Ref_DIR, "bwa_index_log.txt")))

## 3. alignment - bwa
Align_DIR <- paste0(WORK_DIR, "/3.alignment/")
system(paste0("mkdir -p ", Align_DIR))

align <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("bwa mem", 
                paste0(Ref_DIR, "indexed_genome"),
                paste0(Trim_DIR, sample_name, "_1.fastq.gz"),
                paste0(Trim_DIR, sample_name, "_2.fastq.gz"),
                paste0("-t 48 > ", Align_DIR, sample_name, ".sam"),
                sep = " ")
  align <- c(align, temp)
}

for(i in 1:length(align)){
  system(align[i])
}


## 4. BAM - samtools
Bam_DIR <- paste0(WORK_DIR, "/4.BAM/")
system(paste0("mkdir -p ", Bam_DIR))

bam <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("samtools view",
                "-Sb -@ 24",
                paste0(Align_DIR, sample_name, ".sam"),
                paste0("> ", Bam_DIR, sample_name, ".bam"))
  bam <- c(bam, temp)
}

for(i in 1:length(bam)){
  system(bam[i])
}

## sorted BAM - samtools
sortedBam_DIR <- paste0(WORK_DIR, "/5.sorted_BAM/")
system(paste0("mkdir -p ", sortedBam_DIR))

sort_bam <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("samtools sort",
                "-@ 24 -o",
                paste0(sortedBam_DIR, sample_name, "_sorted.bam"),
                paste0(Bam_DIR, sample_name, ".bam"))
  sort_bam <- c(sort_bam, temp)
}

for(i in 1:length(sort_bam)){
  system(sort_bam[i])
}

## duplicated removal - picardMarkduplicate
Dup_DIR <- paste0(WORK_DIR, "/6.dupicated_remove/")
system(paste0("mkdir -p ", Dup_DIR))

picard <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("java -jar /program/picard.jar MarkDuplicates",
                paste0("I=", sortedBam_DIR, sample_name, "_sorted.bam"),
                paste0("O=", Dup_DIR, sample_name, "_marked_dup.bam"),
                paste0("M=", Dup_DIR, sample_name, "_marked_dup.txt"))
  picard <- c(picard, temp)
}

for(i in 1:length(picard)){
  system(picard[i])
}

## bigwig - deeptools
BW_DIR <- paste0(WORK_DIR, "/8.bigwig/")
system(paste0("mkdir -p ", BW_DIR))

#index
index <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("/program/samtools/bin/samtools index -@ 16",
                paste0(Dup_DIR, sample_name, "_marked_dup.bam"))
  index <- c(index, temp)
}
for(i in 1:length(index)){
  system(index[i])
}

bigwig <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("/program/anaconda3/bin/bamCoverage -b",
                paste0(Dup_DIR, sample_name, "_marked_dup.bam"), 
                paste0("--extendReads 200"),
                paste0("--normalizeUsing RPKM --binSize 1"),
                paste0("-o ", BW_DIR, sample_name, ".bw -p 24"))
  bigwig <- c(bigwig, temp)
}

for(i in 1:length(bigwig)){
  system(bigwig[i])
}

## bigwig2 (project 별) - deeptools
BW_DIR2 <- paste0(WORK_DIR, "/8.bigwig2/")
system(paste0("mkdir -p ", BW_DIR2))

# Bam merge
Bam_merge_DIR <- paste0(WORK_DIR, "/8.bigwig2/merged_BAM/")
system(paste0("mkdir -p ", Bam_merge_DIR))

L_4.5 <- c("SRR15584312", "SRR15584313", "SRR15584325")
R_4.5 <- c("SRR15584318", "SRR15584319", "SRR15584320")
L_5.5 <- c("SRR15584315", "SRR15584316", "SRR15584317")
R_5.5 <- c("SRR15584310", "SRR15584311", "SRR15584314")
L_7 <- c("SRR15584309", "SRR15584331", "SRR15584332")
R_7 <- c("SRR15584328", "SRR15584329", "SRR15584330")
L_10 <- c("SRR15584324", "SRR15584326", "SRR15584327")
R_10 <- c("SRR15584321", "SRR15584322", "SRR15584323")

pro <- cbind(L_4.5, R_4.5, L_5.5, R_5.5, L_7, R_7, L_10, R_10)

merge <- c()
for (i in 1:ncol(pro)) {
  merge[i] <- paste("samtools merge -@ 24",
                    paste0(Bam_merge_DIR, colnames(pro)[i], ".bam"),
                    paste0(Dup_DIR, pro[1,i], "_marked_dup.bam"),
                    paste0(Dup_DIR, pro[2,i], "_marked_dup.bam"),
                    paste0(Dup_DIR, pro[3,i], "_marked_dup.bam"))
}

for(i in 1:length(merge)){
  system(merge[i])
}

 
#index
index <- c()
temp <- c()
for (i in 1:ncol(pro)) {
  temp <- paste("/program/samtools/bin/samtools index -@ 16",
                paste0(Bam_merge_DIR, colnames(pro)[i], ".bam"))
  index <- c(index, temp)
}
for(i in 1:length(index)){
  system(index[i])
}

bigwig2 <- c()
temp <- c()
for (i in 1:ncol(pro)) {
  temp <- paste("/program/anaconda3/bin/bamCoverage -b",
                paste0(Bam_merge_DIR, colnames(pro)[i], ".bam"),
                paste0("--extendReads 200"),
                paste0("--normalizeUsing RPKM --binSize 1"),
                paste0("-o ", BW_DIR2, colnames(pro)[i], ".bw -p 24"))
  bigwig2 <- c(bigwig2, temp)
}

for(i in 1:length(bigwig2)){
  system(bigwig2[i])
}


## peak calling - macs2
Peak_DIR <- paste0(WORK_DIR, "/7.peak_calling/")
system(paste0("mkdir -p ", Peak_DIR))

peak <- c()
temp <- c()
for (sample_name in sample_names) {
  temp <- paste("macs2 callpeak",
                paste0("-t ", Dup_DIR, sample_name, "_marked_dup.bam"),
                paste0("--outdir ", Peak_DIR),
                paste0("-n ", sample_name),
                paste0("-g 9.6e8 -q 0.01 -B --SPMR --nomodel --shift -100 --extsize 200"))
  peak <- c(peak, temp)
}

for(i in 1:length(peak)){
  system(peak[i])
}
