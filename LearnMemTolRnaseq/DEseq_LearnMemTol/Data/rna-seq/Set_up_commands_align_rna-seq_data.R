#Set up all the commands we need to do the alignment
#Set up a job array

#one command
#hisat2 --dta -x /group/kinglab/Patricka/base_pop/indexes/bdgp6_tran/genome_tran -1 "/storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_raw/RAPiD-Genomics_HJYM3BBXX_MIZ_117501_P01_WB12_i5-509_i7-75_S577_L002_trimmed_R1_001.fastq.gz" -2 "/storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_raw/RAPiD-Genomics_HJYM3BBXX_MIZ_117501_P01_WB12_i5-509_i7-75_S577_L002_trimmed_R2_001.fastq.gz" -S processed_align_test_2.sam


cat('', file="../rna-seq_data_raw/Align_cmds_LearnMemRNA.txt")

ll<- list.files("/storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_raw/Align_cmds_LearnMemRNA.txt")


base.cmd <- "hisat2 --dta -x /group/kinglab/Patricka/base_pop/indexes/bdgp6_tran/genome_tran -1 /storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_trimmed/"

ll <- gsub("trimmed_R1_001.fastq.gz", "", ll, fixed=TRUE)
ll <- gsub("trimmed_R2_001.fastq.gz", "", ll, fixed=TRUE)
ll <- unique(ll)

for(ff in ll) 
  {
  cmd.cut <- paste(base.cmd,ff,"trimmed_R1_001.fastq.gz -2 ",ff, "trimmed_R2_001.fastq.gz",
                   " -S ../rna-seq_data_aligned/",ff,"aligned_rna-seq.sam",sep="")
  cat(cmd.cut,"\n", file="../../RawData/Align_cmds_LearnMemRNA.txt",append=TRUE)
  
}

#for(i in 1:114)
#{
#  cat(paste("temp", i,".txt", sep=""), "\n", file = "Trim_cmds_temp.txt", append=TRUE)
#  cat(paste("error", i,".txt", sep=""), "\n", file = "Trim_cmds_error.txt", append=TRUE)
#}




