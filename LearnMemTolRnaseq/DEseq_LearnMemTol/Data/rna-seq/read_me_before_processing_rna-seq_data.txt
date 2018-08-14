###Steps to processing RNA-seq data


###Preprocessing raw RNA-seq data

1.to stack raw data for trimming off adaptors 

#Use R to convert raw data to collated data. Saved on Navalis

#Set up all the commands we need to do the trimming
#Set up a job array

#one command
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --output _trimmed_R1.fastq.gz --paired-output _trimmed_R2.fastq.gz  _R1.fastq.gz _R2.fastq.gz


base.cmd <- "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b A{100} -b T{100} --trim-n -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -B A{100} -B T{100} --output"

load(file="~/DSPR/Scripts/rna-seq_processing_raw_data/filelist.rda")

filelist <- gsub("_R1_001.fastq.gz", "", filelist, fixed=TRUE)
filelist <- gsub("_R2_001.fastq.gz", "", filelist, fixed=TRUE)
filelist <- unique(filelist)

for(ff in filelist) 
  {
  cmd.cut <- paste(base.cmd," ",ff, "_trimmed_R1_001.fastq.gz --paired-output ",
                   ff, "_trimmed_R2_001.fastq.gz ", ff, "_R1_001.fastq.gz ",ff, "_R2_001.fastq.gz", sep="")
  cat(cmd.cut,"\n", file="Trim_cmds_LearnMemRNA.txt",append=TRUE)
  
}

#for(i in 1:114)
#{
#  cat(paste("temp", i,".txt", sep=""), "\n", file = "Trim_cmds_temp.txt", append=TRUE)
#  cat(paste("error", i,".txt", sep=""), "\n", file = "Trim_cmds_error.txt", append=TRUE)
#}


2. Import data to Lewis using rsync


#To run a sbatch file


1. Set up nano file. Saved as nano sbatch_trimRNA.sh on Lewis


 
2.Make sure you are in the folder where the raw data is located

3.Script to run sbatach array. sbatch --array=3-114 ../sbatch_files/sbatch_triRNA.sh

 