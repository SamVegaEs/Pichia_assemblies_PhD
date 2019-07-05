# pichia
commands used in the assembly and annotation of Pichia genomes

##Demultiplexing reads using Albacore.
Data was basecalled again using Albacore 2.3.3 on the minion server:


```bash
screen -a  

#Sceen -a opens a new session.

ssh nanopore@nanopore

~/.local/bin/read_fast5_basecaller.py --version

mkdir Pichia_31_01_2019
cd Pichia_31_01_2019

# Oxford nanopore 07/03/17
Organism=P.stipitis
Date=31-01-19
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
RawDatDir=/data/seq_data/minion/2018/Pichia/GA20000/reads
OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3
mkdir -p $OutDir

mkdir -p ~/Pichia_31_01_2019/$Date
cd ~/Pichia_31_01_2019/$Date
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path Pichia_albacore_v2.3.3_demultiplexed \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000 \
  --barcoding

#Run the commmand ls to check that the folders actually exist.

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01/*.fastq |wc -l

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02/*.fastq |wc -l

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03/*.fastq |wc -l

#Run each cat command individually and check that the output exists.

  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode01.fastq.gz
  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode02.fastq.gz
  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode03.fastq.gz

#Set the space where the output data will be.

  OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3
  mkdir -p $OutDir
  chmod +w $OutDir
  cp Pichia_albacore_v2.3.3_barcode0*.fastq.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

  #Doing in in this way was not possible since the message of no space left in the device appeared. So we will try another way: First running the OutDir variable and then running the tar command. So:

  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

  tar -cz -f $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed

```

##Building of directory structure

```bash
ProjDir=/home/groups/harrisonlab/project_files/Pichia
cd $ProjDir

Organism=P.stipitis
Strain=589
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=P.stipitis
Strain=591
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

rganism=P.stipitis
Strain=594
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

```


#ASSEMBLY.

#Removal of the adapters: Splitting reads and trimming adapters using porechop. We only need the FASTA files for this step, so we can do it while the "tar" command is runing.

for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=qc_dna/minion/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
