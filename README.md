# pichia
commands used in the assembly and annotation of Pichia genomes

##Demultiplexing reads using Albacore. 
Data was basecalled again using Albacore 2.3.3 on the minion server:


```bash
screen -a  

#Sceen -a opens a new session. 

ssh nanopore@nanopore

~/.local/bin/read_fast5_basecaller.py --version

mkdir Pichia_14_12_2018
cd Pichia_14_12_2018

# Oxford nanopore 07/03/17
Organism=P.stipitis
Date=14-12-18
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
RawDatDir=/data/seq_data/minion/2018/Pichia/GA20000/reads
OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3
mkdir -p $OutDir

mkdir -p ~/Pichia_14_12_2018/$Date
cd ~/Pichia_14_12_2018/$Date
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