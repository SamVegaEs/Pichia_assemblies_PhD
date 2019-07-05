# pichia
commands used in the assembly and annotation of Pichia genomes

##Demultiplexing reads using Albacore.
Data was basecalled again using Albacore 2.3.3 on the minion server:

#The programs were run in:

cd /home/groups/harrisonlab/project_files/Pichia

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
  OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz
```

##Building of directory structure

```bash
ProjDir=/home/groups/harrisonlab/project_files/Pichia
mkdir $ProjDir
cd $ProjDir

Organism=P.stipitis
Strain=589
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=P.stipitis
Strain=591
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=P.stipitis
Strain=594
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode03.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

```

#Doing it in in this way was not possible since the message of no space left in the device appeared. So we will try another way: First running the OutDir variable and then running the tar command. So:

```bash
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

  tar -cz -f $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed

```

#ASSEMBLY.

#Removal of the adapters: Splitting reads and trimming adapters using porechop. We only need the FASTA files for this step, so we can do it while the "tar" command is runing. To do so, we have to create the directory in our project folder and not the nanopore node.


```bash
for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=qc_dna/minion/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```

#Identify sequencing coverage

#For Minion data:

```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=16
OutDir=$(dirname $RawData)
mkdir -p $OutDir
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done


  for StrainDir in $(ls -d qc_dna/minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
Coverage.

589     186.88
591     145.17
594     158.09
```


#Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.8/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 16m $Strain $OutDir
done
```

#Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.8/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

#Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
# ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon2_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/*.fasta | grep 'round_10'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

#Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
# for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'FON_63' | grep 'racon_min_500bp_renamed'); do
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/P*/*/assembly/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

#Results.

```bash
short_summary_589_smartdenovo.dmo.lay.txt       531     13      306     478    1315
short_summary_589_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_589_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_6.txt 970     4       150     195    1                                                                                        315
short_summary_589_smartdenovo_racon_round_7.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_8.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_9.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo.dmo.lay.txt       501     6       321     493    1                                                                                        315
short_summary_591_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_591_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo.dmo.lay.txt       505     1       321     489    1                                                                                        315
short_summary_594_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_594_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_8.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_9.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315

```
#This results do not seem correct since some files have the value "0" for complete genes, that means that there are 0 genes identified. I checked the files after the run in racon, and they were empty, so I re-run the command. After running it I should check if there is information in them with ls -lh. After repeating it, the files had information on them, and BUSCO was run again. The results obtained after the second run are:


```bash
short_summary_589_smartdenovo.dmo.lay.txt         531     13      306     478    1315
short_summary_589_smartdenovo_racon_round_10.txt  975     5       158     182    1315
short_summary_589_smartdenovo_racon_round_1.txt   945     9       174     196    1315
short_summary_589_smartdenovo_racon_round_2.txt   949     13      169     197    1315
short_summary_589_smartdenovo_racon_round_3.txt   989     9       151     175    1315
short_summary_589_smartdenovo_racon_round_4.txt   977     10      142     196    1315
short_summary_589_smartdenovo_racon_round_5.txt   971     8       160     184    1315
short_summary_589_smartdenovo_racon_round_6.txt   970     4       150     195    1315
short_summary_589_smartdenovo_racon_round_7.txt   984     5       149     182    1315
short_summary_589_smartdenovo_racon_round_8.txt   970     5       145     200    1315
short_summary_589_smartdenovo_racon_round_9.txt   972     4       155     188    1315
short_summary_racon_min_500bp_renamed.txt         975     5       158     182    1315

short_summary_591_smartdenovo.dmo.lay.txt         501     6       321     493    1315
short_summary_591_smartdenovo_racon_round_10.txt  957     4       163     195    1315
short_summary_591_smartdenovo_racon_round_1.txt   944     3       159     212    1315
short_summary_591_smartdenovo_racon_round_2.txt   957     8       160     198    1315
short_summary_591_smartdenovo_racon_round_3.txt   960     7       170     185    1315
short_summary_591_smartdenovo_racon_round_4.txt   962     4       162     191    1315
short_summary_591_smartdenovo_racon_round_5.txt   973     4       161     181    1315
short_summary_591_smartdenovo_racon_round_6.txt   955     4       168     192    1315
short_summary_591_smartdenovo_racon_round_7.txt   962     4       160     193    1315
short_summary_591_smartdenovo_racon_round_8.txt   967     4       155     193    1315
short_summary_591_smartdenovo_racon_round_9.txt   973     4       149     193    1315
short_summary_racon_min_500bp_renamed.txt         957     4       163     195    1315


short_summary_594_smartdenovo.dmo.lay.txt         505     1       321     489    1315
short_summary_594_smartdenovo_racon_round_10.txt  950     1       170     195    1315
short_summary_594_smartdenovo_racon_round_1.txt   936     0       182     197    1315
short_summary_594_smartdenovo_racon_round_2.txt   949     1       166     200    1315
short_summary_594_smartdenovo_racon_round_3.txt   948     1       175     192    1315
short_summary_594_smartdenovo_racon_round_4.txt   951     0       178     186    1315
short_summary_594_smartdenovo_racon_round_5.txt   961     1       168     186    1315
short_summary_594_smartdenovo_racon_round_6.txt   949     2       173     193    1315
short_summary_594_smartdenovo_racon_round_7.txt   948     1       172     195    1315
short_summary_594_smartdenovo_racon_round_8.txt   973     0       165     177    1315
short_summary_594_smartdenovo_racon_round_9.txt   946     0       178     191    1315
short_summary_racon_min_500bp_renamed.txt         950     1       170     195    1315
```

#After problems running the tar command we ran it again from inside the nanopore node.

```bash
screen -a

ssh nanopore@nanopore

cd /home/nanopore/Pichia_31_01_2019/31-01-19/
mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

tar -cz -f $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed

```
#To check if the tar command is still running we have to resume the screen in which it was running.


#ASSEMBLY CORRECTION USING NANOPOLISH.

#Fast5 files are very large and need to be stored as gzipped tarballs. These needed temporarily unpacking but must be deleted after nanpolish has finished running.

#Raw reads were moved onto the cluster scratch space for this step and unpacked:

```bash
screen -a

for Tar in $(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Pichia_albacore_v2.3.3_demultiplexed.tar.gz); do
  ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
  mkdir -p $ScratchDir
  tar -zxvf $Tar -C $ScratchDir
done
```
```bash


# When you try to run this part in a screen session it asks for an update in certain programs required to run nanopolish, so do not run it in screen.

for Assembly in $(ls assembly/SMARTdenovo/*/589/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

for Assembly in $(ls assembly/SMARTdenovo/*/591/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

for Assembly in $(ls assembly/SMARTdenovo/*/594/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

```

# Split the assembly into 50Kb fragments an submit each to the cluster for nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
RawReads=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish_range.txt); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done    
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

# When the code for splitting was running we had to stop it. To take it from where it stopped and not from the beginning:
```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
RawReads=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log2.txt
ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt"); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done    
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log2.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
InDir=$(dirname $Assembly)
python $NanoPolishDir/nanopolish_merge.py $InDir/nanopolish/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```
#Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

#For Pilon assembly correction we need the Illumina data, and we don't have them yet in our working directory. We need to copy them there.

```bash
RawDatDir=/data/seq_data/miseq/2018/RAW/181026_M04465_0088_000000000-C5M4V/Data/Intensities/BaseCalls

  ProjectDir=/home/groups/harrisonlab/project_files/Pichia

# Isolate 589

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/589

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-589_S2_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-589_S2_L001_R2_001.fastq.gz  .

  cd $ProjectDir

# Isolate 591

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/591

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-591_S3_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-591_S3_L001_R2_001.fastq.gz  .

  cd $ProjectDir

# Isolate 594

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/594

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-594_S4_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-594_S4_L001_R2_001.fastq.gz  .

  cd $ProjectDir

```

#QC of MiSeq Data.

# programs: fastqc fastq-mcf kmc
#Data quality was visualised using fastqc:

```bash
for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```
```bash
for StrainPath in $(ls -d raw_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/F/*.fastq*)
    ReadsR=$(ls $StrainPath/R/*.fastq*)
    echo $ReadsF
    echo $ReadsR
    qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```
# Pilon assembly correction
# Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/Pichia/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=$(dirname $Assembly)/../pilon
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

# Contigs were renamed.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_10'); do
echo $Assembly
echo "" > tmp.txt
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

# Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  echo "$Organism - $Strain"
  # OutDir=$(dirname $Assembly)
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  done
  printf "\n"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/SMARTdenovo/*/*/pilon/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```
# Summary results:

```bash
short_summary_pilon_10.txt                      1290    11      5       20      1315
short_summary_pilon_1.txt                       1290    11      5       20      1315
short_summary_pilon_2.txt                       1290    11      5       20      1315
short_summary_pilon_3.txt                       1290    11      5       20      1315
short_summary_pilon_4.txt                       1290    11      5       20      1315
short_summary_pilon_5.txt                       1290    11      5       20      1315
short_summary_pilon_6.txt                       1290    11      5       20      1315
short_summary_pilon_7.txt                       1290    11      5       20      1315
short_summary_pilon_8.txt                       1290    11      5       20      1315
short_summary_pilon_9.txt                       1290    11      5       20      1315
short_summary_pilon_min_500bp_renamed.txt       1290    11      5       20      1315

short_summary_pilon_10.txt                      1284    3       7       24      1315
short_summary_pilon_1.txt                       1284    3       6       25      1315
short_summary_pilon_2.txt                       1285    3       6       24      1315
short_summary_pilon_3.txt                       1284    3       7       24      1315
short_summary_pilon_4.txt                       1284    3       7       24      1315
short_summary_pilon_5.txt                       1284    3       7       24      1315
short_summary_pilon_6.txt                       1284    3       7       24      1315
short_summary_pilon_7.txt                       1284    3       7       24      1315
short_summary_pilon_8.txt                       1284    3       7       24      1315
short_summary_pilon_9.txt                       1284    3       7       24      1315
short_summary_pilon_min_500bp_renamed.txt       1284    3       7       24      1315

short_summary_pilon_10.txt                      1288    1       5       22      1315
short_summary_pilon_1.txt                       1288    1       5       22      1315
short_summary_pilon_2.txt                       1288    1       5       22      1315
short_summary_pilon_3.txt                       1288    1       5       22      1315
short_summary_pilon_4.txt                       1288    1       5       22      1315
short_summary_pilon_5.txt                       1288    1       5       22      1315
short_summary_pilon_6.txt                       1288    1       5       22      1315
short_summary_pilon_7.txt                       1288    1       5       22      1315
short_summary_pilon_8.txt                       1288    1       5       22      1315
short_summary_pilon_9.txt                       1288    1       5       22      1315
short_summary_pilon_min_500bp_renamed.txt       1288    1       5       22      1315
```

```bash

# To get the information about the total length of the genome, and the number of contigs we have we use this command:  

```bash
for File in $(ls assembly/SMARTdenovo/P.stipitis/*/racon2_10/report.txt); do
echo $File
cat $File
done

589

All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        11
# contigs (>= 1000 bp)     11
Total length (>= 0 bp)     15687345
Total length (>= 1000 bp)  15687345
# contigs                  11
Largest contig             2692878
Total length               15687345
GC (%)                     41.02
N50                        1881524
N75                        1677956
L50                        4
L75                        6
# N's per 100 kbp          0.00

591

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        11
# contigs (>= 1000 bp)     11
Total length (>= 0 bp)     15639864
Total length (>= 1000 bp)  15639864
# contigs                  11
Largest contig             3983890
Total length               15639864
GC (%)                     41.11
N50                        1904357
N75                        1319906
L50                        3
L75                        6
# N's per 100 kbp          0.00

594

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        9
# contigs (>= 1000 bp)     9
Total length (>= 0 bp)     15358923
Total length (>= 1000 bp)  15358923
# contigs                  9
Largest contig             3438607
Total length               15358923
GC (%)                     41.02
N50                        1898369
N75                        1767446
L50                        3
L75                        5
# N's per 100 kbp          0.00



```

# HYBRID ASSEMBLY.

# Spades Assembly

```bash
for TrimReads in $(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/*/*.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/Pichia/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=assembly/spades_minion/$Organism/"$Strain"
echo $TrimF1_Read
echo $TrimR1_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
qsub $ProgDir/sub_spades_minion.sh $TrimReads $TrimF1_Read $TrimR1_Read $OutDir
done
```

# Contigs shorter than 500bp were removed from the assembly.

```bash
for Contigs in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```

# Quast and BUSCO

```bash
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
#To do yet:

```bash
 for File in $(ls assembly/spades_minion/*/*/*/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

#The folder polished is not present, so I have to create it before going on with the Merging of both assemblies.


# Repeat the run of Quast and busco were run to assess the effects of pilon on assembly quality, but in this case using the data from  saccharomycetales_odb9, since it's Pichia's order.


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  echo "$Organism - $Strain"
  # OutDir=$(dirname $Assembly)
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  done
  printf "\n"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/saccharomycetales_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/SMARTdenovo/*/*/pilon/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```


# Summary results:

```bash
short_summary_pilon_10.txt                      1683    13      11      17      1711
short_summary_pilon_1.txt                       1682    13      12      17      1711
short_summary_pilon_2.txt                       1683    13      11      17      1711
short_summary_pilon_3.txt                       1683    13      11      17      1711
short_summary_pilon_4.txt                       1683    13      11      17      1711
short_summary_pilon_5.txt                       1683    13      11      17      1711
short_summary_pilon_6.txt                       1683    13      11      17      1711
short_summary_pilon_7.txt                       1683    13      11      17      1711
short_summary_pilon_8.txt                       1683    13      11      17      1711
short_summary_pilon_9.txt                       1683    13      11      17      1711
short_summary_pilon_min_500bp_renamed.txt       1683    13      11      17      1711

short_summary_pilon_10.txt                      1678    9       14      19      1711
short_summary_pilon_1.txt                       1676    9       15      20      1711
short_summary_pilon_2.txt                       1678    9       13      20      1711
short_summary_pilon_3.txt                       1678    9       14      19      1711
short_summary_pilon_4.txt                       1678    9       14      19      1711
short_summary_pilon_5.txt                       1678    9       14      19      1711
short_summary_pilon_6.txt                       1678    9       14      19      1711
short_summary_pilon_7.txt                       1678    9       14      19      1711
short_summary_pilon_8.txt                       1678    9       14      19      1711
short_summary_pilon_9.txt                       1678    9       14      19      1711
short_summary_pilon_min_500bp_renamed.txt       1678    9       14      19      1711

short_summary_pilon_10.txt                      1685    5       12      14      1711
short_summary_pilon_1.txt                       1684    5       13      14      1711
short_summary_pilon_2.txt                       1685    5       12      14      1711
short_summary_pilon_3.txt                       1685    5       12      14      1711
short_summary_pilon_4.txt                       1685    5       12      14      1711
short_summary_pilon_5.txt                       1685    5       12      14      1711
short_summary_pilon_6.txt                       1685    5       12      14      1711
short_summary_pilon_7.txt                       1685    5       12      14      1711
short_summary_pilon_8.txt                       1685    5       12      14      1711
short_summary_pilon_9.txt                       1685    5       12      14      1711
short_summary_pilon_min_500bp_renamed.txt       1685    5       12      14      1711
```

#Remove mit. DNA

#Using an exclusion database with deconseq:

```bash
 for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "Calb_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

#Results were summarised using the commands:

```bash
for Exclude_db in "Calb_mtDNA"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```
```bash
Calb_mtDNA
589     11      0
591     11      0
594     8       1
```
#Quast was run on the removed mtDNA:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_Calb_mtDNA/*_cont.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

#Repeat Masking.
#Repeat masking was done for the Nanopore assembly polished with Illumina data (SMARTdenovo).

```bash
 # for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
  for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_Calb_mtDNA/contigs_min_500bp_filtered_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```

# The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
Output of previous run:

repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
520465
repeat_masked/P.stipitis/591/filtered_contigs/591_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
591492
repeat_masked/P.stipitis/594/filtered_contigs/594_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
539679
```
```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3); do
Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
DDE_1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'DDE_1' | sed "s/^\s*//g" | cut -f1 -d ' ')
Gypsy=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'gypsy' | sed "s/^\s*//g" | cut -f1 -d ' ')
HAT=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'hAT' | sed "s/^\s*//g" | cut -f1 -d ' ')
TY1_Copia=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'TY1_Copia' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep -w 'mariner' | sed "s/^\s*//g" | cut -f1 -d ' ')
Cacta=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'cacta' | sed "s/^\s*//g" | cut -f1 -d ' ')
LINE=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'LINE' | sed "s/^\s*//g" | cut -f1 -d ' ')
MuDR_A_B=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'MuDR_A_B' | sed "s/^\s*//g" | cut -f1 -d ' ')
HelitronORF=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'helitronORF' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner_ant1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'mariner_ant1' | sed "s/^\s*//g" | cut -f1 -d ' ')
ISC1316=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'ISC1316' | sed "s/^\s*//g" | cut -f1 -d ' ')
Crypton=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'Crypton' | sed "s/^\s*//g" | cut -f1 -d ' ')
printf "$Organism\t$Strain\t$DDE_1\t$Gypsy\t$HAT\t$TY1_Copia\t$Mariner\t$Cacta\t$LINE\t$MuDR_A_B\t$HelitronORF\t$Mariner_ant1\t$ISC1316\t$Crypton\n"
done
```
```bash
Output of previous run:

P.stipitis      589             13      1       99      1               13     14       8
P.stipitis      591             12      1       112     1       1       14     14       8
P.stipitis      594             13      1       108     1               14     14       8
```

#Quast and BUSCO

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/saccharomycetales_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/run_*_contigs_unmasked/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```
```bash
Output of Busco

# P.stipitis      589     1683    11      17      1711
# P.stipitis      591     1678    14      19      1711
# P.stipitis      594     1685    12      14      1711
```

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/report.tsv); do
  Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
  Contigs=$(cat $File | grep "contigs (>= 0 bp)" | cut -f2)
  Length=$(cat $File | grep "Total length (>= 0 bp)" | cut -f2)
  Largest=$(cat $File | grep "Largest contig" | cut -f2)
  N50=$(cat $File | grep "N50" | cut -f2)
  echo -e "$Organism\t$Strain\t$Contigs\t$Length\t$Largest\t$N50"
  done
```
```bash
Output:

P.stipitis      589     11      15681772        2691698 1880751
P.stipitis      591     11      15602566        3982952 1903767
P.stipitis      594     8       15273029        3437457 1897630
```

#Promer alignment of Assemblies.

#Alignment of our assemblies against the reference S. stipitis genome (AB580, Y-11545_v2)

```bash
Reference=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Y-11454_v2
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```



#Alignment of all our assemblies against the assembly done for 589

```bash
Reference=$(ls repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_589
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

#Alignment of all our assemblies against the assembly done for 591

```bash
Reference=$(ls repeat_masked/P.stipitis/591/*/591_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_591
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

#Alignment of all our assemblies against the assembly done for 594

```bash
Reference=$(ls repeat_masked/P.stipitis/594/*/594_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_594
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```


#Alignment of S. stipitis raw reads vs reference genome.

#####################################
## All against reference genome.  ###
#####################################

#Alignment of reads from a single run:

```bash
for Reference in $(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked); do
for StrainPath in $(ls -d qc_dna/paired/P.stipitis/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```
