# pichia
commands used in the assembly and annotation of Pichia genomes


##Building of directory structure

```bash
ProjDir=/home/groups/harrisonlab/project_files/alternaria
cd $ProjDir

Organism=A.alternata_ssp_tenuissima
Strain=1166
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=gaisen
Strain=650
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

```