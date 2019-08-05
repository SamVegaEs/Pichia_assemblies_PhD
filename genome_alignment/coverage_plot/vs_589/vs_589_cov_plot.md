  
```bash
 OutDir=vs_589_cov_plot
  mkdir -p $OutDir

  Ref_genome=$(ls repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_hardmasked_repeatmasker_TPSI_appended.fa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/ref_genome.txt
  ```

  # Calculate coverage over 10kb windows for unmasked files. Parse read depth files: 

```bash
  OutDir=vs_589_cov_plot
for File in $(ls alignment/genome_alignment/vs_589_unmasked/*/*_vs_589_unmasked_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f4 -d '/')
  echo $Strain 
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/contig_/589_contig_/g' > $OutDir/${Strain}_vs_ref_hardmasked_scatterplot.tsv
done

```bash
ProgDir=/home/vegasa/git_repos/scripts/pichia/genome_alignment/coverage_plot/vs_589
circos -conf $ProgDir/vs_589_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/vs_589_cov_unmasked_circos.png
mv $OutDir/circos.svg $OutDir/vs_589_cov_unmasked_circos.svg
ls $PWD/$OutDir/vs_589*_circos.png
```

#Reference unmasked. 

```bash
 OutDir=vs_589_cov_plot
  mkdir -p $OutDir

  Ref_genome=$(ls repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/ref_genome_softmasked.txt
  ```

  # Calculate coverage over 10kb windows for unmasked files. Parse read depth files: 

```bash
  OutDir=vs_589_cov_plot
for File in $(ls alignment/genome_alignment/vs_589_unmasked/*/*_vs_589_unmasked_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f4 -d '/')
  echo $Strain 
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/contig_/589_contig_/g' > $OutDir/${Strain}_vs_ref_softmasked_scatterplot.tsv
done

```bash
ProgDir=/home/vegasa/git_repos/scripts/pichia/genome_alignment/coverage_plot/vs_589
circos -conf $ProgDir/vs_589_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/vs_589_cov_softmasked_circos.png
mv $OutDir/circos.svg $OutDir/vs_589_cov_softmasked_circos.svg
ls $PWD/$OutDir/vs_589*_circos.png
```


# Differences observed?

```bash

#1 Strong coverage in a region of strain 591 compared to both 589 and 594.

less alignment/genome_alignment/vs_589_unmasked/P.stipitis_591$ less P.stipitis_591_vs_589_unmasked_depth_10kb.tsv

contig_3        1071000 115     P.stipitis_591
contig_3        1072000 103     P.stipitis_591
contig_3        1073000 161     P.stipitis_591
contig_3        1074000 239     P.stipitis_591
contig_3        1075000 265     P.stipitis_591
contig_3        1076000 410     P.stipitis_591
contig_3        1077000 323     P.stipitis_591
contig_3        1078000 392     P.stipitis_591
contig_3        1079000 387     P.stipitis_591
contig_3        1080000 372     P.stipitis_591
contig_3        1081000 359     P.stipitis_591
contig_3        1082000 362     P.stipitis_591
contig_3        1083000 391     P.stipitis_591
contig_3        1084000 409     P.stipitis_591
contig_3        1085000 417     P.stipitis_591
contig_3        1086000 433     P.stipitis_591
contig_3        1087000 347     P.stipitis_591
contig_3        1088000 385     P.stipitis_591
contig_3        1089000 374     P.stipitis_591
contig_3        1090000 397     P.stipitis_591
contig_3        1091000 362     P.stipitis_591
contig_3        1092000 387     P.stipitis_591
contig_3        1093000 354     P.stipitis_591
contig_3        1094000 362     P.stipitis_591
contig_3        1095000 383     P.stipitis_591
contig_3        1096000 345     P.stipitis_591
contig_3        1097000 382     P.stipitis_591
contig_3        1098000 366     P.stipitis_591
contig_3        1099000 352     P.stipitis_591
contig_3        1100000 358     P.stipitis_591
contig_3        1101000 329     P.stipitis_591
contig_3        1102000 351     P.stipitis_591
contig_3        1103000 357     P.stipitis_591
contig_3        1104000 350     P.stipitis_591
contig_3        1105000 370     P.stipitis_591
contig_3        1106000 387     P.stipitis_591
contig_3        1107000 391     P.stipitis_591
contig_3        1108000 356     P.stipitis_591
contig_3        1109000 400     P.stipitis_591
contig_3        1110000 402     P.stipitis_591
contig_3        1111000 360     P.stipitis_591
contig_3        1112000 363     P.stipitis_591
contig_3        1113000 376     P.stipitis_591
contig_3        1114000 331     P.stipitis_591
contig_3        1115000 322     P.stipitis_591
contig_3        1116000 306     P.stipitis_591
contig_3        1117000 327     P.stipitis_591
contig_3        1118000 407     P.stipitis_591
contig_3        1119000 416     P.stipitis_591
contig_3        1120000 400     P.stipitis_591
contig_3        1121000 367     P.stipitis_591
contig_3        1122000 350     P.stipitis_591
contig_3        1123000 378     P.stipitis_591
contig_3        1124000 362     P.stipitis_591
contig_3        1125000 389     P.stipitis_591
contig_3        1126000 383     P.stipitis_591
contig_3        1127000 380     P.stipitis_591
contig_3        1128000 350     P.stipitis_591
contig_3        1129000 321     P.stipitis_591
contig_3        1130000 404     P.stipitis_591
contig_3        1131000 347     P.stipitis_591
contig_3        1132000 379     P.stipitis_591
contig_3        1133000 353     P.stipitis_591
contig_3        1134000 371     P.stipitis_591
contig_3        1135000 378     P.stipitis_591
contig_3        1136000 359     P.stipitis_591
contig_3        1137000 358     P.stipitis_591
contig_3        1138000 354     P.stipitis_591
contig_3        1139000 354     P.stipitis_591
contig_3        1140000 333     P.stipitis_591
contig_3        1141000 356     P.stipitis_591
contig_3        1142000 333     P.stipitis_591
contig_3        1143000 378     P.stipitis_591
contig_3        1144000 326     P.stipitis_591
contig_3        1145000 400     P.stipitis_591
contig_3        1146000 349     P.stipitis_591
contig_3        1147000 327     P.stipitis_591
contig_3        1148000 331     P.stipitis_591
contig_3        1149000 372     P.stipitis_591
contig_3        1150000 356     P.stipitis_591
contig_3        1151000 356     P.stipitis_591
contig_3        1152000 366     P.stipitis_591
contig_3        1153000 138     P.stipitis_591

#2 Stronger coverage area in both strains 591 and 594 compared to 589.

File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_594/P.stipitis_594_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_2' | less -S

contig_2        1219000 86      P.stipitis_594
contig_2        1220000 143     P.stipitis_594
contig_2        1221000 112     P.stipitis_594
contig_2        1222000 110     P.stipitis_594
contig_2        1223000 219     P.stipitis_594
contig_2        1224000 117     P.stipitis_594
contig_2        1225000 89      P.stipitis_594
contig_2        1226000 91      P.stipitis_594

File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_591/P.stipitis_591_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_2' | less -S

contig_2        1205000 48      P.stipitis_591
contig_2        1206000 66      P.stipitis_591
contig_2        1207000 57      P.stipitis_591
contig_2        1208000 88      P.stipitis_591
contig_2        1209000 98      P.stipitis_591
contig_2        1210000 110     P.stipitis_591
contig_2        1211000 88      P.stipitis_591
contig_2        1212000 122     P.stipitis_591
contig_2        1213000 132     P.stipitis_591
contig_2        1214000 91      P.stipitis_591
contig_2        1215000 69      P.stipitis_591
contig_2        1216000 73      P.stipitis_591
contig_2        1217000 89      P.stipitis_591
contig_2        1218000 82      P.stipitis_591
contig_2        1219000 42      P.stipitis_591


#3 Region with high coverage in 589 and 594 is missing in 591 in contig 1 (between the positions at 0.4 Mbp and 0.5 Mbp).

File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_589/P.stipitis_589_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_1' | less -S

contig_1        399000  123     P.stipitis_589
contig_1        400000  101     P.stipitis_589
contig_1        401000  101     P.stipitis_589
contig_1        402000  103     P.stipitis_589
contig_1        403000  100     P.stipitis_589
contig_1        404000  93      P.stipitis_589
contig_1        405000  85      P.stipitis_589
contig_1        406000  90      P.stipitis_589
contig_1        407000  96      P.stipitis_589
contig_1        408000  93      P.stipitis_589
contig_1        409000  108     P.stipitis_589
contig_1        410000  100     P.stipitis_589
contig_1        411000  103     P.stipitis_589
contig_1        412000  94      P.stipitis_589
contig_1        413000  95      P.stipitis_589
contig_1        414000  87      P.stipitis_589
contig_1        415000  103     P.stipitis_589
contig_1        416000  96      P.stipitis_589
contig_1        417000  88      P.stipitis_589
contig_1        418000  81      P.stipitis_589
contig_1        419000  93      P.stipitis_589
contig_1        420000  85      P.stipitis_589
contig_1        421000  77      P.stipitis_589
contig_1        422000  83      P.stipitis_589
contig_1        423000  92      P.stipitis_589
contig_1        424000  96      P.stipitis_589
contig_1        425000  93      P.stipitis_589
contig_1        426000  95      P.stipitis_589
contig_1        427000  100     P.stipitis_589
contig_1        428000  93      P.stipitis_589
contig_1        429000  99      P.stipitis_589
contig_1        430000  89      P.stipitis_589
contig_1        431000  90      P.stipitis_589
contig_1        432000  68      P.stipitis_589
contig_1        433000  90      P.stipitis_589
contig_1        434000  91      P.stipitis_589
contig_1        435000  90      P.stipitis_589
contig_1        436000  82      P.stipitis_589
contig_1        437000  91      P.stipitis_589
contig_1        438000  90      P.stipitis_589
contig_1        439000  96      P.stipitis_589
contig_1        440000  92      P.stipitis_589
contig_1        441000  99      P.stipitis_589
contig_1        442000  87      P.stipitis_589
contig_1        443000  71      P.stipitis_589
contig_1        444000  78      P.stipitis_589
contig_1        445000  75      P.stipitis_589
contig_1        446000  87      P.stipitis_589
contig_1        447000  85      P.stipitis_589
contig_1        448000  92      P.stipitis_589
contig_1        449000  89      P.stipitis_589
contig_1        450000  96      P.stipitis_589
contig_1        451000  89      P.stipitis_589
contig_1        452000  96      P.stipitis_589
contig_1        453000  97      P.stipitis_589
contig_1        454000  102     P.stipitis_589
contig_1        455000  90      P.stipitis_589
contig_1        456000  94      P.stipitis_589
contig_1        457000  92      P.stipitis_589
contig_1        458000  70      P.stipitis_589
contig_1        459000  85      P.stipitis_589
contig_1        460000  92      P.stipitis_589
contig_1        461000  99      P.stipitis_589
contig_1        462000  92      P.stipitis_589
contig_1        463000  86      P.stipitis_589
contig_1        464000  82      P.stipitis_589
contig_1        465000  82      P.stipitis_589
contig_1        466000  95      P.stipitis_589
contig_1        467000  71      P.stipitis_589
contig_1        468000  98      P.stipitis_589
contig_1        469000  84      P.stipitis_589
contig_1        470000  95      P.stipitis_589
contig_1        471000  95      P.stipitis_589
contig_1        472000  51      P.stipitis_589
contig_1        473000  95      P.stipitis_589
contig_1        474000  173     P.stipitis_589
contig_1        475000  81      P.stipitis_589
contig_1        476000  110     P.stipitis_589
contig_1        477000  97      P.stipitis_589
contig_1        478000  91      P.stipitis_589
contig_1        479000  104     P.stipitis_589
contig_1        480000  109     P.stipitis_589
contig_1        481000  76      P.stipitis_589
contig_1        482000  90      P.stipitis_589
contig_1        483000  105     P.stipitis_589
contig_1        484000  97      P.stipitis_589
contig_1        485000  95      P.stipitis_589
contig_1        486000  93      P.stipitis_589
contig_1        487000  99      P.stipitis_589
contig_1        488000  96      P.stipitis_589
contig_1        489000  108     P.stipitis_589
contig_1        490000  96      P.stipitis_589
contig_1        491000  94      P.stipitis_589
contig_1        492000  97      P.stipitis_589
contig_1        493000  87      P.stipitis_589
contig_1        494000  89      P.stipitis_589
contig_1        495000  83      P.stipitis_589
contig_1        496000  91      P.stipitis_589
contig_1        497000  72      P.stipitis_589
contig_1        498000  81      P.stipitis_589
contig_1        499000  96      P.stipitis_589
contig_1        500000  87      P.stipitis_589
contig_1        501000  80      P.stipitis_589
contig_1        502000  99      P.stipitis_589
contig_1        503000  95      P.stipitis_589

#4 Small peak in contig 4 present in 591 compared to both 589 and 594.

File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_591/P.stipitis_591_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_4' | less -S

contig_4        1409000 53      P.stipitis_591
contig_4        1410000 59      P.stipitis_591
contig_4        1411000 63      P.stipitis_591
contig_4        1412000 73      P.stipitis_591
contig_4        1413000 59      P.stipitis_591
contig_4        1414000 42      P.stipitis_591
contig_4        1415000 63      P.stipitis_591
contig_4        1416000 60      P.stipitis_591
contig_4        1417000 93      P.stipitis_591
contig_4        1418000 102     P.stipitis_591
contig_4        1419000 128     P.stipitis_591
contig_4        1420000 96      P.stipitis_591
contig_4        1421000 57      P.stipitis_591
contig_4        1422000 64      P.stipitis_591
contig_4        1423000 50      P.stipitis_591
contig_4        1424000 78      P.stipitis_591
contig_4        1425000 64      P.stipitis_591
contig_4        1426000 63      P.stipitis_591
contig_4        1427000 74      P.stipitis_591


#5 Small peak present in 589 and 591, but not in 594 in contig 5 (Between the positions at 1.0 Mbp and 1.1 Mbp).
#6 Small peak present in 589 and 594, but not in 591 in contig 5 (Between the positions at 1.1 Mbp and 1.3 Mbp).


File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_589/P.stipitis_589_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_5' | less -S

contig_5        1021000 89      P.stipitis_589
contig_5        1022000 100     P.stipitis_589
contig_5        1023000 103     P.stipitis_589
contig_5        1024000 106     P.stipitis_589
contig_5        1025000 119     P.stipitis_589
contig_5        1026000 96      P.stipitis_589
contig_5        1027000 114     P.stipitis_589
contig_5        1028000 118     P.stipitis_589
contig_5        1029000 81      P.stipitis_589
contig_5        1030000 96      P.stipitis_589
contig_5        1031000 110     P.stipitis_589
contig_5        1032000 127     P.stipitis_589
contig_5        1033000 104     P.stipitis_589
contig_5        1034000 156     P.stipitis_589
contig_5        1035000 124     P.stipitis_589
contig_5        1036000 90      P.stipitis_589
contig_5        1037000 87      P.stipitis_589
contig_5        1038000 108     P.stipitis_589
contig_5        1039000 119     P.stipitis_589
contig_5        1040000 103     P.stipitis_589
contig_5        1041000 106     P.stipitis_589
contig_5        1042000 102     P.stipitis_589
contig_5        1043000 109     P.stipitis_589

contig_5        1258000 102     P.stipitis_589
contig_5        1259000 114     P.stipitis_589
contig_5        1260000 109     P.stipitis_589
contig_5        1261000 111     P.stipitis_589
contig_5        1262000 104     P.stipitis_589
contig_5        1263000 108     P.stipitis_589
contig_5        1264000 103     P.stipitis_589
contig_5        1265000 148     P.stipitis_589
contig_5        1266000 103     P.stipitis_589
contig_5        1267000 107     P.stipitis_589
contig_5        1268000 100     P.stipitis_589
contig_5        1269000 101     P.stipitis_589
contig_5        1270000 118     P.stipitis_589
contig_5        1271000 106     P.stipitis_589
contig_5        1272000 102     P.stipitis_589
contig_5        1273000 114     P.stipitis_589
contig_5        1274000 115     P.stipitis_589
contig_5        1275000 116     P.stipitis_589
contig_5        1276000 134     P.stipitis_589
contig_5        1277000 82      P.stipitis_589
contig_5        1278000 74      P.stipitis_589

#7 Peak at 591 seems widther in contig 9.

File=$(ls alignment/genome_alignment/vs_589_unmasked/P.stipitis_591/P.stipitis_591_vs_589_unmasked_depth_10kb.tsv)
cat $File | grep 'contig_9' | less -S

contig_9        581000  86      P.stipitis_591
contig_9        582000  75      P.stipitis_591
contig_9        583000  84      P.stipitis_591
contig_9        584000  84      P.stipitis_591
contig_9        585000  82      P.stipitis_591
contig_9        586000  80      P.stipitis_591
contig_9        587000  91      P.stipitis_591
contig_9        588000  104     P.stipitis_591
contig_9        589000  66      P.stipitis_591
contig_9        590000  145     P.stipitis_591
contig_9        591000  122     P.stipitis_591
contig_9        592000  79      P.stipitis_591
contig_9        593000  129     P.stipitis_591
contig_9        594000  158     P.stipitis_591
contig_9        595000  82      P.stipitis_591
contig_9        596000  69      P.stipitis_591
contig_9        597000  125     P.stipitis_591
contig_9        598000  170     P.stipitis_591
contig_9        599000  148     P.stipitis_591
contig_9        600000  53      P.stipitis_591
contig_9        601000  70      P.stipitis_591
contig_9        602000  72      P.stipitis_591
contig_9        603000  82      P.stipitis_591
contig_9        604000  79      P.stipitis_591
contig_9        605000  80      P.stipitis_591
contig_9        606000  88      P.stipitis_591
contig_9        607000  71      P.stipitis_591
contig_9        608000  92      P.stipitis_591
contig_9        609000  81      P.stipitis_591
contig_9        610000  71      P.stipitis_591
contig_9        611000  85      P.stipitis_591


```
