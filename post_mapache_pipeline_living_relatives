#!/bin/bash
source /home/##username##/miniconda3/etc/profile.d/conda.sh #provide the proper path

###This script should be run in the project main directory (e.g. here is named Adele_project), but be deposited in another directory within the former (e.g. here is named Adele_project/scripts).
###The main directory should also include the directory in which mapache was run (e.g. here is named mapache), which in turn will include the config (e.g. here is named mapache/config) and the results of mapache (e.g. here is named mapache/results_adele_modern), as well as an empty directory named post_analyses.

###script runs as: bash post_mapache_scripts_living_relatives 
###It requires installed software and conda environments and editing the script to set the proper paths manually.
###It requires the indexed hs37d5.fa reference genome downloaded from the 1000 Genomes project ftp server.
###For downstream analyses please manually make the sample_IDs_mapache_males.txt files"

dataset="name_of_dataset" #provide the name of the suffix given in the mapache config file for output directory name: eg adele_modern for results_adele_modern
software_path="sofware_directory_path" #provide the path where the following software is installed: angsd, contamMix, EigenStratDatabaseTools, predict_haplogroup.py)
hs37d5_path="hs37d5_fasta_path" #provide the path of hs37d5.fa reference genome
HapMapChrX_path="HapMapChrX.gz_path" #provide the path of HapMapChrX.gz (see STEP_05 for details)
5M="Koptekin.5M.auto_path" #provide the path of Koptekin.5M.auto files (snp, pos, bed) are located (see STEP_09 for details)
1240K="1240K_path" #provide the path of 1240K files (snp, pos, bed) are located (see STEP_10 for details)

###STEP_01
###This script creates the sample/library IDs for post-mapache analyses.
###Requires a list of sample_level_IDs (sample_IDs_mapache.txt; one name per line) located in mapache/post_analyses (directory post_analyses should be manually made)
###Check the lists to remove potential duplicates.

cd mapache/post_analyses
tail -n +2 ../config/samples_$dataset.tsv | awk 'BEGIN{FS=" "} {printf $1 "\n"}' > sample_IDs_mapache_temp.txt
sort -u sample_IDs_mapache_temp.txt > sample_IDs_mapache.txt
rm sample_IDs_mapache_temp.txt
for sample_level_ID in $(< sample_IDs_mapache.txt); do
  grep "$sample_level_ID" ../config/samples_$dataset.tsv | awk 'BEGIN{FS=" "} {printf $2 "\n"}' >> $sample_level_ID.library_IDs_mapache_temp.txt
  sort -u $sample_level_ID.library_IDs_mapache_temp.txt > $sample_level_ID.library_IDs_mapache.txt
  rm $sample_level_ID.library_IDs_mapache_temp.txt
  grep "$sample_level_ID" ../config/samples_$dataset.tsv |grep ",DS_"| awk 'BEGIN{FS=" "} {printf $2 "\n"}' >> $sample_level_ID.library_IDs_mapache_DS_temp.txt ### This requires providing the "DS" group tag in mapache samplelist file.
  sort -u $sample_level_ID.library_IDs_mapache_DS_temp.txt > $sample_level_ID.library_IDs_mapache_DS.txt
  rm $sample_level_ID.library_IDs_mapache_DS_temp.txt  
  done
echo "!!! Important: For downstream analyses please manually make the ${population}_IDs_READ.txt, ${population}_IDs_KIN.txt and sample_IDs_mapache_males.txt files" 

###STEP_02
###The following script performs genetic sex determination following the Ry method that is based on the ratio of reads mapped to Y and to X and Y chromosomes
###The results (txt files) can be found in the directory mapache/post_analyses/00_2_genetic_sex_Ry.
###It requires ry_compute.py from sup.info of https://www.sciencedirect.com/science/article/pii/S0305440313002495
###It requires the conda environment of mapache installed

echo "!!!Program starts: Genetic sex inference with Ry method!!!"
mkdir 00_2_genetic_sex_Ry
cd 00_2_genetic_sex_Ry
conda activate mapache
for sample_level_ID in $(< ../sample_IDs_mapache.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  samtools view -q 30 ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam | python2 ~/software/Ry_compute/ry_compute.py > $sample_level_ID.Ry.txt
  echo "complete"
  cd ../
  done
conda deactivate
echo "!!!End of program: Genetic sex inference with Ry method!!!"

###STEP_03
###This script is used for contamination estimation using a method based on mtDNA (contamMix).
###The results (txt files) can be found in the directory mapache/post_analyses/01_1_contamination_mtDNA_contamMix.
###It requires having contamMix installed and having mt311.fa downloaded from https://github.com/mpieva/mapping-iterative-assembler/blob/master/misc/mt311.fa
###It requires having ANGSD installed
###It requires the conda environment of mapache installed
###It requires a conda env (here named as adna1) with mafft installed

echo "!!!End of program: mtDNA contamination estimation with contamMix!!!"
mkdir 01_1_contamination_mtDNA_contamMix_samplelevel
cd 01_1_contamination_mtDNA_contamMix_samplelevel
for sample_level_ID in $(< ../sample_IDs_mapache.txt); do
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  conda activate mapache
  samtools view -hb ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam MT -o $sample_level_ID.MT.bam
  samtools index $sample_level_ID.MT.bam
  conda deactivate  
  ~/$software_path/angsd/angsd -i $sample_level_ID.MT.bam -doCounts 1 -minMapQ 30 -minQ 30 -doFasta 2 -out $sample_level_ID.MT.depth1.majcons 
  gunzip $sample_level_ID.MT.depth1.majcons.fa.gz
  cat ~/$software_path/contamMix/mt311.fa $sample_level_ID.MT.depth1.majcons.fa > $sample_level_ID.311.fasta #set proper path for location of mt311.fa
  conda activate mapache
  bwa index $sample_level_ID.MT.depth1.majcons.fa
  conda deactivate
  conda activate adna1 #a conda env with mafft installed
  mafft --auto $sample_level_ID.311.fasta > $sample_level_ID.311.MAFFT.fasta
  conda deactivate
  Rscript /usr/local/lib/R/site-library/contamMix/exec/estimate.R --samFn $sample_level_ID.MT.bam  --malnFn $sample_level_ID.311.MAFFT.fasta --consId MT --figure $sample_level_ID.contamMix_fig | tee $sample_level_ID.contamMixout.txt
  readsused=`grep "consist of" "$sample_level_ID.contamMixout.txt" | awk '{print $3}'`
  map=`grep "MAP authentic" "$sample_level_ID.contamMixout.txt" | cut -d":" -f2`
  quantiles=`awk 'NR==9 {print $0}' "$sample_level_ID.contamMixout.txt"`
  pr=`awk 'NR==4 {print $0}' "$sample_level_ID.contamMixout.txt" | cut -d"(" -f2 | cut -d")" -f1`	 
  err=`grep "error rate" "$sample_level_ID.contamMixout.txt" | awk '{print $9}' | sed 's/).//g'`
  echo -e "Reads used \t authentic \t quantiles \t probability \t error rate" >> $sample_level_ID.tempcmix
  echo -e "$readsused \t $map \t $quantiles \t $pr \t $err" >> $sample_level_ID.tempcmix
  rm $sample_level_ID.311.fasta
  cd ../
  done
cat */*.tempcmix > all_samplelevel.tempcmix
echo "!!!End of program: mtDNA contamination estimation with contamMix!!!"

###STEP_04
###This script is used for contamination estimation using the method based on chromosome X in males (software ANGSD)
###The results (txt files) can be found in the directory mapache/post_analyses/01_3_contamination_chrXmales_angsd.
###It requires having ANGSD installed
###It requires the HapMap file for chrX (HapMapChrX.gz) provided by ANGSD.
###It requires a list (sample_IDs_mapache_males.txt) manually made by copying sample_IDs_mapache and retaining only the males based on sex estimation results.

echo "!!!Program starts: chrX-in-males contamination estimation with angsd!!!"
mkdir 01_3_contamination_chrXmales_angsd
cd 01_3_contamination_chrXmales_angsd
cp ~/$HapMapChrX_path/HapMapChrX.gz ./
for sample_level_ID in $(< ../sample_IDs_mapache_males.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  ~/$software_path/angsd/angsd -i ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam -r X:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out $sample_level_ID
  ~/$software_path/angsd/misc/contamination -a $sample_level_ID.icnts.gz -h ../HapMapChrX.gz -p 32 -s 32 2> $sample_level_ID.contamination.txt
  echo "complete"
  cd ../
  done
rm HapMapChrX.gz
echo "!!!End of program: chrX-in-males contamination estimation with angsd!!!"

###STEP_05
###This script constructs the consensus mtDNA
###The results (txt and pdf) can be found in the directory mapache/post_analyses/02_1_mtDNA_haplogroup.
###It requires having ANGSD installed
###It requires the conda environment of mapache installed

echo "!!!Program starts: mtDNA consensus generation!!!"
mkdir 02_1_mtDNA_haplogroup
cd 02_1_mtDNA_haplogroup
for sample_level_ID in $(< ../sample_IDs_mapache.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  conda activate mapache
  samtools view -hb ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam MT -o $sample_level_ID.MT.bam
  conda deactivate
  ~/$software_path/angsd/angsd -P 100 -i $sample_level_ID.MT.bam -doCounts 1 -minMapQ 30 -minQ 30 -setMinDepth 1 -doFasta 2 -out ${sample_level_ID}_MT_depth1.majcons
  gunzip ${sample_level_ID}_MT_depth1.majcons.fa.gz
  sed -i -e "s/MT/${sample_level_ID}_MT_depth1.majcons/g" ${sample_level_ID}_MT_depth1.majcons.fa
  echo "complete"
  cd ../
  done
echo "fasta files can be used for haplogrpup assignments with the online tools of HaploCart and HaploGrep"  
echo "!!!End of program: mtDNA haplogroup assignment using HaploGrep!!!"
    
###STEP_06
###The script below is using the software Yleaf to assign the male human samples to chromosome-Y haplogroups.
###The results (txt files) can be found in the directory mapache/post_analyses/03_1_chrY_haplogroup_Yleaf.
###It requires having yleaf installed via conda
###It also requires predict_haplogroup.py of yleaf software

echo "!!!Program starts: chrY haplogroup assignment using Yleaf!!!"
mkdir 03_1_chrY_haplogroup_Yleaf
cd 03_1_chrY_haplogroup_Yleaf
conda activate yleaf
for sample_level_ID in $(< ../sample_IDs_mapache_males.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  Yleaf -bam ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam -o $sample_level_ID.depth1 -q 30 -r 1 -b 90 -t 16 -rg hg19 -dh
  python ~/$software_path/predict_haplogroup.py -i ./$sample_level_ID.depth1 -o $sample_level_ID.depth1.hg
  echo "complete"
  cd ../
  done
conda deactivate
echo "End of program: chrY haplogroup assignment using Yleaf!!!"

###STEP_07
###The script below is using the software Yhaplo to assign the male human samples to chromosome-Y haplogroups.
###The results (vcf and txt files) can be found in the directory mapache/post_analyses/03_1_chrY_haplogroup_Yhaplo.
###It requires the conda environment of mapache installed
###It requires the indexed hs37d5.fa reference genome downloaded from NCBI.

echo "!!!Program starts: chrY haplogroup assignment using Yhaplo!!!"
mkdir 03_1_chrY_haplogroup_Yhaplo
cd 03_1_chrY_haplogroup_Yhaplo
for sample_level_ID in $(< ../sample_IDs_mapache_males.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  conda activate mapache
  bcftools mpileup -q 30 -Q 30 -C 50 -Ou -f ~/$hs37d5_path/hs37d5.fa ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam -r Y --ignore-RG | bcftools call --ploidy 1 -m -Ov -o $sample_level_ID.Y.vcf
  bcftools norm -f ~/$hs37d5_path/hs37d5.fa $sample_level_ID.Y.vcf -Ov -o $sample_level_ID.Y.norm.depth1.vcf
  conda deactivate
  yhaplo -i $sample_level_ID.Y.norm.depth1.vcf -o $sample_level_ID.Y.norm.depth1 -aao
  echo "complete"
  cd ../
  done
echo "End of program: chrY haplogroup assignment using Yhaplo!!!"

###STEP_08
### This script performs SNP calling (with SAMtools) for the positions of the 5M panel (Koptekin et al. 2023) and pseudohaploidization (with pileupcaller)
###The results (EIGENSTRAT files) can be found in the directory mapache/post_analyses/06_1_pseudohaplo_call_5M.
###It requires the conda environment of mapache installed
###It requires a conda env (here named as adna1) with pileupCaller and mergeit installed
###It requires having EigenStratDatabaseTools installed.
###It requires the indexed hs37d5.fa reference genome downloaded from NCBI.
###It requires the 5M_auto SNPs list as a position file (Koptekin.5M.auto.pos) but also in EIGENSTRAT format (Koptekin.5M.auto.snp; can be downloaded from their study https://www.sciencedirect.com/science/article/pii/S0960982222018243)

echo "!!!Program starts: 5M call and pseudohaploidization!!!"
mkdir 06_1_pseudohaplo_call_5M
cd 06_1_pseudohaplo_call_5M
for sample_level_ID in $(< ../sample_IDs_mapache.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  conda activate mapache 
  samtools mpileup -B -R -q 30 -Q 30 -f ~/$hs37d5_path/hs37d5.fa ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam -l ~/$5M/Koptekin.5M.auto.pos -o $sample_level_ID.5M.pileup
  conda deactivate
  conda activate adna1 #a conda env with pileupCaller and mergeit installed
  pileupCaller --randomHaploid --sampleNames $sample_level_ID --samplePopName $sample_level_ID -f ~/$5M/Koptekin.5M.auto.snp -e $sample_level_ID.5M < $sample_level_ID.5M.pileup
  rm $sample_level_ID.5M.pileup
  echo "complete"
  cd ../
  done
tail -n +2 ../sample_IDs_mapache.txt > temp.all_individuals_but_first.txt
first_individual=$(head -n 1 ../sample_IDs_mapache.txt)
cp $first_individual/$first_individual.5M.snp $dataset.5M.snp
cp $first_individual/$first_individual.5M.geno $dataset.5M.geno
cp $first_individual/$first_individual.5M.ind $dataset.5M.ind
for individual in $(< temp.all_individuals_but_first.txt); do
  echo "using Sample" $individual
  touch temp.merge.par
  echo "geno1:" $individual/$individual.5M.geno > temp.merge.par
  echo "snp1:" $individual/$individual.5M.snp >> temp.merge.par
  echo "ind1:" $individual/$individual.5M.ind >> temp.merge.par
  echo "geno2:" $dataset.5M.geno >> temp.merge.par
  echo "snp2:" $dataset.5M.snp >> temp.merge.par
  echo "ind2:" $dataset.5M.ind >> temp.merge.par 
  echo "genooutfilename:" temp.geno >> temp.merge.par
  echo "snpoutfilename:" temp.snp >> temp.merge.par
  echo "indoutfilename:" temp.ind >> temp.merge.par
  mergeit -p temp.merge.par
  mv temp.geno $dataset.5M.geno
  mv temp.snp $dataset.5M.snp
  mv temp.ind $dataset.5M.ind
  rm temp.merge.par
  done
touch temp.convert.par
echo "genotypename:"	$dataset.5M.geno > temp.convert.par
echo "snpname:"	$dataset.5M.snp >> temp.convert.par
echo "indivname:"	$dataset.5M.ind >> temp.convert.par
echo "outputformat:	PACKEDPED" >> temp.convert.par
echo "genotypeoutname:"	$dataset.5M.bed >> temp.convert.par
echo "snpoutname:"	$dataset.5M.bim >> temp.convert.par 
echo "indivoutname:"	$dataset.5M.fam >> temp.convert.par
convertf -p temp.convert.par
rm temp.convert.par
touch temp.convert.par
echo "genotypename:" $dataset.5M.geno > temp.convert.par
echo "snpname:"     $dataset.5M.snp >> temp.convert.par
echo "indivname:"   $dataset.5M.ind >> temp.convert.par
echo "outputformat:"    EIGENSTRAT >> temp.convert.par
echo "genotypeoutname:" $dataset.5M.unpacked.geno >> temp.convert.par
echo "snpoutname:"      $dataset.5M.unpacked.snp >> temp.convert.par
echo "indivoutname:"    $dataset.5M.unpacked.ind >> temp.convert.par
convertf -p temp.convert.par
rm temp.convert.par
conda deactivate
~/$software_path/EigenStratDatabaseTools/eigenstrat_snp_coverage.py -g $dataset.5M.unpacked.geno -s $dataset.5M.unpacked.snp -i $dataset.5M.unpacked.ind -o $dataset.5M
rm *unpacked*
echo "!!! Important!: Fill the sex field in the ind file using the genetic sex inference results. !!!"
echo " "
echo "!!!End of program: 5M call and pseudohaploidization!!!"

###STEP_09
### This script performs SNP calling (with SAMtools) for the positions of the 1240K panel and pseudohaploidization (with pileupcaller)
###The results (EIGENSTRAT files) can be found in the directory mapache/post_analyses/06_1_pseudohaplo_call_1240K.
###It requires the conda environment of mapache installed
###It requires a conda env (here named as adna1) with pileupCaller and mergeit installed
###It requires having EigenStratDatabaseTools installed.
###It requires the indexed hs37d5.fa reference genome downloaded from NCBI.
###It requires the 1240K SNPs list as a position file (1240K.pos) but also in EIGENSTRAT format (1240K.snp; can be downloaded from their study https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10104067/)


echo "!!!Program starts: 1240K call and pseudohaploidization!!!"
mkdir 06_1_pseudohaplo_call_1240K
cd 06_1_pseudohaplo_call_1240K
for sample_level_ID in $(< ../sample_IDs_mapache.txt); do  
  mkdir $sample_level_ID
  cd $sample_level_ID
  echo "using sample" $sample_level_ID
  conda activate mapache 
  samtools mpileup -B -R -q 30 -Q 30 -f ~/$hs37d5_path/hs37d5.fa ../../../results_$dataset/03_sample/03_final_sample/01_bam/$sample_level_ID.hs37d5.bam -l ~/$1240K/1240K.pos -o $sample_level_ID.1240K.pileup
  conda deactivate
  conda activate adna1 #a conda env with pileupCaller and mergeit installed
  pileupCaller --randomHaploid --sampleNames $sample_level_ID --samplePopName $sample_level_ID -f ~/$1240K/1240K.snp -e $sample_level_ID.1240K < $sample_level_ID.1240K.pileup
  rm $sample_level_ID.1240K.pileup
  echo "complete"
  cd ../
  done
tail -n +2 ../sample_IDs_mapache.txt > temp.all_individuals_but_first.txt
first_individual=$(head -n 1 ../sample_IDs_mapache.txt)
cp $first_individual/$first_individual.1240K.snp $dataset.1240K.snp
cp $first_individual/$first_individual.1240K.geno $dataset.1240K.geno
cp $first_individual/$first_individual.1240K.ind $dataset.1240K.ind
for individual in $(< temp.all_individuals_but_first.txt); do
  echo "using Sample" $individual
  echo "geno1:" $individual/$individual.1240K.geno > temp.merge.par
  echo "snp1:" $individual/$individual.1240K.snp >> temp.merge.par
  echo "ind1:" $individual/$individual.1240K.ind >> temp.merge.par
  echo "geno2:" $dataset.1240K.geno >> temp.merge.par
  echo "snp2:" $dataset.1240K.snp >> temp.merge.par
  echo "ind2:" $dataset.1240K.ind >> temp.merge.par 
  echo "genooutfilename:" temp.geno >> temp.merge.par
  echo "snpoutfilename:" temp.snp >> temp.merge.par
  echo "indoutfilename:" temp.ind >> temp.merge.par
  mergeit -p temp.merge.par
  mv temp.geno $dataset.1240K.geno
  mv temp.snp $dataset.1240K.snp
  mv temp.ind $dataset.1240K.ind
  rm temp.merge.par
  done
touch temp.convert.par
echo "genotypename:"	$dataset.1240K.geno > temp.convert.par
echo "snpname:"	$dataset.1240K.snp >> temp.convert.par
echo "indivname:"	$dataset.1240K.ind >> temp.convert.par
echo "outputformat:	PACKEDPED" >> temp.convert.par
echo "genotypeoutname:"	$dataset.1240K.bed >> temp.convert.par
echo "snpoutname:"	$dataset.1240K.bim >> temp.convert.par 
echo "indivoutname:"	$dataset.1240K.fam >> temp.convert.par
convertf -p temp.convert.par
rm temp.convert.par
touch temp.convert.par
echo "genotypename:" $dataset.1240K.geno > temp.convert.par
echo "snpname:"     $dataset.1240K.snp >> temp.convert.par
echo "indivname:"   $dataset.1240K.ind >> temp.convert.par
echo "outputformat:"    EIGENSTRAT >> temp.convert.par
echo "genotypeoutname:" $dataset.1240K.unpacked.geno >> temp.convert.par
echo "snpoutname:"      $dataset.1240K.unpacked.snp >> temp.convert.par
echo "indivoutname:"    $dataset.1240K.unpacked.ind >> temp.convert.par
convertf -p temp.convert.par
rm temp.convert.par
conda deactivate
~/$software_path/EigenStratDatabaseTools/eigenstrat_snp_coverage.py -g $dataset.1240K.unpacked.geno -s $dataset.1240K.unpacked.snp -i $dataset.1240K.unpacked.ind -o $dataset.1240K
rm *unpacked*
echo "!!! Important!: Fill the sex field in the ind file using the genetic sex inference results. !!!"
echo " "
echo "!!!End of program: 1240K call and pseudohaploidization!!!"

cd ../../../../