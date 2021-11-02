# Pancreatic Cancer Microbiome: Call Microbial SNVs

*Project*: PanCan
*Author*: Sebastian Schmidt
*Contact*: sebastian.schmidt@embl.de
*Date*: 2018-09-25
*Version*: v01

```
export v=v01
```

## Introduction

To quantify the extent of oral-fecal microbial transmission among pancreatic cancer patients and controls, and to demarcate potential sub-species in the dataset, a first step is to call microbial SNVs on both oral and fecal samples. This can be achieved using `[metaSNV](http://metasnv.embl.de)` by following these steps, starting from unique mappings to reference genomes:

* compute the horizontal (*breadth*) and vertical (*depth*) coverage for each genome in each sample.
* perform SNV calls
* filter called SNVs by lenient criteria, splitting calls by genome

Before starting, create a metaSNV folder within the project folder.

```
folderBase=/g/bork1/tschmidt/pancreatic_cancer_microbiome/
cd $folderBase

/g/bork1/tschmidt/bin/metaSNP/metaSNP_NEW metaSNV_$v
```

## Compute per-Sample Coverages

First, merge relevant bam files that are technical replicates (re-sequencing of samples)

```
module add samtools

echo "MMPW67722444ST_MMPC63233635ST
MMPW80230725ST_MMPC21059026ST
MMPW14699153ST_MMPC49013463ST
MMPW92389414ST_MMPC68480658ST
MMPW83144761ST_MMPC43926327ST
MMPW39199901ST_MMPC53747560ST
MMPW76972986ST_MMPC49776422ST
MMPW24306397ST_MMPC96048560ST
MMPW28983031ST_MMPC45272194ST
MMPW21548424ST_MMPC73667956ST
MMPW31043630ST_MMPC23673679ST
MMPW97100755ST_MMPC10431544ST
MMPW29397139ST_MMPC14211565ST
MMPW82409611ST_MMPC61634104ST" | while read pair; do
  mmpc=`echo $pair | cut -d "_" -f 2`
  mmpw=`echo $pair | cut -d "_" -f 1`

  mmpcName=/g/bork1/tschmidt/pancreatic_cancer_microbiome/data/freeze11.repgenomes/freeze11.UL/$mmpc.ULRepGenomesv11.unique.sorted.bam
  mmpwName=/g/bork5/mocat/ngless-processing/repgenomes11UL/2019-07-05_MMPW/outputs/$mmpw.fr11Repv2UL.unique.sorted.bam
  mergedName=/g/bork1/tschmidt/pancreatic_cancer_microbiome/data/freeze11.repgenomes/freeze11.UL/merged/$mmpc.merged.unique.sorted.bam

  if [ -f $mmpwName ]; then
    echo "Merge $mmpc and $mmpw to $mergedName"

    samtools merge $mergedName $mmpcName $mmpwName && ln -s $mergedName /g/bork1/tschmidt/pancreatic_cancer_microbiome/data/freeze11.repgenomes/freeze11.UL/ && rm $mmpcName &
  fi
done
```

First, define the list of mapping `bam` files to operate on:

```
file_bam_list=$folderBase/parameters/metaSNV.bam_list.$v.txt

ls -d $folderBase/data/freeze11.repgenomes/freeze11.UL/*.bam > $file_bam_list

wc -l $file_bam_list
240 /g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaSNV.bam_list.v01.txt
```

Next, prepare and submit coverage computation run.

```
echo "% module add metaSNV" > metaSNV_$v/run.coverage
/g/bork1/tschmidt/bin/metaSNP/metaSNP_COV metaSNV_$v $file_bam_list | perl -pe 's/\/g\/bork1\/tschmidt\/bin\/metaSNP\//\/g\/scb2\/bork\/mocat\/software\/metaSNV\/git-49fec2dc\//g' >> metaSNV_$v/run.coverage

/g/bork1/tschmidt/embltilities/bin/submitjob -n PanCan.cov_comp -l logs/$v.compute_coverage.log metaSNV_$v/run.coverage
```

Collapse and summarise coverage.

```
/g/bork1/tschmidt/bin/metaSNV/wrap.compute_summary.py metaSNV_v01/
```

Split into bed files per genome.

```
cat /g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.txt | while read id; do grep "^$id\." /g/bork1/tschmidt/db/to_freezer/freezes/fr11/v2/prok-refdb-v11.0.0_specI-v2_representatives-v2UL_contigs-v1.bed > /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/split_by_tax_id/$id.bed; done
```

## Call Microbial SNVs

SNV calling is performed one genome at a time. First, define SNV calling parameters.

```
#Minimum coverage at position (across all samples)
PARAM_c=10

#Minimum number of non-reference reads to call SNV
PARAM_t=2

#Minimum frequency of non-ref reads to call a "population SNV"
#=> this effectively retains all SNVs as "population SNVs"
PARAM_p=0.001

#Maximum depth for mpileup
PARAM_d=2500
```

Generate job scripts to submit for processing.

```
file_job=$folderBase/metaSNV_$v/run.snp_calls
echo "% module add samtools/1.5" > $file_job
echo "% module add metaSNV/1.0.3" >> $file_job

ref_db=/scratch/rossum/freezer/prok-refdb/v11/v2UL/prok-refdb-v11.0.0_specI-v2_representatives-v2UL_contigs-v1.fna
ref_anno=/scratch/rossum/freezer/prok-refdb/v11/freeze11.annotations.txt
#  -a /g/bork1/tschmidt/bin/metaSNP/db/freeze9.annotations.txt \

cat /g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.txt | while read id;
do
  echo "samtools mpileup -l $folderBase/metaSNV_$v/split_by_tax_id/$id.bed -f $ref_db -B -b $file_bam_list | /g/scb2/bork/mocat/software/metaSNV/1.0.3/src/snpCaller/snpCall -f $ref_db  -i $folderBase/metaSNV_$v/snpCaller/$id.called_indiv -c $PARAM_c -t $PARAM_t -p $PARAM_p > $folderBase/metaSNV_$v/snpCaller/$id.called_SNVs" >> $file_job
done
```

Submit jobs.

```
module add submitjobs

submitjob -a 32 -m 32 -n PanCan.SNV -l $folderBase/logs/$v.snv_calls.per_genome.fixed.log $file_job
```


```
/g/scb2/bork/rossum/metaSNV-debug/myTest/runAllBams/runWithSebastiansParams/metaSNV-1.0.3-sebaParams/metaSNV_sebaParams.py --print-commands --n_splits 12 --threads 20 /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_v01/ /g/bork1/tschmidt/pancreatic_cancer_microbiome//parameters/metaSNV.bam_list.v01.txt /scratch/rossum/freezer/prok-refdb/v11/v2UL/prok-refdb-v11.0.0_specI-v2_representatives-v2UL_contigs-v1.fna

```

## Filter Results for Relevant SNVs

Only very lenient filters are applied, as most filtering is done during downstream analyses. Set these parameters:

```
#Vertical coverage
#=> avg depth per genome per sample
allCov=/g/scb2/bork/rossum/metaSNV-debug/myTest/runAllBams/runWithSebastiansParams/outputs/outputs.all_cov.tab

#Horizontal coverage
#=> percentage of genome covered at �1x per sample
percCov=/g/scb2/bork/rossum/metaSNV-debug/myTest/runAllBams/runWithSebastiansParams/outputs/outputs.all_perc.tab

#BAM files list
file_sample_list=$folderBase/parameters/metaSNV.sample_list.$v.txt
cat $file_bam_list | cut -d "/" -f 9 | perl -pi -e 's/.ULRepGenomesv11.unique.sorted.bam//' > $file_sample_list

#Horizontal coverage cutoff
#=> for each genome in each sample, set the minimum required horizontal coverage
b=5

#Vertical coverage cutoff
#=> for each genome in each sample, set the minimum required average depth
d=0.001

#Minimum sample criterion
#=> for each genome, in how many samples are above criteria required to be satisfied?
m=2

#SNP coverage cutoff
#=> for each position across a genome, set a minimum coverage cutoff (per sample)
c=1

#SNP sample incidence cutoff
#=> for each position across a genome, set a minimum ratio of samples of interest in which SNP is required to be observed
incidence=0

#Set path to filtering script
FILTERING=/g/bork1/tschmidt/bin/metaSNP/metaSNP_filtering.py
```

Submit jobs to filter called SNVs.

```
joblist=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/$v.filter_SNV.jobs
echo "" > $joblist
for i in $(seq 0 39); do
  echo "/g/bork3/x86_64/bin/python $FILTERING -b $b -d $d -m $m -c $c -p $incidence $percCov $allCov /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/snpCaller/called_SNVs.best_split_$i $file_bam_list /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/filtered.pop/" >> $joblist
  echo "/g/bork3/x86_64/bin/python $FILTERING -b $b -d $d -m $m -c $c -p $incidence $percCov $allCov /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/snpCaller/indiv_called.best_split_$i $file_bam_list /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/filtered.indiv/" >> $joblist
done

/g/bork1/tschmidt/embltilities/bin/submitjob -m 4 -n PanCan.filter_SNV -l logs/$v.snv_filtering.log $joblist
```

Tidy up.
```
gzip /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/snpCaller/called_SNVs* &
gzip /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/snpCaller/indiv_called* &
gzip /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/filtered.pop/* &
gzip /g/bork1/tschmidt/pancreatic_cancer_microbiome/metaSNV_$v/filtered.indiv/* &
```
