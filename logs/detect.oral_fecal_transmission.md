# Pancreatic Cancer Microbiome: Call Microbial SNVs

*Project*: PanCan
*Author*: Sebastian Schmidt
*Contact*: sebastian.schmidt@embl.de
*Date*: 2018-10-06
*Version*: v01

```
export v=v01
```

## Prepare SNV data

In a first step ([oral_fecal.prepare.abundance_data.Rmd](Rmd/oral_fecal.prepare.abundance_data.Rmd)), abundance and coverage data was prefiltered and a list of taxa to process for SNVs was defined. Submit jobs for SNV data prepration for all taxa on that list:

```
module add submitjobs

fileTaxList=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.txt
fileR=/g/bork1/tschmidt/pancreatic_cancer_microbiome/R/job.prepare.data.SNV.R
joblist=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/$v.prepare_SNV_data.joblist.txt

echo "% source /g/bork3/home/tschmidt/.bashrc" > $joblist

#Generate list of jobs
while read taxid; do
  echo "$taxid"
  echo "/g/bork3/home/tschmidt/miniconda3/bin/Rscript $fileR --args $taxid" >> $joblist
done < $fileTaxList

#Submit for parallel processing
submitjob -a 32 -m 32 -n PanCan.prep_SNV -l /g/bork1/tschmidt/pancreatic_cancer_microbiome/logs/$v.prepare_SNV.log $joblist
```

## Detect Oral-Fecal Transmission

The next step is to detect microbial SNV signatures of oral-fecal transmission for all taxa, across all individuals.

```
#Get full list of taxa and filter against already finished runs
fileTaxListPre=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.txt
fileTaxListDone=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/v01.detect_transmission.done
fileTaxList=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.to_run.txt
grep -v -f $fileTaxListDone $fileTaxListPre > $fileTaxList

#Define R script
fileR=/g/bork1/tschmidt/pancreatic_cancer_microbiome/R/job.detect_transmission.alleles.R

#Preallocate job list
joblist=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/$v.detect_transmission.joblist.txt

echo "% source /g/bork3/home/tschmidt/.bashrc" > $joblist

#Generate list of jobs
while read taxid; do
  echo "$taxid"
  echo "/g/bork3/home/tschmidt/miniconda3/bin/Rscript $fileR --args $taxid" >> $joblist
done < $fileTaxList

#Submit for parallel processing
/g/bork1/tschmidt/embltilities/bin/submitjob -s SLURM -a 72 -c 1 -m 32 -n PanCan.detect -k 240:00:00 -l /g/bork1/tschmidt/pancreatic_cancer_microbiome/logs/$v.detect_transmission.log $joblist

/g/bork1/tschmidt/embltilities/bin/submitjob -a 64 -m 32 -n PanCan.detect -l /g/bork1/tschmidt/pancreatic_cancer_microbiome/logs/$v.detect_transmission.log $joblist
```

Consolidate results of per-taxon runs.

```
/g/bork3/home/tschmidt/miniconda3/bin/Rscript /g/bork1/tschmidt/pancreatic_cancer_microbiome/R/consolidate.detect_transmission.alleles.R
```

## Calculate SNV Distances between Samples

```
#fileTaxListPre=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/metaG.oral_fecal.tax_ids_abd.txt

#Define R script
fileR=/g/bork1/tschmidt/pancreatic_cancer_microbiome/R/job.SNV_distances.R

#Preallocate job list
joblist=/g/bork1/tschmidt/pancreatic_cancer_microbiome/parameters/$v.SNV_distances.joblist.txt

echo "% source /g/bork3/home/tschmidt/.bashrc" > $joblist

#Generate list of jobs
while read taxid; do
  echo "$taxid"
  echo "/g/bork3/home/tschmidt/miniconda3/bin/Rscript $fileR --args $taxid" >> $joblist
done < $fileTaxList

#Submit for parallel processing
#/g/bork1/tschmidt/embltilities/bin/submitjob -s SLURM -m 24 -n PanCan.detect -k 72:00:00 -l /g/bork1/tschmidt/pancreatic_cancer_microbiome/logs/$v.detect_transmission.cluster_re.log $joblist

submitjob -m 32 -n PanCan.dist -l /g/bork1/tschmidt/pancreatic_cancer_microbiome/logs/$v.SNV_distances.log $joblist
```
