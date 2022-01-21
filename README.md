# Structurome Analysis with Lead-Seq


Lead-seq is a transcriptome-wide structure probing approach, which combines structure probing with lead(II)-ions with high-throughput sequencing.



The work has been published as [Twittenhoff, Brandenburg, Righetti et al., *NAR* 2021](https://doi.org/10.1093/nar/gkaa404), where Lead-Seq was used to analyze the structurome of *Yersinia pseudotuperculosis* at 25 °C and 37 °C.  The raw sequencing data can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140649).



## sample IDs


| ID | sample      |
|----|-------------- |
| 01 | 25 °C &ensp; control   |
| 02 | 25 °C &ensp; treatment |
| 03 | 37 °C &ensp; control   |
| 04 | 37 °C &ensp; treatment |



# Workflow

The results of Step 2.(Count 5'ends) and all following steps are provided here.
## 1. Mapping

- Download the raw sequence reads from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140649).
- Clean the reads, e.g. with [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

- Prepare the reference transcriptome.  
The *Y. pseudotub.* transcriptome can be found in *data/transcriptome.fa*

- Map the sequences to the reference transcriptome, e.g. with [bowtie](http://bowtie-bio.sourceforge.net/index.shtml). 



## 2. Generate Raw Counts 

To generate raw counts by counting the 5'ends of mapped reads, run the script *sam2tab.pl* for each mapping.
The script is published by [Wan et al., 2013](https://doi.org/10.1038/nprot.2013.045) and can be downloaded [here](https://s3.amazonaws.com/changbackup/ywan/PARS_Nature_Protocols/sample_data.tar.gz).  

The parameters *X Y Z* are defined as the number of bases to trim from the 5'end (*X*) and 3'end (*Y*), and the number of threads (*Z*). 

**For Lead-Seq, *X* must be 0.**

The resulting *.tab* can be found in *./tab_files/*

## 3. Calculate Lead-Scores

To normalize the raw counts and calculate Lead-Scores from them, run *calc_scores.py*.


For the *Y. pseudotub.* samples, instead of using the read depth as integer, you can use these shortcuts:

| sample-ID| 01 | 02 | 03 | 04 |
| ---      | ---| ---| ---| ---|
|shortcut  | sd1| sd2| sd3| sd4|

```
calc_scores.py [control.tab] [treatment.tab] [read depth control] [read depth treatment] [outfile]

>> calc_scores.py tab_files/01.tab tab_files/02.tab sd1 sd2 l25.txt
>> calc_scores.py tab_files/03.tab tab_files/04.tab sd3 sd4 l37.txt
```


## 4. Validate on tRNAs

You can run a validation of the lead-scores on known tRNA structures from the file *data/tRNA_testdata.txt* by running *validate_on_tRNAs.py*

```
validate_on_tRNAs.py [file scores 25°C] [file scores 37°C] [file tRNA testdata]

>> validate_on_tRNAs.py l25.txt l37.txt data/tRNA_testdata.txt
```
The results will be plotted to ./plots/.


## 5. Search for RNA thermometer candidates 

To search for RNA thermometer candidates, first calculate Δscores and Mann-Whitney-U with *calc_delta_mwu.py* 

```
calc_delta_mwu.py [file scores 25°C] [file scores 37°C] [outfile Δscores] [outfile p-values]

>> calc_delta_mwu.py l25.txt l37.txt delta.txt mwu.txt
```

The results around the Ribosome Binding site can then be extracted with the 5'UTR positions notated in *data/5UTR_annotation.txt* and the script *get_SDS.py*

```
get_SDS.py [file Δscores] [file mwu] [file 5'-UTRs] [outfile]

>> get_SDS.py delta.txt mwu.txt data/5UTR_annotation.txt SDS.txt
```

A plot with p-values against Δscores will be plotted to ./plots/.


