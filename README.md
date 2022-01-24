# Structurome Analysis with Lead-Seq

Lead-seq is a transcriptome-wide structure probing approach, which combines structure probing with lead(II)-ions with high-throughput sequencing.

The work has been published as [Twittenhoff, Brandenburg, Righetti et al., _NAR_ 2021](https://doi.org/10.1093/nar/gkaa404), where Lead-Seq was used to analyze the structurome of _Yersinia pseudotuperculosis_ at 25 °C and 37 °C. The raw sequencing data can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140649).

## sample IDs

| ID  | sample                 |
| --- | ---------------------- |
| 01  | 25 °C &ensp; control   |
| 02  | 25 °C &ensp; treatment |
| 03  | 37 °C &ensp; control   |
| 04  | 37 °C &ensp; treatment |

# Workflow

The results of Step 2.(Count 5'ends) and all following steps are provided here.

<br>

## Mapping

* Download the raw sequence reads from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140649).
* Clean the reads, e.g. with [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

* Prepare the reference transcriptome, eg with split2transcript.py.   The _Y. pseudotub._ transcriptome can be found in data/transcriptome.fa.   
  
  If you want to map to the genome instead of the transcriptome, you can use split2transcrips.py to **split your lead-scores into transcriptome-format**. 

```
usage: split2transcripts.py genome_or_scores gff outfile [options]

options:
  -f  gff-ID to filter for. Default is 'CDS'
  -n  name to use from gff attributes. Default is 'locus_tag'
```

* Map the sequences to the reference transcriptome, e.g. with [bowtie](http://bowtie-bio.sourceforge.net/index.shtml).

<br>

## Generate Raw Counts

To generate raw counts by counting the 5'ends of mapped reads, run the script _sam2tab.pl_ for each mapping. The script is published by [Wan et al., 2013](https://doi.org/10.1038/nprot.2013.045) and can be downloaded [here](https://s3.amazonaws.com/changbackup/ywan/PARS_Nature_Protocols/sample_data.tar.gz).

The parameters _X Y Z_ are defined as the number of bases to trim from the 5'end (_X_ ) and 3'end (_Y_ ), and the number of threads (_Z_ ).
**For Lead-Seq, _X_ must be 0.** The resulting .tab-files can be found in tab_files/.

<br>

## Calculate Lead-Scores

To normalize the raw counts and calculate Lead-Scores from them, run calc_scores.py.

For the _Y. pseudotub._ samples, instead of using the read depth as integer, you can use these shortcuts:

| sample-ID | 01  | 02  | 03  | 04  |
| --------- | --- | --- | --- | --- |
| shortcut  | sd1 | sd2 | sd3 | sd4 |

```
usage: calc_scores.py control.tab treatment.tab read_depth_control read_depth_treatment outfile

example:
>> calc_scores.py tab_files/01.tab tab_files/02.tab sd1 sd2 l25.txt
>> calc_scores.py tab_files/03.tab tab_files/04.tab sd3 sd4 l37.txt
```

<br>

## Calculate differences between samples

Calculate Δscores and Mann-Whitney-U with calc_delta_mwu.py.

```
usage: calc_delta_mwu.py score_cond1 score_cond2 outfile_Δscores outfile_p-values

example:
>> calc_delta_mwu.py l25.txt l37.txt delta.txt mwu.txt
```
Lead-scores have to be in transcriptome-format for this. If you did not map to the transcriptome, you can use split2transcrips.py to split your lead-scores into transcriptome-format. 


<br>

## Search for RNA thermometer candidates

The results around the Ribosome Binding site can then be extracted with the 5'UTR positions notated in data/5UTR_annotation.txt and get_SDS.py.
```
usage: get_SDS.py Δscores mwu 5'UTRs outfile

example:
>> get_SDS.py delta.txt mwu.txt data/5UTR_annotation.txt SDS.txt
```
A plot with p-values against Δscores will be plotted to plots/.

<br>

If needed, the 5'UTR-file can be prepared from a gff-file with gff2utr.py

```
usage: gff2utr.py gff outfile_UTR [options]

options:
  -f  gff-ID to filter for. Default is '5UTR'
  -n  name to use from gff attributes. Default is 'locus_tag'
```

<br>

## Validate on tRNAs

You can run a validation of the lead-scores on known tRNA structures from the file data/tRNA_testdata.txt by running validate_on_tRNAs.py.

```
usage: validate_on_tRNAs.py score_cond1 scores_cond2 tRNA_testdata

example:
>> validate_on_tRNAs.py l25.txt l37.txt data/tRNA_testdata.txt
```

The results will be plotted to plots/.
