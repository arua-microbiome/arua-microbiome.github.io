---
title: "Amplicon analysis with QIIME2"
excerpt: "An example workflow using QIIME2 version 2024.2"
layout: single
author: "George Kalogiannis"

---
By George Kalogiannis, Designed from the official [QIIME2 tutorials](https://docs.qiime2.org/2024.2/tutorials/)
{% include toc %}

# Data sets

In this workshop, we will study microbial genes collected from the soil surrounding the roots of Arabidopsis thaliana plants. These plants were genetically modified to alter how they respond to phosphate stress. Phosphate (P) is a cornerstone nutrient for plants, driving ATP-based energy transfer, nucleic-acid synthesis, membrane formation, and phosphorylation-driven signalling. Because phosphate binds tightly to soil particles, it is frequently scarce, so plants switch on an integrated phosphate-starvation response (PSR): reshaping their roots, releasing phosphatases, and recruiting helpful microbes in the rhizosphere. These changes ripple outward, altering the surrounding microbial community and, in turn, feeding back on plant health.

Our workshop dataset comes from _Arabidopsis thaliana_ mutants studied by Castrillo et al. (2017). Each mutant disables a gene that links phosphate sensing to immunity, creating a “natural experiment” for microbiome shifts. To determine the role of phosphate starvation response in controlling microbiome composition, we will analyse five mutants related to the Pi-transport system (pht1;1, pht1;1;pht1;4, phf1, nla and pho2) and two mutants directly involved in the transcriptional regulation of the P-starvation response (phr1 and spx1;spx2).

The microbiome dataset itself comprises 146 paired-end Illumina MiSeq libraries (2 × 250 bp) that target the hyper-variable V4 region of the bacterial 16S rRNA gene, covering root-zone, bulk-soil and endophytic compartments for each plant line. Technical replicates have been merged, low-quality reads removed, and every sample has been rarefied to 10 000 high-quality reads so that comparisons across genotypes are on an equal footing. 

Alongside the raw FASTQ files you’ll find a manifest (listing file paths and read orientation) and a metadata table that records genotype, compartment, replicate and block information. Together these resources provide a tidy, balanced starting point for generating amplicon-sequence variants (ASVs), profiling community diversity, and pinpointing the microbial taxa that respond to each phosphate-stress mutation.

|Download the data here:|
| ------------- |
| [Raw Sequences]() |
| [Manifest File]() | 


> Castrillo, G., Teixeira, P.J.P.L., Paredes, S.H., Law, T.F., de Lorenzo, L., Feltcher, M.E., Finkel, O.M., Breakfield, N.W., Mieczkowski, P., Jones, C.D., Paz-Ares, J., Dangl, J.L., 2017. Root microbiota drive direct integration of phosphate stress and immunity. Nature 543, 513–518. [doi:10.1038/nature21417](https://dx.doi.org/10.1038/nature21417)


# Genetics Primer

## 1 · What actually comes off a sequencer?

Modern DNA sequencers, Illumina (short, very accurate reads), PacBio HiFi and Oxford Nanopore (long, increasingly accurate reads), do not give you a “genome” or a list of species. They give you a bag of millions of little DNA fragments called reads. These reads are written to disk as plain-text files, almost always compressed and ending in .fastq.gz.

- FASTA files ( .fasta, .fa ) hold only the nucleotide strings; they are mainly used once reads have been cleaned and merged.
- FASTQ files ( .fastq, .fq ) hold both the sequence and its Phred quality scores, telling you how confident the instrument was at each base.

Once you have the downloads in place, peek inside one of the real files, say GC1GC1_R1.fastq.gz (forward reads from replicate GC1, experiment GC1), to convince yourself what a FASTQ record looks like and to check that the quality scores are sensible before you hand the data to QIIME 2. You can do this in bash like this:

```bash
# ── Bash one-liner: show the first two reads (8 lines) ───────────
zcat GC1GC1_R1.fastq.gz | head -n 8
```
Breaking that command down, ```zcat``` prints the unzipped file, ```head -n 8``` takes the top 8 lines of the file.


## 2 · “How many reads did I get?”

A core sanity-check is to confirm that every sample reached the rarefaction target of 10 000 reads. If one library is badly under-sequenced it can skew diversity estimates or even drop out of analyses entirely. Because one read occupies four lines in a FASTQ file, total-lines ÷ 4 gives the read count.

```bash
# count reads in one compressed file
zcat GC1GC1_R1.fastq.gz | wc -l | awk '{print $1/4 " reads"}'

# loop across every forward read file and list counts
for f in *_R1.fastq.gz; do
  echo -n "$f  "
  zcat "$f" | wc -l | awk '{printf "%d reads\n",$1/4}'
done
```

## 3 · Checking read quality with Q-scores
Sequencers give every base a Phred quality score (Q) that converts machine‐signal strength into an error estimate. In Illumina-style sequencing, the machine builds each read one base at a time and each round of chemical incorporation is called a cycle. Looking at per-cycle numbers lets you spot where quality begins to drop off toward the end of the reads, and lets you set the number you wish to truncate the reads at:

$Q = -10 \log_{10} P(\text{error})$

| Score |	Error rate | Meaning |
| --- | --- | --- |
| Q20	| 1 error in 100 bases (99% accurate) |	Acceptable but starting to get noisy |
| Q30	| 1 error in 1000 bases (99.9% accurate)	| Clean data - the Illumina gold standard |

When a base falls below Q20 it can introduce false sequence variants, so we normally trim reads where Q20/Q30 statistics start to degrade.

The conda package ```seqtk``` provides useful commands for inspecting fastq base quality scores. This command summarises the quality scores of each base:
```bash
seqtk fqchk GC1GC1_R1.fastq.gz
```

You will see output like this:
```
POS  #bases  %A  %C  %G  %T  %N  avgQ  errQ  %low  %high
ALL  956250  26.7 19.3 35.1 18.9 0.0 35.2 19.3 5.5 94.5
1    3825    ...   ...  ...  ... 0.0 31.3 ...  10.3 89.7
⋮
250  3825    ...   ...  ...  ... 0.0 28.4 ...  22.0 74.0
```

For the posiion in ```POS```, ```%low``` indicates the fraction of bases under Q20, while ```%high``` indicates the fraction of bases over Q30. Although somewhat arbitrary, we want to pick a position to truncate at where we maximise accuracy - this could be 235 across all files, for example.







# Understanding QIIME2 files

> QIIME 2 uses two types of files: `.qza` files for storing data (called artifacts) and `.qzv` files for storing visualizations. These file formats are designed to bundle not just the raw or processed data, but also metadata about how the data was generated. This ensures that the full history (or 'provenance') of a dataset is preserved. As a result, any researcher examining your work can see exactly what steps you took and what settings you used. and trace. This helps track every step of your analysis, ensuring transparency and reproducibility, which are essential in science.

QIIME 2 manages your data and results using two main file types: `.qza` (QIIME Zipped Artifact) for data, and `.qzv` for visualizations. These files are more than just containers — they store metadata about how the data was created, and help ensure your workflow is reproducible. For example, if you import a set of DNA sequences and run a quality control step, QIIME 2 automatically logs the commands and parameters used. You can open `.qzv` files in a web browser using [https://view.qiime2.org](https://view.qiime2.org) to explore interactive plots and summaries. This system makes microbiome analysis more transparent and shareable than traditional scripts or spreadsheets.

# Import your paired-end sequences

> DNA sequencing generates data in the form of FASTQ files. Each FASTQ file contains a collection of DNA sequences (called reads) along with quality scores that indicate how confident the sequencer was in calling each base (A, T, C, or G). In paired-end sequencing, two separate FASTQ files are produced for each sample—one for the forward read and one for the reverse read. These reads come from opposite ends of the same DNA fragment and can be merged later to reconstruct the full sequence. Importing these files into QIIME 2 allows the software to begin processing them in a standardized format. Importing means telling QIIME where your raw data is and converting it into its internal format.

For this project the reads were sequences using Illumina paired-end, 250 base pair reads with forward and reverse reads in separate files. The fastq is imported in to a QIIME2 data artifact ending in ```.qza```

```bash
time qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /project/microbiome_workshop/amplicon/data/manifest.csv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
Time to run: 2 minutes

What's this `time` thing? You can add the `time` command to any command line task
 to see how long it took to run.
 {: .notice--info}

Output:
* ```demux.qza```

# Examine the quality of the data

> Before analyzing the sequences, it's important to assess their quality. Sequencing machines aren't flawless—they may misread certain bases, especially toward the ends of reads. A quality check allows us to visualize how reliable each base position is across all reads. If certain positions show low quality, we can trim them off. This ensures that we only use high-confidence sequences in downstream analyses, which leads to more accurate and meaningful results.

We can view the characteristics of the dataset and the quality scores of the data by creating a QIIME2 visualization artifact.

```bash
time qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
 ```
 Time to run: 1 minute

 Output:
 * ```demux.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdemux.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/demux.qzv)

This will create a visualization file. You can download the file to your local computer. From a new terminal window on your local computer copy the file:

```bash
scp <user.name>@login.scinet.science:/path/to/data .
```

Now you can view the file on your local computer using the [QIIME2 visualization server](https://view.qiime2.org).  Alternatively you can view the precomputed file on that server using the button above.

When viewing the data look for the point in the forward and reverse reads where quality scores decline below 25-30. We will need to trim reads to this point to create high quality sequence variants.

# Selecting Sequence Variants

> A sequence variant is a unique DNA sequence found in your dataset. Unlike older methods that grouped similar sequences into Operational Taxonomic Units (OTUs) based on a similarity threshold, modern methods like DADA2 or Deblur aim to resolve individual, error-corrected sequences—called Amplicon Sequence Variants (ASVs). ASVs allow for higher-resolution analysis and more reproducible results. They represent actual biological sequences from the sample, not clusters influenced by arbitrary thresholds. Identifying ASVs is a core step in modern microbiome research. Rather than grouping similar sequences (as older methods did), newer tools try to identify exact sequences to better reflect the real microbes present.


The process of selecting sequence variants is the core processing step in amplicon analysis. This takes the place of "OTU picking" a method of clustering similar data together that was the common method for dealing with sequencing errors until last year.  Three different methods have been published to select sequence variants, [Dada2](https://dx.doi.org/10.1038/nmeth.3869) uses and statistical error correction model, [Deblur](https://dx.doi.org/10.1128/mSystems.00191-16) takes an information theoretic approach and [UNOISE2](https://doi.org/10.1101/081257) applies a heuristic. Each of these methods attempt to remove or correct reads with sequencing errors and then remove chimeric sequences originating from different DNA templates.
For the next step you can select either the Dada2 method or the Deblur method.

**A note on parallel processing**. Both Dada2 and Deblur can independently process each sample. This becomes very important as experiments grow to thousands of samples. In this tutorial we are taking advantage of " course grained parallelism" provided by the multiple cores in our processor. However, Deblur and Dada2 (after the error model learning step) can select sequence variants in a sample totally independently, making the problem "embarrassingly parallel." This means that samples can be split into many different files and processed on arbitrarily many machines simultaneously.
{: .notice--info}

Sequence variant selection is the slowest step in the tutorial. For that reason it is best to submit this step using the SLURM Sbatch scheduler.

## Option 1: Dada2 (Slower)

```bash
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --o-table table-dada2 \
  --o-representative-sequences rep-seqs-dada2 \
  --p-trim-left-f 9 \
  --p-trim-left-r 9 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-n-threads 40 \
  --p-n-reads-learn 200000
```
To submit this command using Sbatch:
```bash
sbatch /project/microbiome_workshop/amplicon/example/qiime2-phosphate-tutorial/dada2.sh
```
Time to run: 35 minutes

Output:
* ```rep-seqs-dada2.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frep-seqs-dada2.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rep-seqs-dada2.qza)
* ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)

## Option 2: Deblur (Faster)
Deblur only uses forward reads at this time. You could get around this by merging your data with an outside tool like [BBmerge](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) then importing your data as single ended. For simplicity, in this tutorial we will just use the forward reads.

```bash
 time qiime deblur denoise-16S \
   --i-demultiplexed-seqs demux.qza \
   --p-trim-length 220 \
   --output-dir deblurresults \
   --p-jobs-to-start 36
```

To submit this command using Sbatch:
```bash
sbatch /project/microbiome_workshop/amplicon/example/qiime2-phosphate-tutorial/deblur.sh
```
Time to run: 4 minutes

Output:
* ```deblurresults/representative_sequences.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Frepresentative_sequences.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/representative_sequences.qza)
* ```deblurresults/stats.qza```
[View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Fstats.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/stats.qza)
* ```deblurresults/table.qza```
[View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Ftable.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/table.qza)

Okay, we have just done the hard part of amplicon sequence analysis.  At this point we have our BIOM count table, the representative sequence variants and a stats file for Deblur.

We have just called sequence variants two different ways. In a real workflow you would only use one method.  From here on out we will use the output of dada2 only: ```table-dada2.qza```.  

# Adding metadata and examining count tables

```bash
time qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt
  ```
  Time to run: 30 seconds

 Output:
 * ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)

# Phylogenetics

> A phylogenetic tree is a diagram that represents the evolutionary relationships among different organisms—or in this case, different microbial DNA sequences. By aligning and comparing the sequences, we can infer how closely related they are and construct a tree that reflects their shared evolutionary history. This tree is essential for certain types of diversity metrics, such as UniFrac, which measure not just the presence or absence of microbes, but how evolutionarily different the communities are. Understanding these relationships helps us interpret microbial functions and ecological roles more effectively. It's useful for certain types of diversity comparisons between samples.

There are a number of diversity metrics like unifrac distance that require the construction of a phylogenetic tree.

## Multiple sequence alignment
First Mafft is used to align the sequences

```bash
time qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Faligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/aligned-rep-seqs.qza)

## Masking sites
Some sites in the alignment are not phylogenetically informative. These sites are masked.

```bash
time qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```masked-aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fmasked-aligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/masked-aligned-rep-seqs.qza)

## Creating a tree
Fastree is used to generate a phylogenetic tree from the masked alignment.
```bash
time qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
```
Time to run: 1 minute

Output:
* ```unrooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Funrooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/unrooted-tree.qza)

## Midpoint rooting
Fastree creates an unrooted tree. We can root the tree at it's midpoint with this command:
```bash
time qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
 ```
Time to run: 5 seconds

Output:
* ```rooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rooted-tree.qza)

# Taxonomic analysis
> **Why taxonomy?**
> Just knowing you have different sequences isn't enough — you want to know what organisms they come from. Taxonomy assigns biological names to the sequences.

Sequence variants are of limited usefulness by themselves. Often we are interested in what kinds of organisms are present in our sample, not just the diversity of the sample.  To identify these sequence variants two things are needed:  a reference database and an algorithm for identifying the sequence using the database.

The primary databases are:

Database | Description | License
---------|-------------|--------
[Greengenes](http://greengenes.secondgenome.com/) | A curated database of archaea and bacteria - static since 2013 | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[Silva](https://www.arb-silva.de/) | The most up-to-date and extensive  database of prokaryotes and eukaryotes, several versions | [Free academic](https://www.arb-silva.de/silva-license-information) / [Paid commercial license](http://www.ribocon.com/silva_licenses)
[The RDP database](https://rdp.cme.msu.edu/) | A large collection of archaeal bacterial and fungal sequences | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[UNITE](https://unite.ut.ee/) | The primary database for fungal ITS and 28S data | Not stated

There are several methods of taxonomic classification available. The most commonly used classifier is the [RDP classifier](https://rdp.cme.msu.edu/classifier/classifier.jsp). Other software includes [SINTAX](http://www.drive5.com/usearch/manual/cmd_sintax.html) and [16S classifier](http://metabiosys.iiserb.ac.in/16Sclassifier/). We will be using the QIIME2's built-in naive Bayesian classifier (which is built on Scikit-learn but similar to RDP), noting that the method, while fast and powerful, has a tendency  [over-classify](http://www.drive5.com/usearch/manual/tax_err.html) reads.

There are two steps to taxonomic classification: [training the classifier](https://docs.qiime2.org/2024.2/tutorials/feature-classifier/) (or using a [pre-trained](https://docs.qiime2.org/2024.2/data-resources/) dataset) and classifying the sequence variants.  Generally it is best to train the classifier on the exact region of the 16S, 18S or ITS you sequenced.
For this tutorial we will be using a classifier model trained on the Silva 99% database trimmed to the V4 region.

```bash
time qiime feature-classifier classify-sklearn \
  --i-classifier  /project/microbiome_workshop/amplicon/data/taxonomy/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza
```
Time to run: 4 minutes

Output:
* ```taxonomy.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qza)

```bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  ```
Time to run: 1 second

* ```taxonomy.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qzv)


Create a bar plot visualization of the taxonomy data:
```bash
qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt \
  --o-visualization taxa-bar-plots.qzv
```
Time to run: 1 minute

Output:
* ```taxa-bar-plots.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxa-bar-plots.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxa-bar-plots.qzv)


Looking at the the ```taxonomy.qzv``` file using https://view/qiime2.org We can see the data presented at different taxonomic levels and grouped by different experimental factors. If we drill down to taxonomic level 5 something looks a bit odd. There's lots of "Rickettsiales;f__mitochondria".  This is really  plant mitochondrial contamination. Some of these samples also have chloroplast contamination.  This kind of Taxonomic filtering isn't available in QIIME2 yet but it can be be done manually.

## Filtering contaminants
By viewing the ```taxonomy.qzv``` file in the browser we can easily search for the sequence ID's that do or do not match "mitochondria or chloroplast". We can also do this via the command line:

First create a list of all sequence variant ids that are not mitochondria or chloroplasts.

```bash
# Export taxonomy data to tabular format
qiime tools export --output-dir taxonomy-export taxonomy.qza

# search for matching lines with grep then select the id column
grep -v -i "mitochondia|chloroplast|Feature" taxonomy-export/taxonomy.tsv | cut  -f 1 > no-chloro-mito-ids.txt
```

Convert our Qiime data artifact to the underlying [biom](http://biom-format.org/) file
```bash
# Export data to biom format
qiime tools export --output-dir dada2-table-export table-dada2.qza
# Move into the directory
cd dada2-table-export

# Convert the HDF5 biom file to a tsv biom file
biom subset-table \
  --input-hdf5-fp feature-table.biom \
  --axis observation \
  --ids ../no-chloro-mito-ids.txt \
  --output-fp feature-table-subset.biom

# Create a new QIIME2 data artifact with the filtered Biom file
qiime tools import \
  --input-path feature-table-filtered.biom \
  --output-path ../table-dada2-filtered.qza \
  --type FeatureTable[Frequency]

cd ..
```
Time to run: 2 minutes

Output:
* ```table-dada2-filtered.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2-filtered.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2-filtered.qza)


Since we have altered the qza file we can create a new bar plots:

```bash
qiime taxa barplot \
  --i-table table-dada2-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt \
  --o-visualization taxa-bar-plots-filtered.qzv

```
Time to run: 1 minute

Output:
* ```taxa-bar-plots-filtered.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxa-bar-plots-filtered.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxa-bar-plots-filtered.qzv)



# Ecological Analysis of Microbial Communities in Agriculture

Microbial ecology is the study of how microbes interact with each other and their environment. In agricultural settings, these interactions can profoundly affect crop health, nutrient cycling, and disease suppression. Ecological analysis allows us to characterize microbial diversity, compare communities between samples (like different soil types or plant treatments), and identify patterns that inform how microbial communities function and respond to changes.

### Alpha Diversity: Within-Sample Diversity

Alpha diversity refers to the number and evenness of microbial species in a single sample. Common metrics include:

- **Observed Features**: The count of unique sequence variants or species.
- **Shannon Diversity Index**: Takes into account both richness and evenness; higher values indicate more diverse communities.
- **Faith’s Phylogenetic Diversity**: Considers how evolutionarily diverse a community is, using a phylogenetic tree.

In agriculture, alpha diversity can reflect the health or disturbance of a soil environment — for example, more diverse microbial communities are often associated with healthier or less disturbed soils.

### Beta Diversity: Between-Sample Differences

Beta diversity compares microbial composition across samples. It helps you answer questions like “Do different fertilizers cause different microbial communities?” or “Are microbes in root zones different from those in bulk soil?”

Common metrics include:

- **Bray-Curtis Dissimilarity**: Based on differences in counts.
- **Jaccard Distance**: Based on shared presence/absence of species.
- **UniFrac (weighted/unweighted)**: Measures how phylogenetically different two communities are.

QIIME 2 calculates these distances and uses techniques like Principal Coordinates Analysis (PCoA) to visualize differences.

### Ordination and Visualization

Ordination is a statistical method that reduces complex distance matrices into 2 or 3 dimensions for easier visualization. QIIME 2 uses **Emperor** to create 3D plots of PCoA results. These plots help identify clusters of similar samples, which can reflect treatment effects, environmental conditions, or sample types.

### Differential Abundance Analysis

Identifying which microbes are more or less abundant between groups (e.g., treated vs. untreated plants) is crucial. While traditional statistical methods like t-tests aren't suitable for microbiome data (because it's compositional and sparse), specialized tools like **Songbird**, **ANCOM**, or **DESeq2** can model these differences.

As an example, if a plant genotype promotes beneficial microbes that suppress disease, a differential abundance test might identify those beneficial taxa as more abundant in samples from that genotype.

### Why It Matters for Agriculture

Microbiome analysis in agriculture helps us:

- Understand how farming practices impact soil health.
- Improve crop yields by promoting beneficial microbes.
- Monitor soil and plant health non-invasively.
- Design more sustainable systems that rely less on chemical inputs.

By linking microbial community structure with environmental or treatment variables, ecological analysis becomes a powerful tool for both research and practical decision-making in agriculture.


> **What is ordination?**
> It's a way to simplify and visualize complex data. For example, it helps you see patterns in how samples cluster based on their microbial composition.

Ordination is a dimensionality reduction technique that enables the visualization of sample differences. QIIME has a plugin called emperor that calculates a Bray-Curtis dissimilarity matrix and uses principal coordinates analysis (PCoA). you could also export the pcoa data and plot it yourself in the package of your choice.

 ```bash
 qiime emperor plot --i-pcoa core-diversity/bray_curtis_pcoa_results.qza --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt --o-visualization pcoa-visualization.qzv
 ```
 Time to run: 15 seconds

 Output:
 * ```pcoa-visualization.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fpcoa-visualization.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/pcoa-visualization.qzv)


# Differential abundance of sequence variants

> ⚠️ **Note:** The `gneiss` plugin has been deprecated in recent QIIME 2 versions. Consider using `songbird` or external tools like `ALDEx2` for differential abundance analysis.

If you are doing an experimental manipulation rather than just observing an environment you will likely have an experimental design with treatments and want to know which bacteria respond to these treatments. The best way to go this is an active area of applied statistics research.

The problem is challenging for several reasons:
* The data is compositional so the abundance of each taxa affects every other taxa.
* The data is over-dispersed count data, fitting (arguably) a negative binomial model.
* The data is sparse and some 0's mean a taxa is not present while other zeros mean an organism is present at a level below the limit of detection for the sequences sampled.
* Each library is sampled to a different depth so the issue of how to standardize the data comes up (simply dividing by the sum does not work).

**To rarefy or not to rarefy?** It is a common but controversial practice to downsample count data to the lowest count in your dataset to get around the issue of differential sequencing depth. In their paper titled "Waste not want not, why rarifying microbiome data is inadmissible" [McMurdie et al. (2014)](https://doi.org/10.1371/journal.pcbi.1003531) point out that this is a large waste of data and statistical power, and advocate for using differential expression software like DESeq2 that uses special normalizations and a negative binomial distribution to model data. The software uses a generalized linear model so it has a very flexible experimental design interface.  [Weiss et al. (2017)](https://doi.org/10.1186/s40168-017-0237-y) argue that the assumptions underling both the normalization and the distribution used by DEseq2 and other normalization methods are inappropriate for microbiome data. [Ancom](https://dx.doi.org/10.3402%2Fmehd.v26.27663) uses a zero inflated Gaussian model but only allows for simple two-way comparisons not richer statistical models. Gneiss [(paper)](http://doi.org/10.1128/mSystems.00162-16), [(tutorial)](https://forum.qiime2.org/t/gneiss-tutorial/932) is currently the only compositional method available in QIIME2 (support for ANCOM was dropped).  The only thing everyone agrees on is that agrees you can't just do a T-test.
{: .notice--warning}



## Gneiss analysis (deprecated)
Gneiss applies a method for compositional data borrowed from geology to "sidestep" the question of the absolute changes of sequences and "instead look at the balance between particular subsets of the microbial community," [(Morton et al. 2017)](https://doi.org/10.1128/mSystems.00162-16). For more on the concept of balances see this [post](https://github.com/biocore/gneiss/blob/master/ipynb/balance_trees.ipynb). Sequence variants are hierarchically clustered based on environmental gradients or co-occurrence. Then the isometric log ratios are calculated and compared for subsets of taxa. While individual taxa cannot be compared, groups of taxa responding to environmental effects can be compared. Statistical tests of for differences can also be applied.

Since we don't have a particular environmental gradient that is structuring lets start by doing a hierarchical clustering based on co-occurrence patterns.

First add pseudocounts to each cell in the matrix. This is done so that log transformations can be taken across the table.

```bash
time qiime gneiss add-pseudocount \
    --i-table table-dada2-filtered.qza \
    --p-pseudocount 1 \
    --o-composition-table composition.qza
```
Time to run: 6 seconds

Output:
* ```composition.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcomposition.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/composition.qza)


Perform [Ward's agglomerative clustering](https://arxiv.org/abs/1111.6285)
```bash    
time qiime gneiss correlation-clustering \
    --i-table composition.qza\
    --o-clustering hierarchy.qza
```
Time to run: 5 minutes

Output:
* ```hierarchy.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fhierarchy.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/hierarchy.qza)

A tree has now been generated that can be used for making comparisons of sample groups.

Calculate the isometric log transforms on each internal node of the tree

```bash
time qiime gneiss ilr-transform \
    --i-table composition.qza \
    --i-tree hierarchy.qza \
    --o-balances balances.qza
```
Time to run: 15 seconds

Output:
* ```balances.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fbalances.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/balances.qza)

The balances are normally distributed an can now be analyzed using mixed linear
models  We can perform a regression on the three categorical data types, Genotype, Fraction (soil or endophytic compartment) or Soil).  Themodel explains about 10% of the total variation at all nodes of the trees. This is typical for these complex experiments.  The amount that can be explained increases as we move up the covariance tree. Overall the most predictive factor is Genotype which is encouraging.

```bash
time qiime gneiss ols-regression \
    --p-formula "Genotype+Soil+Fraction" \
    --i-table balances.qza \
    --i-tree hierarchy.qza \
    --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt \
    --o-visualization regression_summary.qzv
```

Output:
* ```regression_summary.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fregression_summary.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/regression_summary.qzv)


One of the assumptions if the ordinary least squares model is that the fixed factors are random, in other words the authors randomly arrived at the genotypes they knocked out. Of course that's not true, the genotypes were selected because they had an impact on the phosphorus stress response.  They are Fixed factors. A mixed linear model can account for fixed and random factors and effects. Gneiss offerers a linear mixed model regression too but the interface seems to be in development  so there is not much I can say about it but we can try it now. Statistical modeling is done by the [statsmodels](http://www.statsmodels.org/stable/mixed_linear.html) python package.

```bash
qiime gneiss lme-regression \
  --p-formula "Genotype" \
  --i-table balances.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt \
  --p-groups Soil \
  --o-visualization linear_mixed_effects_model.qzv \
```

Output:
* ```linear_mixed_effects_model.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Flinear_mixed_effects_model.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/linear_mixed_effects_model.qzv)


We can look at the most statistically significant balances and examine what taxa make up those partitions.

```bash
qiime gneiss balance-taxonomy \
    --i-balances balances.qza \
    --i-tree hierarchy.qza \
    --i-taxonomy taxonomy.qza \
    --p-taxa-level 2 \
    --p-balance-name 'y0' \
    --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt  \
    --m-metadata-category Genotype \
    --o-visualization y0_taxa_summary.qzv
```


Output:
* ```y0_taxa_summary.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fy0_taxa_summary.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/y0_taxa_summary.qzv)


In this case the y0 balance is a split between samples that have plants in them and raw soil. It makes sense that this is the largest effect.  What happens if you run balance y2 or decrease the taxonomic level?

# Other analyses

This tutorial covered a range of analyses that can be done with microbiome data but there are other types on analyses that can be done too.

* Functional analysis - Several packages attempt to impute function from taxonomy including [PiCrust](https://picrust.github.io/picrust/), [Tax4fun](http://tax4fun.gobics.de/), [Piphillin](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166104)
* Inferring ecological interaction networks -[SPIEC-EASI](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). Co-Variance is an issue, but SPIEC-EASI attempts to model conditional independence.
* Data management tools - [Qiita](https://qiita.ucsd.edu/), and from ARS scientist Dan Manter [myPhyloDB](http://www.myphylodb.org/)
* Set analysis - From ARS Scientist Devin Coleman-Derr [MetaComet](https://probes.pw.usda.gov/MetaCoMET/)
* Other general analysis tools - [Mothur](https://www.mothur.org/) and the R-based [Phyloseq](https://joey711.github.io/phyloseq/)
