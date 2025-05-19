---
title: "Amplicon analysis with QIIME2"
excerpt: "An example workflow using QIIME2 version 2024.2"
layout: single
mathjax: true
author: "George Kalogiannis"

---
By George Kalogiannis, Designed from the official [QIIME2 tutorials](https://docs.qiime2.org/2024.2/tutorials/)
{% include toc %}

# Learning Outcomes
In this workshop, you’ll learn how to process and analyse microbial amplicon sequencing data using QIIME 2. Starting from raw FASTQ files, you’ll import your data, assess its quality, and use tools like DADA2 or Deblur to generate high-confidence sequence variants (ASVs). You’ll also learn how to assign taxonomy to those sequences and create summaries of microbial community structure.

By the end of the session, you will be able to:
- Import and quality-check amplicon data in QIIME 2
- Identify ASVs and generate a feature table
- Assign taxonomy using a reference database
- View and interpret results through interactive visualisations
- Submit resource-intensive steps to an HPC cluster (CHPC) for efficient processing

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

Modern DNA sequencers, Illumina (short, very accurate reads), PacBio HiFi and Oxford Nanopore (long, increasingly accurate reads), do not give you a “genome” or a list of species. They give you a bag of millions of little DNA fragments called reads. 

Reads can be created in two ways: single or paired. In single-end sequencing, the machine reads each DNA fragment from one end only, starting at the 5′ end and working its way across the fragment until the read is complete. In paired-end sequencing, the machine reads from both ends of the same fragment: one read starts at the 5′ end and moves toward the middle (the forward read), and the other starts at the opposite 3′ end and also moves toward the middle (the reverse read). This gives you two reads per fragment that come from opposite directions, which helps improve accuracy and recover more of the original sequence.

![Single vs Paired-End Sequencing](../../assets/images/single_vs_paired.png)

These reads are written to disk as plain-text files, almost always compressed and ending in .fastq.gz.

- FASTA files ( .fasta, .fa ) hold only the nucleotide strings; they are mainly used once reads have been cleaned and merged.
- FASTQ files ( .fastq, .fq ) hold both the sequence and its Phred quality scores, telling you how confident the instrument was at each base.

>Once you have the downloads in place, peek inside one of the real files, say GC1GC1_R1.fastq.gz (forward reads from replicate GC1, experiment GC1), to convince yourself what a FASTQ record looks like and to check that the quality scores are sensible before you hand the data to QIIME 2. You can do this in bash like this:
>```bash
># ── Bash one-liner: show the first two reads (8 lines) ───────────
>zcat GC1GC1_R1.fastq.gz | head -n 8
>```
>Breaking that command down, ```zcat``` prints the unzipped file, ```head -n 8``` takes the top 8 lines of the file.


## 2 · “How many reads did I get?”

Before you begin any analysis, it’s important to check that your files actually contain the number of reads you expect. This step can catch problems like incomplete downloads, interrupted sequencing runs, or incorrect file handling. 

>Since each read in a FASTQ file takes up four lines (a header, sequence, separator, and quality line), you can calculate the number of reads by dividing the total number of lines by 4. This gives you a quick way to verify that your samples are complete before moving on.
>```bash
>## count reads in one compressed file
>zcat GC1GC1_R1.fastq.gz | wc -l | awk '{print $1/4 " reads"}'
>
># loop across every forward read file and list counts
>for f in *_R1.fastq.gz; do
>  echo -n "$f  "
>  zcat "$f" | wc -l | awk '{printf "%d reads\n",$1/4}'
>done
>```

## 3 · Checking read quality with Q-scores
Sequencers give every base a Phred quality score (Q) that converts machine‐signal strength into an error estimate. In Illumina-style sequencing, the machine builds each read one base at a time and each round of chemical incorporation is called a cycle. Looking at per-cycle numbers lets you spot where quality begins to drop off toward the end of the reads, and lets you set the number you wish to truncate the reads at:

| Score |	Error rate | Meaning |
| --- | --- | --- |
| Q20	| 1 error in 100 bases (99% accurate) |	Acceptable but starting to get noisy |
| Q30	| 1 error in 1000 bases (99.9% accurate)	| Clean data - the Illumina gold standard |

When a base falls below Q20 it can introduce false sequence variants, so we normally trim reads where Q20/Q30 statistics start to degrade.

> The conda package ```seqtk``` provides useful commands for inspecting fastq base quality scores. This command summarises the quality scores of each base:
>```bash
>seqtk fqchk GC1GC1_R1.fastq.gz
>```
>
>You will see output like this:
>```
>POS  #bases  %A  %C  %G  %T  %N  avgQ  errQ  %low  %high
>ALL  956250  26.7 19.3 35.1 18.9 0.0 35.2 19.3 5.5 94.5
>1    3825    ...   ...  ...  ... 0.0 31.3 ...  10.3 89.7
>⋮
>250  3825    ...   ...  ...  ... 0.0 28.4 ...  22.0 74.0
>```
>For the posiion in ```POS```, ```%low``` indicates the fraction of bases under Q20, while ```%high``` indicates the fraction of bases over Q30. Although somewhat arbitrary, we want to pick a position to truncate at where we maximise accuracy - this could be 235 across all files, for example.


## 4 · “Do my paired files really match?”

This section checks that your paired-end FASTQ files truly match, which is essential before you run any microbiome pipeline like QIIME 2 or DADA2. Each pair of files—one forward read file (*_R1.fastq.gz) and one reverse read file (*_R2.fastq.gz) should contain the same number of reads, and those reads should be from the same DNA fragments. If the files get out of sync (due to download errors, interruptions, or corrupted files), your pipeline will fail or produce misleading results.

> This compares the number of reads in each file. ```zcat``` unzips the .fastq.gz file directly in memory. ```wc -l``` counts the total number of lines. FASTQ files use 4 lines per read, so dividing by 4 gives you the read count.
>```bash
>echo "R1: $(($(zcat GC1GC1_R1.fastq.gz | wc -l)/4))   R2: $(($(zcat GC1GC1_R2.fastq.gz | wc -l)/4))"
>```

> This checks that the read IDs match between files. ```sed -n '1~4p'``` pulls out every 4th line starting from line 1 (the read headers). ```cut -d' ' -f1``` trims each header to just the read ID, ignoring barcode info. ```paste``` prints them side by side so you can compare them easily.
>```bash
>paste \
>  <(zcat GC1GC1_R1.fastq.gz | sed -n '1~4p' | cut -d' ' -f1 | head -5) \
>  <(zcat GC1GC1_R2.fastq.gz | sed -n '1~4p' | cut -d' ' -f1 | head -5)
>```

## 5 · Where is genetic data stored?
*needs text here*

## 4 · BLAST searches

If you’ve ever used BLAST (Basic Local Alignment Search Tool), you know the basic idea: you give it a DNA or protein sequence, and it compares that sequence to a reference database to find the best matches. It reports which known sequences are most similar, how long the matching region is, and how confident the match is. This is the foundation of how we assign names or functions to unknown sequences.

In amplicon analysis, tools like QIIME 2 do something very similar. After denoising your reads into amplicon sequence variants (ASVs), which are essentially cleaned-up, exact sequences, it matches each ASV to a reference database such as SILVA, Greengenes, or GTDB. Each ASV is compared to known 16S rRNA gene sequences in the database to find its closest match, which is then used to assign it a taxonomic label (e.g. Bacillus subtilis, Pseudomonas, etc.).

So while QIIME doesn’t run BLAST directly by default, it’s doing the same type of work: comparing your sequences to a trusted reference and reporting the best biological match. It’s just automated, scaled up, and optimised for microbial marker genes like 16S.

> Lets run a blast search for the first sequence in our file. Print the first sequence and its information and paste it into the BLAST Nucleotide website: [blast.ncbi.nlm.nih.gov/](blast.ncbi.nlm.nih.gov/).
>```bash
>zcat GC1GC1_R1.fastq.gz | head -n 4
>```

 **What do you see? Is it a bacterium that you have identified?** Hold onto this thought until we discuss contamination later.

# Intro to HPC

High-Performance Computing (HPC) refers to the use of supercomputers or parallel computing techniques to perform complex computations quickly. HPC systems aggregate computing resources to solve problems that would be unfeasible or time-consuming using conventional methods. These systems play a critical role in research fields like bioinformatics, physics, climate modeling, and engineering.

But how do HPCs achieve better results in less time? This is because of parallel computing, which is the ability to divide large tasks into smaller pieces that can be completed simultaneously. You can think of this like a group of people working together to build a house - each person has a specific task which they work on at the same time, finishing the job faster.

A frequent misunderstanding is that executing code on a cluster will automatically result in improved performance. Clusters do not inherently accelerate code execution; for enhanced performance, the code must be explicitly adapted for parallel processing, a task that falls under your responsibility as a programmer.

HPC's typically have thousands of cores, made up by separate computers (nodes) that are linked and can transfer data to each other. The Lengau CHPC we will be using today has a total of 33,436 cores available to its users.

![Sequencial vs Parallel Computing](../../assets/images/parallel.png)

## Qiime2 on HPC
- CPU Power: Many QIIME2 plugins, such as dada2, vsearch, and feature-classifier, can utilize multithreading to significantly reduce processing time. On an HPC, users can request dozens of CPU cores per job, greatly accelerating large analyses involving hundreds or thousands of samples.

- RAM Availability: Some QIIME2 steps, especially taxonomy classification and large distance matrix calculations, can consume large amounts of memory (10–100+ GB). HPC systems often provide nodes with hundreds of gigabytes of RAM, preventing crashes and memory overflows that might occur on a personal computer.

- Storage Capacity: Intermediate QIIME2 artifacts and results (e.g., .qza, .qzv files) accumulate quickly and may exceed local storage limits. HPC environments usually offer shared high-capacity storage systems (terabytes to petabytes) optimized for fast read/write speeds, which helps manage these large datasets efficiently.

- Job Scheduling: HPCs use job scheduling systems (e.g., SLURM, PBS) to manage and queue computational tasks. This enables long or resource-heavy QIIME2 workflows to run in the background without tying up your local machine.

## Connecting to the CHPC
*insert text here*

## Interacting with the CHPC
*sftp/scp/ssh*

*conda envs*

*directories*

# QIIME2

## Understanding QIIME2 Files
QIIME 2 organises your analysis using two main file types: .qza files for data (called QIIME Zipped Artifacts) and .qzv files for visualisations (QIIME Zipped Visualisations). These files are more than simple containers, as they include detailed information about how each was generated, such as the input files, commands, and parameters used. This is known as provenance tracking, and it allows others to trace every step of your workflow. You can open .qzv files in your browser at [https://view.qiime2.org](https://view.qiime2.org) to explore interactive summaries like quality plots, taxonomic profiles, and diversity results. This approach helps make your microbiome analysis more transparent, easier to share, and reproducible from start to finish.

## Import your paired-end sequences
DNA sequencing generates data in the form of FASTQ files. Importing these files into QIIME 2 allows the software to begin processing them in a standardized format.

QIIME2 expects file locations to be written in a manifest file (here ```wednesday_manifest.tsv```), whose internal format is specified in the ```--input-format```. There is the option of ```PairedEndFastqManifestPhred33V2``` or ```SingleEndFastqManifestPhred33V2```, depending on if we supply only forward reads or paired forward-reverse reads. We would format our manifest file accordingly:

The fastq.gz absolute filepaths may contain environment variables (e.g., $HOME or $PWD). The following example illustrates a simple fastq manifest file for paired-end read data for four samples.
```
sample-id     forward-absolute-filepath       reverse-absolute-filepath
sample-1      $PWD/some/filepath/sample0_R1.fastq.gz  $PWD/some/filepath/sample1_R2.fastq.gz
sample-2      $PWD/some/filepath/sample2_R1.fastq.gz  $PWD/some/filepath/sample2_R2.fastq.gz
sample-3      $PWD/some/filepath/sample3_R1.fastq.gz  $PWD/some/filepath/sample3_R2.fastq.gz
sample-4      $PWD/some/filepath/sample4_R1.fastq.gz  $PWD/some/filepath/sample4_R2.fastq.gz
```

The following example illustrates a simple fastq manifest file for fastq single-end read data for two samples.
```
sample-id     absolute-filepath
sample-1      $PWD/some/filepath/sample1_R1.fastq
sample-2      $PWD/some/filepath/sample2_R1.fastq
```

>Importing means telling QIIME where your raw data is and converting it into its internal format (when we run a command with time in front of it, it shows us how long it took to run):
>```bash
>time qiime tools import \
>  --type 'SampleData[PairedEndSequencesWithQuality]' \
>  --input-path wednesday_manifest.tsv \
>  --output-path wednesday_outputs/demux.qza \
>  --input-format PairedEndFastqManifestPhred33V2
>```
>
>Time to run: 2 minutes
>
>Output: ```demux.qza```

## Examine the quality of the data
Before analyzing the sequences, it's important to assess their quality like we did before in the [Genetics Primer](#3--checking-read-quality-with-q-scores). A quality check using QIIME2 allows us to visualize how reliable each base position is across all reads. If certain positions show low quality, we can trim them off. 

>We can view the characteristics of the dataset and the quality scores of the data by creating a QIIME2 visualization artifact.
>
>```bash
>time qiime demux summarize \
>  --i-data wednesday_outputs/demux.qza \
>  --o-visualization wednesday_outputs/demux.qzv
> ```
>
> Time to run: 1 minute
>
>Output: * ```demux.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdemux.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/demux.qzv)

>This will create a visualization file. You can download the file to your local computer. From a new terminal window on your local computer copy the file:
>
>```bash
>scp <user.name>@lengau.chpc.ac.za:/path/to/data .
>```
>Now you can view the file on your local computer using the [QIIME2 visualization server](https://view.qiime2.org). Alternatively you can view the precomputed file on that server using the button above.
>
>When viewing the data look for the point in the forward and reverse reads where quality scores decline below 25-30. We will need to trim reads to this point to create high quality sequence variants.

## Selecting Sequence Variants

Sequencing isn’t perfect and errors, noise, and artefacts are common. To make sense of the data, the first core step is to group similar reads together and identify the true biological sequences that were present in the sample. This process is called **sequence variant picking**.

Modern tools like [Dada2](https://dx.doi.org/10.1038/nmeth.3869), [Deblur](https://dx.doi.org/10.1128/mSystems.00191-16), and [UNOISE2](https://doi.org/10.1101/081257) aim to identify these exact sequences by using statistical error correction models, taking an information theoretic approach and applying a heuristic. This is a major improvement over older methods like OTU picking, where sequences were grouped based on a similarity threshold. ASVs (Amplicon Sequence Variants) provide higher resolution, greater reproducibility, and better biological accuracy.

For the next step you can select either the Dada2 method or the Deblur method. Sequence variant selection is the slowest step in the tutorial.

### Option 1: Dada2 (Slower)

>```bash
>time qiime dada2 denoise-paired \
>  --i-demultiplexed-seqs demux.qza \   # the imported FASTQ data (paired-end reads).
>  --o-table table-dada2 \   # the output file containing the ASV count table.
>  --o-representative-sequences rep-seqs-dada2 \   # the actual DNA sequences of the ASVs.
>  --p-trim-left-f 9 \   # trims 9 bases from the start of each forward/reverse read (e.g. to remove primers).
>  --p-trim-left-r 9 \
>  --p-trunc-len-f 220 \   # truncates reads to 220/200 bases (based on where quality drops off).
>  --p-trunc-len-r 200 \
>  --p-n-threads 24 \   # number of CPU threads to use. Adjust based on your system.
>  --p-n-reads-learn 200000    # number of reads used to learn the error model. You can lower this on small datasets.
>```
>
>Time to run: 35 minutes
>
>Output:
> - ```rep-seqs-dada2.qza```: These are your representative sequences: the exact, cleaned-up ASV sequences found in your dataset. They are matched against reference databases later (like SILVA or Greengenes) to assign taxonomy. [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frep-seqs-dada2.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rep-seqs-dada2.qza)
> - ```table-dada2.qzv```: This is your feature table, stored in QIIME 2’s .qza format. Internally, it follows the BIOM (Biological Observation Matrix) standard. It records how many times each ASV appears in each sample, similar to a species count table. This is the key file you’ll use for diversity analysis, statistical testing, and visualisation. [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)

### Option 2: Deblur (Faster)
> Deblur only uses forward reads at this time. You could get around this by merging your data with an outside tool like [BBmerge](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) then importing your data as single ended. For simplicity, in this tutorial we will just use the forward reads.
> 
>```bash
> time qiime deblur denoise-16S \
>   --i-demultiplexed-seqs demux.qza \   # the imported FASTQ data (paired-end reads).
>   --p-trim-length 220 \    # truncates reads to 220 bases (based on where quality drops off).
>   --output-dir deblurresults \   # the output directory containing the output files.
>   --p-jobs-to-start 24    # number of CPU threads to use. Adjust based on your system.
>```
>
> Time to run: 4 minutes
>
> Output: Okay, we have just done the hard part of amplicon sequence analysis.  At this point we have our BIOM count table, the representative sequence variants and a stats file for Deblur.
>
> - ```deblurresults/representative_sequences.qza```: A list of the representative sequences (representative_sequences.qza), which are the actual DNA sequences of the ASVs. [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Frepresentative_sequences.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/representative_sequences.qza)
> - ```deblurresults/stats.qza```: A statistics file which logs how many reads passed quality control and how many were removed at each filtering step. [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Fstats.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/stats.qza)
> - ```deblurresults/table.qza```: The same type of table as the one outputted by Dada2 above.
[View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Ftable.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/table.qza)


We have just called sequence variants two different ways. In a real workflow you would only use one method.  From here on out we will use the output of dada2 only: ```table-dada2.qza```.  

## Taxonomic analysis
Sequence variants (ASVs) are high-resolution markers of microbial diversity, but on their own, they don’t tell us which organisms are present. While they help define the structure of the community, we often want to go a step further and identify the actual microbes behind each sequence. To do that, we need to assign taxonomy: linking each ASV to a known group of bacteria or archaea.

This process requires two key components: a reference database of well-annotated 16S rRNA sequences, and a classification algorithm that can compare our ASVs to that database. QIIME 2 handles this automatically using pretrained classifiers. The classifier searches for the closest match in the reference and labels each ASV accordingly, down to the genus or species level if possible.

Earlier, when you ran a BLAST search on one of your raw reads, you may have seen a match to Arabidopsis thaliana. This illustrates why taxonomy assignment is more than just a formality. Host DNA, especially from chloroplasts or mitochondria, can end up in your sequences — and without a curated, microbial-only reference, those contaminants might remain in your results. By using the right database, we ensure that only relevant microbial taxa are kept for downstream analysis.
The primary databases are:

Database | Description | License
---------|-------------|--------
[Greengenes](http://greengenes.secondgenome.com/) | A curated database of archaea and bacteria - static since 2013 | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[Silva](https://www.arb-silva.de/) | The most up-to-date and extensive  database of prokaryotes and eukaryotes, several versions | [Free academic](https://www.arb-silva.de/silva-license-information) / [Paid commercial license](http://www.ribocon.com/silva_licenses)
[The RDP database](https://rdp.cme.msu.edu/) | A large collection of archaeal bacterial and fungal sequences | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[UNITE](https://unite.ut.ee/) | The primary database for fungal ITS and 28S data | Not stated

There are several methods of taxonomic classification available. The most commonly used classifier is the [RDP classifier](https://rdp.cme.msu.edu/classifier/classifier.jsp). Other software includes [SINTAX](http://www.drive5.com/usearch/manual/cmd_sintax.html) and [16S classifier](http://metabiosys.iiserb.ac.in/16Sclassifier/). We will be using the QIIME2's built-in naive Bayesian classifier (which is built on Scikit-learn but similar to RDP), noting that the method, while fast and powerful, has a tendency  [over-classify](http://www.drive5.com/usearch/manual/tax_err.html) reads.

There are two steps to taxonomic classification: [training the classifier](https://docs.qiime2.org/2024.2/tutorials/feature-classifier/) (or using a [pre-trained](https://docs.qiime2.org/2024.2/data-resources/) dataset) and classifying the sequence variants.  Generally it is best to train the classifier on the exact region of the 16S, 18S or ITS you sequenced.

>For this tutorial we will be using a classifier model trained on the Silva 99% database trimmed to the V4 region.
>
>```bash
>time qiime feature-classifier classify-sklearn \
>  --i-classifier  /project/microbiome_workshop/amplicon/data/taxonomy/gg-13-8-99-515-806-nb-classifier.qza \
>  --i-reads rep-seqs-dada2.qza \
>  --o-classification taxonomy.qza
>```
>
>Time to run: 4 minutes
>
>Output: ```taxonomy.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qza)


>Create a bar plot visualization of the taxonomy data:
>```bash
>qiime taxa barplot \
>  --i-table table-dada2.qza \
>  --i-taxonomy taxonomy.qza \
>  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt \
>  --o-visualization taxa-bar-plots.qzv
>```
>Time to run: 1 minute
>
>Output: ```taxa-bar-plots.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxa-bar-plots.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxa-bar-plots.qzv)

