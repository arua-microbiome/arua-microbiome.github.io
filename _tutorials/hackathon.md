---
title: "Friday Hackathon"
excerpt: "An Amplicon Sequencing Hackathon (in Groups)"
layout: single
mathjax: true
authors: "George Kalogiannis & Balig Panossian"

---

By George Kalogiannis & Balig Panossian, Designed from the official [QIIME2 tutorials](https://docs.qiime2.org/2024.2/tutorials/)
{% include toc %}


## Scenario  
Teams receive paired-end 16S rRNA FASTQ files from faecal samples of **three human communities** that differ in agricultural practice life styles (hunter-gatherer, agriculture/hunter-gatherer, agriculture). Companion metadata and manifest sheets are provided.  

Your mission, **by 12:30 on Friday**, is to transform these raw reads into an ecological narrative and defend it in a **four-slide deck** (cover + three content slides, 5 minutes total). You must work around the scheduled times in the agenda.

**Optional personal dataset:** If your group already has a 16S amplicon dataset that you’d like to showcase, you may substitute it for the workshop dataset—provided it fits the same format we have illustrated (paired-end FASTQ files plus a valid sample-metadata TSV). The grading rubric and four-slide limit remain unchanged, so choose wisely: your own data must still let you demonstrate the full pipeline and draw clear ecological conclusions within the allotted time.


## What We’ll Be Grading (implicitly)  

- **Judgement, not button-mashing**  
  *Why those trim points? Why that diversity metric?*  

- **Story compression**  
  *Can you merge QC, processing, and results into three coherent, uncluttered slides?*  

- **Biological connection**  
  *Do your stats and graphics convincingly relate microbiome shifts to agricultural practice?*  

- **Provenance awareness**  
  *Reference file names/paths or show the provenance stamp to prove reproducibility.*  

## The Files

All files can be found at ```/mnt/lustre/groups/WCHPC/friday_data``` on the CHPC.

>**Hints:** Even though the QC plots _may_ look good, denoising with too high of a length can remove samples that don't meet certain length requirements. Think about setting the truncation length to near 150-175. When denoising, use ```100000``` reads to learn the errors from.


