---
title: "Ecological analysis with QIIME2"
excerpt: "An example workflow using QIIME2 version 2024.2"
layout: single
mathjax: true
author: "George Kalogiannis"

---
You’ve reached the point in the workflow where you can confidently download raw amplicon reads, clean them up, and generate reliable ASV tables and taxonomic summaries. In other words, you now know what sequences are present and how many of each occur in every sample. The next step is to turn those counts into ecological insight: measuring within-sample diversity, comparing whole communities across treatments, and linking patterns to phylogenetic relationships. The sections that follow introduce the core ecological tools you’ll use—alpha and beta diversity metrics, phylogeny-based distances, and statistical tests that reveal how environmental or experimental factors shape the microbiome.
{% include toc %}

# Learning Outcomes
By the end of this session you will be able to:
- Compute alpha- and beta-diversity metrics and visualise them in interactive .qzv files
- Construct a phylogenetic tree and run UniFrac-based ordination and PERMANOVA tests
- Export QIIME artifacts to R with qiime2R and produce custom ggplot figures

## Adding metadata and examining count tables
Once you have generated your ASV count table, the next step is to summarise it and explore how many reads were retained in each sample. This command uses qiime feature-table summarize to produce an interactive .qzv file that includes per-sample read counts, total feature counts, and an overview of how balanced your dataset is.

>```bash
>time qiime feature-table summarize \
>  --i-table table-dada2.qza \
>  --o-visualization table-dada2.qzv \
>  --m-sample-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt
>  ```
>  Time to run: 30 seconds
>
> Output: ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)


## Linking Taxa to Metadata

```bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  ```
Time to run: 1 second

* ```taxonomy.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qzv)

## Filtering contaminants
Looking at the the ```taxonomy.qzv``` file using https://view/qiime2.org We can see the data presented at different taxonomic levels and grouped by different experimental factors. If we drill down to taxonomic level 5 something looks a bit odd. There's lots of "Rickettsiales;f__mitochondria".  This is really  plant mitochondrial contamination. Some of these samples also have chloroplast contamination.  This kind of Taxonomic filtering isn't available in QIIME2 yet but it can be be done manually.

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




## Phylogenetics

> A phylogenetic tree is a diagram that represents the evolutionary relationships among different organisms—or in this case, different microbial DNA sequences. By aligning and comparing the sequences, we can infer how closely related they are and construct a tree that reflects their shared evolutionary history. This tree is essential for certain types of diversity metrics, such as UniFrac, which measure not just the presence or absence of microbes, but how evolutionarily different the communities are. Understanding these relationships helps us interpret microbial functions and ecological roles more effectively. It's useful for certain types of diversity comparisons between samples.

There are a number of diversity metrics like unifrac distance that require the construction of a phylogenetic tree.

### Multiple sequence alignment
First Mafft is used to align the sequences

```bash
time qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Faligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/aligned-rep-seqs.qza)

### Masking sites
Some sites in the alignment are not phylogenetically informative. These sites are masked.

```bash
time qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```masked-aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fmasked-aligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/masked-aligned-rep-seqs.qza)

### Creating a tree
Fastree is used to generate a phylogenetic tree from the masked alignment.
```bash
time qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
```
Time to run: 1 minute

Output:
* ```unrooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Funrooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/unrooted-tree.qza)

### Midpoint rooting
Fastree creates an unrooted tree. We can root the tree at it's midpoint with this command:
```bash
time qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
 ```
Time to run: 5 seconds

Output:
* ```rooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rooted-tree.qza)


# Ecological Analysis of Microbial Communities in Agriculture

Microbial ecology is the study of how microbes interact with each other and their environment. In agricultural settings, these interactions can profoundly affect crop health, nutrient cycling, and disease suppression. Ecological analysis allows us to characterize microbial diversity, compare communities between samples (like different soil types or plant treatments), and identify patterns that inform how microbial communities function and respond to changes.

## Alpha Diversity: Within-Sample Diversity

Alpha diversity refers to the number and evenness of microbial species in a single sample. Common metrics include:

- **Observed Features**: The count of unique sequence variants or species.
- **Shannon Diversity Index**: Takes into account both richness and evenness; higher values indicate more diverse communities.
- **Faith’s Phylogenetic Diversity**: Considers how evolutionarily diverse a community is, using a phylogenetic tree.

In agriculture, alpha diversity can reflect the health or disturbance of a soil environment — for example, more diverse microbial communities are often associated with healthier or less disturbed soils.

## Beta Diversity: Between-Sample Differences

Beta diversity compares microbial composition across samples. It helps you answer questions like “Do different fertilizers cause different microbial communities?” or “Are microbes in root zones different from those in bulk soil?”

Common metrics include:

- **Bray-Curtis Dissimilarity**: Based on differences in counts.
- **Jaccard Distance**: Based on shared presence/absence of species.
- **UniFrac (weighted/unweighted)**: Measures how phylogenetically different two communities are.

QIIME 2 calculates these distances and uses techniques like Principal Coordinates Analysis (PCoA) to visualize differences.

## Ordination and Visualization

Ordination is a statistical method that reduces complex distance matrices into 2 or 3 dimensions for easier visualization. QIIME 2 uses **Emperor** to create 3D plots of PCoA results. These plots help identify clusters of similar samples, which can reflect treatment effects, environmental conditions, or sample types.

## Differential Abundance Analysis

Identifying which microbes are more or less abundant between groups (e.g., treated vs. untreated plants) is crucial. While traditional statistical methods like t-tests aren't suitable for microbiome data (because it's compositional and sparse), specialized tools like **Songbird**, **ANCOM**, or **DESeq2** can model these differences.

As an example, if a plant genotype promotes beneficial microbes that suppress disease, a differential abundance test might identify those beneficial taxa as more abundant in samples from that genotype.

## Why It Matters for Agriculture

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
