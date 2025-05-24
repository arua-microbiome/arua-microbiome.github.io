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

# Adding metadata and examining count tables
Once your ASV table has been generated, it needs to be connected to your sample metadata and taxonomic assignments. This step is essential for creating meaningful visual summaries and performing any ecological comparisons. Without metadata, your sequences are just anonymous counts. Metadata provides biological context — which genotype the sample comes from, which compartment (e.g. root zone, bulk soil), and which treatment it received. Taxonomy links each ASV to a known organism or group of organisms, letting us ask not just how communities vary, but who is driving that variation.

> QIIME 2 allows you to integrate the count table with your mapping file to visualise how read counts distribute across samples.
> ```bash
> qiime feature-table summarize \
>  --i-table wednesday_outputs/table-dada2.qza \
>  --m-sample-metadata-file wednesday_data/wednesday_metadata.csv \
>  --o-visualization thursday_outputs/table-dada2.qzv
> ```
> 
> Output: ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()


> To visualise taxonomic composition by group or treatment:
>```bash
> qiime taxa barplot \
>  --i-table wednesday_outputs/table-dada2.qza \
>  --i-taxonomy wednesday_outputs/taxonomy.qza \
>  --m-metadata-file wednesday_data/wednesday_metadata.tsv \
>  --o-visualization thursday_outputs/taxa-bar-plots.qzv
>```
>
> Output: ```taxa-bar-plots.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()

# Filtering contaminants
Sequencing from soil or root material often includes host DNA like mitochondria or chloroplasts. These non-microbial sequences can obscure patterns in microbial community composition and inflate diversity estimates. Filtering out these contaminants is a standard practice in microbiome workflows, especially when studying plant-associated microbes. This helps ensure that the dataset reflects only the true microbial community of interest.

Looking at the the ```taxonomy.qzv``` file using https://view/qiime2.org We can see the data presented at different taxonomic levels and grouped by different experimental factors. If we drill down to taxonomic level 5 something looks a bit odd. There's lots of "Rickettsiales;f__mitochondria".  This is really  plant mitochondrial contamination. Some of these samples also have chloroplast contamination.

>First remove unwanted taxa (like mitochondria and chloroplasts) from the ASV count table: 
>```bash
>qiime taxa filter-table \
>  --i-table wednesday_outputs/table-dada2.qza \
>  --i-taxonomy wednesday_outputs/taxonomy.qza \
>  --p-exclude mitochondria,chloroplast \
>  --o-filtered-table thursday_outputs/table-dada2-filtered.qza
>```
>
>Output: ```table-dada2-filtered.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()


>Then remove those same taxa from the actual DNA sequences of the ASVs, ensuring both abundance data and sequence data are clean for downstream analysis.
>```bash
>qiime taxa filter-seqs \
>  --i-sequences wednesday_outputs/rep-seqs-dada2.qza \
>  --i-taxonomy wednesday_outputs/taxonomy.qza \
>  --p-exclude mitochondria,chloroplast \
>  --o-filtered-sequences thursday_outputs/rep-seqs-dada2-filtered.qza
>```
>
>Output: ```rep-seqs-dada2-filtered.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()


>Since we have altered the qza file we can create a new bar plots:
>```bash
>qiime taxa barplot \
>  --i-table thursday_outputs/table-dada2-filtered.qza \
>  --i-taxonomy wednesday_outputs/taxonomy.qza \
>  --m-metadata-file wednesday_data/wednesday_metadata.tsv \
>  --o-visualization thursday_outputs/taxa-bar-plots-filtered.qzv
>```
>
>Output: ```taxa-bar-plots-filtered.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()

# Phylogenetics

A phylogenetic tree is a diagram that represents the evolutionary relationships among different organisms. By aligning and comparing the sequences, we can infer how closely related they are and construct a tree that reflects their shared evolutionary history. Some communities may appear different in terms of ASV presence but are composed of closely related taxa. To assess this, we construct a phylogenetic tree of all ASVs. 

This tree is essential for certain types of diversity metrics, such as UniFrac, which measure not just the presence or absence of microbes, but how evolutionarily different the communities are. Understanding these relationships helps us interpret microbial functions and ecological roles more effectively. It's also useful for certain types of diversity comparisons between samples.

## Multiple sequence alignment
We begin by using Mafft, which aligns all representative ASV sequences so that homologous nucleotide positions are lined up across sequences:

>```bash
>qiime alignment mafft \
>  --i-sequences rep-seqs-dada2-filtered.qza \
>  --o-alignment aligned-rep-seqs.qza
>```
>
>Output: ```aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

## Masking sites
Masking is the process of removing highly variable, gappy, or uninformative positions from a multiple sequence alignment. These positions often arise from sequencing noise, misalignments, or non-homologous regions, and can distort phylogenetic inference by introducing noise into the tree-building process.

>In QIIME 2, this command takes the aligned ASV sequences and filters out alignment columns (positions) that don't contain reliable, conserved sequence information—leaving a cleaner dataset for more accurate and robust tree construction:
>
>```bash
>qiime alignment mask \
>  --i-alignment aligned-rep-seqs.qza \
>  --o-masked-alignment masked-aligned-rep-seqs.qza
>```
>
>Output: ```masked-aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

## Creating a tree
FastTree builds a phylogenetic tree from the masked, aligned ASV sequences using an approximate maximum-likelihood method. This tree reflects the evolutionary relationships among ASVs and is essential for phylogeny-based diversity metrics like UniFrac.

>The resulting tree is unrooted, meaning it shows relationships but not direction of ancestry:
>
>```bash
> qiime phylogeny fasttree \
>  --i-alignment masked-aligned-rep-seqs.qza \
>  --o-tree unrooted-tree.qza
>```
>
>Output: ```unrooted-tree.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

## Midpoint rooting
Rooting the tree defines a starting point for evolutionary comparisons. Since our tree is initially unrooted (no known ancestor), we use midpoint rooting, which places the root at the midpoint of the longest distance between any two tips—giving a balanced view of divergence.

>```bash
>qiime phylogeny midpoint-root \
>  --i-tree unrooted-tree.qza \
>  --o-rooted-tree rooted-tree.qza
> ```
>
>Output: ```rooted-tree.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

## Visualisation
Once you’ve generated a rooted phylogenetic tree and cleaned, filtered abundance data, it’s time to visualise these relationships in an interactive, intuitive way. QIIME 2’s Empress plugin allows you to explore phylogenetic trees alongside sample metadata and taxonomic annotations. This is particularly powerful for understanding not only which taxa are present, but how community shifts map onto the evolutionary history of your organisms.

> The tree-plot command displays the phylogenetic tree and overlays taxonomic information. This is useful for checking the structure of the tree and seeing how ASVs group by lineage. It’s also a helpful tool to explore diversity patterns across evolutionary lineages.
>```bash
> qiime empress tree-plot \
>  --i-tree thursday_outputs/rooted-tree.qza \
>  --m-feature-metadata-file wednesday_outputs/taxonomy.qza \
>  --o-visualization thursday_outputs/empress-tree-tax.qzv
>```
>
>Output: ```empress-tree-tax.qzv``` [View]() \| [Download]()

>The community-plot command takes things further by integrating the phylogenetic tree, ASV abundance table, sample metadata, and taxonomy. This interactive plot lets you explore which lineages dominate particular samples or treatments, trace shifts in community composition, and highlight specific branches or taxa across groups.
>```bash
> qiime empress community-plot \
>  --i-tree thursday_outputs/rooted-tree.qza \
>  --i-feature-table thursday_outputs/table-dada2-filtered.qza \
>  --m-sample-metadata-file wednesday_data/wednesday_metadata.tsv \
>  --m-feature-metadata-file wednesday_outputs/taxonomy.qza \
>  --o-visualization thursday_outputs/community-empress-tree-tax.qzv
>```
>
>>Output: ```community-empress-tree-tax.qzv``` [View]() \| [Download]()



# Diversity

Microbial ecology is the study of how microbes interact with each other and their environment. In agricultural settings, these interactions can profoundly affect crop health, nutrient cycling, and disease suppression. Ecological analysis allows us to characterize microbial diversity, compare communities between samples (like different soil types or plant treatments), and identify patterns that inform how microbial communities function and respond to changes.

## Alpha Diversity: Within-Sample Diversity

Alpha diversity refers to the number and evenness of microbial species in a single sample. Three common metrics that Qiime2 can calculate are these:

>**Observed Features**: The count of unique sequence variants or species.
>
>```bash
>qiime diversity alpha \
>  --i-table thursday_outputs/table-dada2-filtered.qza \
>  --p-metric observed_features \
>  --o-alpha-diversity thursday_outputs/alpha-observed-features.qza
>```
>
>Output: ```alpha-observed-features.qza``` [View]() \| [Download]()

 
>**Shannon Diversity Index**: Takes into account both richness and evenness; higher values indicate more diverse communities.
>
>```bash
>qiime diversity alpha \
>  --i-table thursday_outputs/table-dada2-filtered.qza \
>  --p-metric shannon \
>  --o-alpha-diversity thursday_outputs/alpha-shannon.qza
>```
>
>Output: ```alpha-shannon.qza``` [View]() \| [Download]()

>**Faith’s Phylogenetic Diversity**: Considers how evolutionarily diverse a community is, using a phylogenetic tree.
>
>```bash
>qiime diversity alpha-phylogenetic \
>  --i-table thursday_outputs/table-dada2-filtered.qza \
>  --i-phylogeny thursday_outputs/rooted-tree.qza \
>  --p-metric faith_pd \
>  --o-alpha-diversity thursday_outputs/alpha-faith-pd.qza
>```
>
>Output: ```alpha-faith-pqd.qza``` [View]() \| [Download]()


As we have seen throughout the day, Qiime ```.qza``` don't contain data in an interpretable format. We thus want to export the alpha diversity metrics into three separate tables, using the ```export``` function.

>Repeat this for each sample:
>
>```bash
>qiime tools export \
>  --input-path thursday_outputs/alpha-<metric>.qza \
>  --output-path thursday_outputs/alpha-<metric>
>```
>
>Output: three directories containing alpha-diversity outputs.

## Beta Diversity: Between-Sample Differences

Beta diversity compares microbial composition across samples. It helps you answer questions like “Do different fertilizers cause different microbial communities?” or “Are microbes in root zones different from those in bulk soil?”

Common metrics include:

- **Bray-Curtis Dissimilarity**: Based on differences in counts.
- **Jaccard Distance**: Based on shared presence/absence of species.
- **UniFrac (weighted/unweighted)**: Measures how phylogenetically different two communities are.

QIIME 2 calculates these distances and uses techniques like Principal Coordinates Analysis (PCoA) to visualize differences.

# Ordination and Visualization

Ordination is a statistical method that reduces complex distance matrices into 2 or 3 dimensions for easier visualization. QIIME 2 uses **Emperor** to create 3D plots of PCoA results. These plots help identify clusters of similar samples, which can reflect treatment effects, environmental conditions, or sample types.

> **What is ordination?**
> It's a way to simplify and visualize complex data. For example, it helps you see patterns in how samples cluster based on their microbial composition.

Ordination is a dimensionality reduction technique that enables the visualization of sample differences. QIIME has a plugin called emperor that calculates a Bray-Curtis dissimilarity matrix and uses principal coordinates analysis (PCoA). you could also export the pcoa data and plot it yourself in the package of your choice.

 ```bash
 qiime emperor plot --i-pcoa core-diversity/bray_curtis_pcoa_results.qza --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt --o-visualization pcoa-visualization.qzv
 ```
 Time to run: 15 seconds

 Output:
 * ```pcoa-visualization.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()


# Differential Abundance Analysis

Identifying which microbes are more or less abundant between groups (e.g., treated vs. untreated plants) is crucial. While traditional statistical methods like t-tests aren't suitable for microbiome data (because it's compositional and sparse), specialized tools like **Songbird**, **ANCOM**, or **DESeq2** can model these differences.

As an example, if a plant genotype promotes beneficial microbes that suppress disease, a differential abundance test might identify those beneficial taxa as more abundant in samples from that genotype.

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
* ```composition.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()


Perform [Ward's agglomerative clustering](https://arxiv.org/abs/1111.6285)
```bash    
time qiime gneiss correlation-clustering \
    --i-table composition.qza\
    --o-clustering hierarchy.qza
```
Time to run: 5 minutes

Output:
* ```hierarchy.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

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
* ```balances.qza``` [View](https://view.qiime2.org/?src=) \| [Download]()

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
* ```regression_summary.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()


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
* ```linear_mixed_effects_model.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()


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
* ```y0_taxa_summary.qzv``` [View](https://view.qiime2.org/?src=) \| [Download]()


In this case the y0 balance is a split between samples that have plants in them and raw soil. It makes sense that this is the largest effect.  What happens if you run balance y2 or decrease the taxonomic level?

# Other analyses

This tutorial covered a range of analyses that can be done with microbiome data but there are other types on analyses that can be done too.

* Functional analysis - Several packages attempt to impute function from taxonomy including [PiCrust](https://picrust.github.io/picrust/), [Tax4fun](http://tax4fun.gobics.de/), [Piphillin](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166104)
* Inferring ecological interaction networks -[SPIEC-EASI](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). Co-Variance is an issue, but SPIEC-EASI attempts to model conditional independence.
* Data management tools - [Qiita](https://qiita.ucsd.edu/), and from ARS scientist Dan Manter [myPhyloDB](http://www.myphylodb.org/)
* Set analysis - From ARS Scientist Devin Coleman-Derr [MetaComet](https://probes.pw.usda.gov/MetaCoMET/)
* Other general analysis tools - [Mothur](https://www.mothur.org/) and the R-based [Phyloseq](https://joey711.github.io/phyloseq/)
