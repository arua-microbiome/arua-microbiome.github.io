---
title: "Ecological analysis with QIIME2"
excerpt: "An example workflow using QIIME2 version 2024.2"
layout: single
mathjax: true
authors: "George Kalogiannis & Balig Panossian"

---

By George Kalogiannis & Balig Panossian, Designed from the official [QIIME2 tutorials](https://docs.qiime2.org/2024.2/tutorials/)
{% include toc %}

You’ve reached the point in the workflow where you can confidently inspect amplicon reads, clean them up, and generate reliable ASV tables and taxonomic summaries. In other words, you now know what sequences are present and how many of each occur in every sample. The next step is to turn those counts into ecological insight: measuring within-sample diversity, comparing whole communities across treatments, and linking patterns to phylogenetic relationships. The sections that follow introduce the core ecological tools you’ll use—alpha and beta diversity metrics, phylogeny-based distances, and statistical tests that reveal how environmental or experimental factors shape the microbiome.

# Learning Outcomes
By the end of this session you will be able to:
- Compute alpha- and beta-diversity metrics and visualise them in interactive .qzv files
- Construct a phylogenetic tree and run UniFrac-based ordination and PERMANOVA tests
- Export QIIME artifacts to R with and run statistical analyses on 

# Directories
Yesterday you worked in the directory called ```wednesday_outputs```. We will be using the files you created for today's ecological analysis and will output files to a directory called ```thursday_outputs```.

>```bash
>mkdir thursday_outputs
>```

# Phylogenetics

A phylogenetic tree is a diagram that represents the evolutionary relationships among different organisms. By aligning and comparing the sequences, we can infer how closely related they are and construct a tree that reflects their shared evolutionary history. Some communities may appear different in terms of ASV presence but are composed of closely related taxa. To assess this, we construct a phylogenetic tree of all ASVs. 

This tree is essential for certain types of diversity metrics, such as UniFrac, which measure not just the presence or absence of microbes, but how evolutionarily different the communities are. Understanding these relationships helps us interpret microbial functions and ecological roles more effectively. It's also useful for certain types of diversity comparisons between samples.

## Multiple sequence alignment
We begin by using Mafft, which aligns all representative ASV sequences so that homologous nucleotide positions are lined up across sequences:

>```bash
>qiime alignment mafft \
>  --i-sequences wednesday_outputs/rep-seqs-dada2-filtered.qza \
>  --o-alignment thursday_outputs/aligned-rep-seqs.qza
>```
>
>Output: ```aligned-rep-seqs.qza``` 

## Masking sites
Masking is the process of removing highly variable, gappy, or uninformative positions from a multiple sequence alignment. These positions often arise from sequencing noise, misalignments, or non-homologous regions, and can distort phylogenetic inference by introducing noise into the tree-building process.

>In QIIME 2, this command takes the aligned ASV sequences and filters out alignment columns (positions) that don't contain reliable, conserved sequence information—leaving a cleaner dataset for more accurate and robust tree construction:
>
>```bash
>qiime alignment mask \
>  --i-alignment thursday_outputs/aligned-rep-seqs.qza \
>  --o-masked-alignment thursday_outputs/masked-aligned-rep-seqs.qza
>```
>
>Output: ```masked-aligned-rep-seqs.qza```

## Creating a tree
FastTree builds a phylogenetic tree from the masked, aligned ASV sequences using an approximate maximum-likelihood method. This tree reflects the evolutionary relationships among ASVs and is essential for phylogeny-based diversity metrics like UniFrac.

>The resulting tree is unrooted, meaning it shows relationships but not direction of ancestry:
>
>```bash
> qiime phylogeny fasttree \
>  --i-alignment thursday_outputs/masked-aligned-rep-seqs.qza \
>  --o-tree thursday_outputs/unrooted-tree.qza
>```
>
>Output: ```unrooted-tree.qza```

## Midpoint rooting
Rooting the tree defines a starting point for evolutionary comparisons. Since our tree is initially unrooted (no known ancestor), we use midpoint rooting, which places the root at the midpoint of the longest distance between any two tips—giving a balanced view of divergence.

>```bash
>qiime phylogeny midpoint-root \
>  --i-tree thursday_outputs/unrooted-tree.qza \
>  --o-rooted-tree thursday_outputs/rooted-tree.qza
> ```
>
>Output: ```rooted-tree.qza``` 

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
>Output: ```empress-tree-tax.qzv```

>The community-plot command takes things further by integrating the phylogenetic tree, ASV abundance table, sample metadata, and taxonomy. This interactive plot lets you explore which lineages dominate particular samples or treatments, trace shifts in community composition, and highlight specific branches or taxa across groups.
>```bash
> qiime empress community-plot \
>  --i-tree thursday_outputs/rooted-tree.qza \
>  --i-feature-table wednesday_outputs/table-dada2-filtered.qza \
>  --m-sample-metadata-file /mnt/lustre/groups/WCHPC/wednesday_data/wednesday_metadata.tsv \
>  --m-feature-metadata-file wednesday_outputs/taxonomy.qza \
>  --o-visualization thursday_outputs/community-empress-tree-tax.qzv
>```
>
>Output: ```community-empress-tree-tax.qzv```


# Diversity

Microbial ecology is the study of how microbes interact with each other and their environment. In agricultural settings, these interactions can profoundly affect crop health, nutrient cycling, and disease suppression. Ecological analysis allows us to characterize microbial diversity, compare communities between samples (like different soil types or plant treatments), and identify patterns that inform how microbial communities function and respond to changes.

## Alpha Diversity: Within-Sample Diversity

Alpha diversity refers to the number and evenness of microbial species in a single sample. Three common metrics that Qiime2 can calculate are these:

>**Observed Features**: The count of unique sequence variants or species.
>
>```bash
>qiime diversity alpha \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --p-metric observed_features \
>  --o-alpha-diversity thursday_outputs/alpha-observed-features.qza
>```
>
>Output: ```alpha-observed-features.qza```

 
>**Shannon Diversity Index**: Takes into account both richness and evenness; higher values indicate more diverse communities.
>
>```bash
>qiime diversity alpha \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --p-metric shannon \
>  --o-alpha-diversity thursday_outputs/alpha-shannon.qza
>```
>
>Output: ```alpha-shannon.qza```

>**Faith’s Phylogenetic Diversity**: Considers how evolutionarily diverse a community is, using a phylogenetic tree.
>
>```bash
>qiime diversity alpha-phylogenetic \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --i-phylogeny thursday_outputs/rooted-tree.qza \
>  --p-metric faith_pd \
>  --o-alpha-diversity thursday_outputs/alpha-faith-pd.qza
>```
>
>Output: ```alpha-faith-pqd.qza```


As we have seen throughout the day, Qiime ```.qza``` files don't contain data in an interpretable format. We thus want to export the alpha diversity metrics into three separate tables, using the ```export``` function.

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

Beta diversity compares microbial composition across samples. It helps you answer questions like “Do different fertilizers cause different microbial communities?” or “Are microbes in root zones different from those in bulk soil?”. Common metrics include:

>**Bray-Curtis Dissimilarity**: Based on differences in counts.
>This abundance-based metric calculates dissimilarity between two communities based on the counts of shared features.
>```bash
>qiime diversity beta \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --p-metric braycurtis \
>  --o-distance-matrix thursday_outputs/beta-braycurtis.qza
>```
>
>Output: ```beta-braycurtis.qza```

>**Jaccard Distance**: Based on shared presence/absence of species.
>This metric compares samples based only on presence or absence of features, ignoring their abundance.
>```bash
>qiime diversity beta \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --p-metric jaccard \
>  --o-distance-matrix thursday_outputs/beta-jaccard.qza
>```
>
>Output: ```beta-jaccard.qza```

>**UniFrac (weighted/unweighted)**: Measures how phylogenetically different two communities are.
>
>a. Unweighted UniFrac
>```bash
>qiime diversity beta-phylogenetic \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --i-phylogeny thursday_outputs/rooted-tree.qza \
>  --p-metric unweighted_unifrac \
>  --o-distance-matrix thursday_outputs/beta-unweighted-unifrac.qza
>```
>
>Output: ```beta-unweighted-unifrac.qza```
>
>b. Weighted UniFrac
>```bash
>qiime diversity beta-phylogenetic \
>  --i-table wednesday_outputs/table-dada2-filtered.qza \
>  --i-phylogeny thursday_outputs/rooted-tree.qza \
>  --p-metric weighted_unifrac \
>  --o-distance-matrix thursday_outputs/beta-weighted-unifrac.qza
>```
>
>Output: ```beta-weighted-unifrac.qza```

>Lets export the Weighted UniFrac beta-diversity so we can then use it in R
>```bash
>qiime tools export \
>  --input-path friday_outputs/beta-weighted-unifrac.qza \
>  --output-path friday_outputs/beta-weighted-unifrac
>```
>
>Output: ```beta-braycurtis/distance-matrix.tsv```


# Statistical Analysis of Diversity

It's time to move the data outputs into the ARUA directory we were working in on Tuesday and **open up R**. Place the data into the ```data``` directory, we will be running our code from the ```code``` directory. You can bring the data onto your local machine using ```sftp```.

>```r
># Load libraries
>library(dplyr)
>library(vegan)
>library(readr)
>library(ggplot2)
>set.seed(719885)
>
># Load data
>metadata <- read.csv("../data/wednesday_metadata.tsv", sep = "\t")
>alpha <- read.csv("../data/alpha-diversity.tsv", sep = "\t")
>distance_matrix <- read.csv("../data/distance-matrix.tsv", sep = "\t")
>
># Prepare and merge metadata
>colnames(alpha)[1] <- "sample.id"
>
>metadata <- metadata %>%
>  filter(Genotype != "Soil") %>%
>  inner_join(alpha, by = "sample.id") %>%
>  mutate(
>    GenotypeGroup = case_when(
>      Genotype == "phr1" ~ "PHR1",
>      Genotype == "SPX1/SPX2" ~ "SPX",
>      TRUE ~ "non_psr"
>    ),
>    GenotypeGroup = factor(GenotypeGroup),
>    Genotype = factor(Genotype),
>    Experiment = factor(Experiment)
>  )
>
># Match sample IDs with distance matrix
>valid_ids <- intersect(rownames(distance_matrix), metadata$sample.id)
>distance_matrix <- distance_matrix[valid_ids, valid_ids]
>metadata <- metadata %>% filter(sample.id %in% valid_ids) %>%
>  arrange(factor(sample.id, levels = rownames(distance_matrix)))
>
># Alpha diversity model
>summary(lm(shannon_entropy ~ GenotypeGroup + Experiment, data = metadata))
>
># Beta diversity (PERMANOVA)
>adonis_result <- adonis2(as.dist(distance_matrix) ~ GenotypeGroup + Experiment, data = metadata)
>print(adonis_result)
>
># Plot: Alpha diversity boxplot
>ggplot(metadata, aes(x = GenotypeGroup, y = shannon_entropy, fill = GenotypeGroup)) +
>  geom_boxplot(alpha = 0.8) +
>  labs(title = "Alpha Diversity (Shannon Index)",
>       x = "Genotype Group",
>       y = "Shannon Entropy") +
>  theme_minimal() +
>  theme(legend.position = "none",
>        axis.text.x = element_text(angle = 45, hjust = 1))
>```







# Other analyses

This tutorial covered a range of analyses that can be done with microbiome data but there are other types on analyses that can be done too.

* Functional analysis - Several packages attempt to impute function from taxonomy including [PiCrust](https://picrust.github.io/picrust/), [Tax4fun](http://tax4fun.gobics.de/), [Piphillin](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166104)
* Inferring ecological interaction networks -[SPIEC-EASI](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). Co-Variance is an issue, but SPIEC-EASI attempts to model conditional independence.
* Data management tools - [Qiita](https://qiita.ucsd.edu/), and from ARS scientist Dan Manter [myPhyloDB](http://www.myphylodb.org/)
* Set analysis - From ARS Scientist Devin Coleman-Derr [MetaComet](https://probes.pw.usda.gov/MetaCoMET/)
* Other general analysis tools - [Mothur](https://www.mothur.org/) and the R-based [Phyloseq](https://joey711.github.io/phyloseq/)
