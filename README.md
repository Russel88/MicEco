[![Travis Build
Status](https://travis-ci.org/Russel88/MicEco.svg)](https://travis-ci.org/Russel88/MicEco)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Package-License](http://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

MicEco: Various functions for analysis for microbial community data
-------------------------------------------------------------------

### Installation

    library(devtools)
    install_github("Russel88/MicEco")

### Citation
[![DOI](https://zenodo.org/badge/83547545.svg)](https://zenodo.org/badge/latestdoi/83547545)

## Phyloseq extensions
#### ps_prune

Prune taxa (ASVs, OTUs) from a phyloseq object based on their abundance and/or prevalence

#### ps_venn

Make Venn diagram of shared taxa (ASVs, OTUs) across sample groups from a phyloseq object. Overlap can be weighted by relative abundance

#### ps_euler

Make Euler diagram of shared taxa (ASVs, OTUs) across sample groups from a phyloseq object. Overlap can be weighted by relative abundance

#### ps_pheatmap

Make pretty heatmap directly from a phyloseq object. Built-in agglomoration, filtering, ordering, scaling, transformation, and annotation.

#### rcurve

Rarefaction curve (theoretical and fast) from a phyloseq object. Output ready for plotting in ggplot2

#### ps_tax_clean

Clean tax_table such that NAs are replaced with names of the most specific known taxonomy prefixed with the rank.

## Miscellaneous functions
#### clr

CLR transformation of community matrix, with multiplicative zero replacement

#### adonis_OmegaSq

Calculate the unbiased effect size estimation (partial) omega-squared for adonis (PERMANOVA) models. Note that the calculation is similar to a standard ANOVA and is not based on a theoretical foundation specifically for PERMANOVA.

#### WdS.test

Wd* - robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9). This code is taken from https://github.com/alekseyenko/WdStar/. An alternative to PERMANOVA.

#### UniFrac.multi

With unrooted phylogenies UniFrac sets the root randomly on the tree. 
The position of the root affects the results. 
This function runs UniFrac multiple times in parallel, with different roots, and takes the average to smooth potential bias.

#### proportionality

Calculate proportionality on a phyloseq object or otu-table. Proposed by
Lovell et al. 2016 Proportionality: a valid alternative to correlation
for relative data
(<http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075>)

## 16S rRNA gene copy number analyses
#### community\_rrna

Calculate the average 16S rRNA copy number of the OTUs in a
community/sample

Please cite [Schostag et al. 2019 ISMEJ](https://doi.org/10.1038/s41396-019-0351-x) if you use this function.

#### rarefy\_rrna

Combines rarefaction with 16S rRNA copy number correction. It rarefies
counts with a probability of the inverse 16S rRNA copy number, such that
besides rarefying the read counts the otu-table will be corrected for
the varying 16S rRNA copy numbers of the OTUs.

## Neutral model
#### neutral.fit

Fit neutral model developed by Sloan et al. (2006, Environ Microbiol 8(4):732-740) and implemented by Burns et al. (2015, ISME J 10(3):655-664).

#### neutral.rand

Fit neutral model developed by Sloan et al. (2006, Environ Microbiol 8(4):732-740) and implemented by Burns et al. (2015, ISME J 10(3):655-664) several times on ramdomly picked samples and with 16S rRNA gene copy number corrected rarefaction (rarefy_rrna).

## Beta diversity null models
#### ses.UniFrac

Standardized effect size of UniFrac, based on null models created with
permatfull/permatswap from the vegan package, or simple shuffling of
phylogenetic tree.

#### ses.comdist

Standardized effect size of MPD (mean pairwise distance) separating taxa
in two communities, a measure of phylogenetic beta diversity (also
called betaNRI and betaMPD). This is a combination of ses.mpd
(Standardized effect size of MPD in single communities) and comdist (MPD
between two communities) from the picante package.

#### ses.comdistnt

Standardized effect size of MNTD (mean nearest taxon distance)
separating taxa in two communities, a measure of phylogenetic beta
diversity (also called betaNTI and betaMNTD). This is a combination of
ses.mntd (Standardized effect size of MNTD in single communities) and
comdistnt (MNTD between two communities) from the picante package.

#### ses.comdist2

As `ses.comdist`, but null models are created with permatfull/permatswap
from the vegan package

#### ses.comdistnt2

As `ses.comdistnt`, but null models are created with
permatfull/permatswap from the vegan package

#### comdist.par

A parallel version of the comdist function from the picante package for significant speedup on large datasets

#### comdistnt.par

A parallel version of the comdistnt function from the picante package for significant speedup on large datasets

#### ses.mpd.par

A parallel version of the ses.mpd function from the picante package for significant speedup on large datasets

#### ses.mntd.par

A parallel version of the ses.mntd function from the picante package for significant speedup on large datasets

#### ses.permtest

Permutation test of z-matrix from `ses.comdist`, `ses.comdist2`, `ses.comdistnt`, `ses.comdistnt2` and `ses.UniFrac`.

## Copyright notice

`rarefy_rrna`: Some code is from vegan licensed under GPL-2
(<https://github.com/vegandevs/vegan>)

`ses.mpd.par`, `ses.mntd.par`, `comdist.par`, `comdistnt.par`, `ses.comdist`, `ses.comdist2`, `ses.comdistnt` and `ses.comdistnt2`: Some
code is from picante licensed under GPL-2
(<https://github.com/skembel/picante>)
