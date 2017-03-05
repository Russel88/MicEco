MicEco: Various functions for analysis for microbial community data
-------------------------------------------------------------------

### Installation

    library(devtools)
    install_github("Russel88/MicEco")

#### rarefy\_rrna

This functions combines rarefaction with normalization. It rarefies
counts with a probability of the inverse 16S rRNA copy number, such that
besides rarefying the read counts the otu-table will be corrected for
the varying 16S rRNA copy numbers of the OTUs.

#### proportionality

Calculate proportionality on a phyloseq object or otu-table. Proposed by
Lovell et al. 2016 Proportionality: a valid alternative to correlation
for relative data
(<http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075>)

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
