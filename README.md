MicEco: Various functions for analysis for 16S rRNA amplicon data, or similar.
------------------------------------------------------------------------------

#### rarefy\_rrna

This functions combines rarefaction with normalization. It simply
rarefies reads with a probability of the inverse 16S rRNA copy number,
such that besides rarefying the read counts the otu-table will be
corrected for the varying 16S rRNA copy numbers.

#### proportionality

Calculate proportionality on a phyloseq object or otu-table. Proposed by
Lovell et al. 2016 Proportionality: a valid alternative to correlation
for relative data
(<http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004075>)
