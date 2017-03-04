MicEco: Various functions for analysis for 16S rRNA amplicon data, or similar.
------------------------------------------------------------------------------

#### rarefy\_rrna

This functions combines rarefaction with normalization. It simply
rarefies reads with a probability of the inverse 16S rRNA copy number,
such that besides rarefying the read counts will be corrected for the
varying 16S rRNA copy numbers. `{r}rarefy_rrna`
