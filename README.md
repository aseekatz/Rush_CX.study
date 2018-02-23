## Rush Cross-sectional LTACH Study

### Code description

This repository includes processed sequence files, raw data files, sequence processing steps, and code used in the manuscript: "Gut microbiota and clinical features distinguish *Klebsiella pneumoniae* carbapenemase-producing *Klebsiella pneumoniae* colonization at the time of admission to a long-term acute care hospital", 2018 (submitted), authored by Anna M. Seekatz, Christine Bassis, Louis Fogg, Nicholas M Moore, Yoona Rhee, Karen Lolans, Robert A Weinstein, Michael Y Lin, Vincent B Young, Mary K Hayden, and for the Centers for Disease Control and Prevention Epicenters Program.

#### Importance:
We identified fecal microbiota features that distinguished *Klebsiella pneumoniae* carbapanemase-producing *K. pneumoniae* colonized and non-colonized patients who were newly admitted to a long-term acute care hospital. These features may warrant further investigation as candidates for inclusion in microbiome disruption indices.

#### Abstract:

**Background:** Identification of gut microbiota features associated with antibiotic-resistant bacterial colonization may reveal new infection prevention targets.

**Methods:** We conducted a matched, case-control study of long-term acute care hospital (LTACH) patients to identify gut microbiota and clinical features associated with colonization by *Klebsiella pneumoniae* carbapenemase-producing *Klebsiella pneumoniae* (KPC-Kp), an urgent antibiotic resistance threat. Fecal or rectal swab specimens were collected and tested for KPC-Kp; 16S rRNA gene-based sequencing was performed. Comparisons were made between cases and controls in calibration and validation subsamples using microbiota similarity indices, logistic regression, and unit-weighted predictive models. 

**Results:** Case (n=32) and control (n=99) patients had distinct fecal microbiota communities, but neither microbiota diversity nor inherent clustering into community types distinguished case and control specimens. Comparison of differentially abundant operational taxonomic units (OTUs) revealed one OTU associated with case status in both calibration (n=51) and validation (n=80) subsets that matched the canonical KPC-Kp strain ST258.  Permutation analysis using presence or absence of OTUs and hierarchical logistic regression identified two OTUs (belonging to genera Desulfovibrio and Ruminococcus) associated with KPC-Kp colonization. Among clinical variables, presence of a decubitus ulcer alone was independently and consistently associated with case status. Combining presence of OTUs *Desulfovibrio* and *Ruminococcus* with decubitus ulcer increased the likelihood of KPC-Kp colonization to >38% in a unit-weighted predictive model.

**Conclusions:** We identified microbiota and clinical features that distinguished KPC-Kp gut colonization in LTACH patients, a population particularly susceptible to KPC-Kp infection. These features may warrant further investigation as markers of risk for KPC-Kp colonization.

### Raw sequence data

Raw sequence data for this project can be retrieved from SRA (BioProject PRJNA428477, BioSample SAMN8292036-8292166)

### Contributors

Anna M. Seekatz, Christine Bassis, and Louis Fogg contributed to this code. Please contact aseekatz@umich.edu for questions.

### Repository descriptions

**Code:** Code for figures 2-4 in manuscript.

**Data:** Data files generated from R code:
- mothur-generated files
- LEfSe input/output files (also ran with mothur)
  - note: email for larger files: 
    - rushfinal.final.list 
    - rushfinal.trim.contigs.good.unique.good.filter.unique.fasta
    - rushfinal.final.0.03.fasta

**Figures:** Figures and supplementary files included in manuscript.
