# TBB
R code for TBB correlative studies analyses
TBB = Talazoparib Beyond BRCA

This is a README file for the Repository. 

These are R code script files used to generate analyses for Gruber et al. Nature Cancer, _in press_ 
Article title: Talazoparib monotherapy in BRCA1 and BRCA2 wild-type patients with a mutation in other HR genes (Talazoparib Beyond BRCA)

Please cite our article if you use code or data from this Repository.  

Direct requests:

Joshua Gruber, MD, PhD,
Department of Internal Medicine,
Cecil H. and Ida Green Center for Reproductive Biology Sciences,
Simmons Comprehensive Cancer Center,
UT Southwestern Medical Center

email: joshua.gruber@UTsouthwestern.edu

lab website: https://www.gruberlabutsw.org

Article Abstract

Talazoparib, a PARP inhibitor, is active in germline (g) BRCA1/2 mutant advanced HER2-negative breast cancer, but its activity beyond gBRCA1/2 is poorly understood. We conducted an investigator-initiated single institution open-label, non-randomized phase II trial with a 2-stage design to evaluate talazoparib in patients with pre-treated advanced HER2-negative breast cancer or other solid tumors with germline or somatic alterations in homologous recombination (HR) pathway genes, not including BRCA1/2. Patients received talazoparib 1 mg orally daily until disease progression or unacceptable toxicity. The primary endpoint was objective response rate per RECIST 1.1. Twenty patients were enrolled in this cohort; 13 with breast and 7 non-breast cancer (pancreas, colon, uterine, testicular, parotid salivary). In patients with breast cancer, 4 had a RECIST partial response (ORR = 31%) and 3 additional patients had stable disease (SD) greater than 6 months (clinical benefit rate; CBR = 54%). All patients with gPALB2 mutations had treatment-associated tumor regression. Tumor HR deficiency score was positively correlated with magnitude of treatment response and was elevated in all gPALB2 tumors. In addition, ctHRD, a novel metric of HR deficiency calculated from plasma whole exome sequencing, was also correlated with treatment outcome, as was a gPALB2-associated mutational signature. Thus, talazoparib demonstrated activity in HER2-negative advanced breast cancer patients with gPALB2 mutations, showing its activity in HR pathway gene mutations beyond germline BRCA1/2. 

Sequencing data for this study, including clinical variables, is housed at dbGAP:

https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002803.v1.p1

All processed source data used to generate figures is also provided as Supplementary Table 5 in the publication.

Understanding the R scripts

There are 4 R scripts deposited, which each contain separate analyses used to generate conclusions and figures presented in the manuscript. Each script requires separate data input files.  These data input files are available for download in the repository. To follow the general flow of the paper the scripts should be run in the following order: 1) TBB_clinical_analysis.R, 2) Fig2.R, 3) TBB_myriad.R, 4) TBB_natera_ctDNA_github.R





