# mutation_exclusivity
I wanted to test a statement I saw in "The Biology of Cancer" by Weinberg. It says that in some cancers RAS and TP53 mutations are exclusive, which sounded bizzare to me
## Data
I wasn't able to upload some files onto github. I will share them here
- Database for TP53 https://tp53.cancer.gov/get_tp53data look for file data_mutations.txt
- Database for KRas https://www.oncokb.org/ ended up sorting manually, look for [KRAS mutations](unique_KRAS_mutations.csv)
- Tumor cohort data https://www.cbioportal.org/study/clinicalData?id=coadread_dfci_2016

## Results so far (as per my first commit)
Analysis comes from Fisher exact test
First batch of data is from CRC (colorectal cancer)
- odds ratio: 0.8262032085561497
- p-value: 0.3150610985143039
Summary: so far inconclusive

## Future plans
- [ ] Figure out a suitable statistical model to test other types of cancer and summarize all of the data
- [ ] Get more data and build a pipeline  
