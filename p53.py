import scipy.stats as stats
import pandas as pd
import numpy as np

# 1. Load the cBioPortal mutation data
df = pd.read_csv(
    "data_mutations.txt",
    sep="\t",
    comment="#",
    low_memory=False
)

# 2. Collect unique TP53 mutations for each tumor
tp53 = df[df["Hugo_Symbol"] == "TP53"]

tumor_list = [
    (tumor, group["HGVSp_Short"].dropna().unique().tolist())
    for tumor, group in tp53.groupby("Tumor_Sample_Barcode")
]

# 3. Load TP53 functional annotations
tp53_data = pd.read_csv("MutationView_r21.csv")

tp53_small = tp53_data[
    ["ProtDescription", "TransactivationClass"]
].dropna(subset=["ProtDescription", "TransactivationClass"])

# 4. Inspect the available classifications
#print(tp53_small["TransactivationClass"].value_counts())

# 5. Select loss-of-function classifications
lof_classes = {"non-functional", "partially functional"}

lof_mutations = set(
    tp53_small.loc[
        tp53_small["TransactivationClass"]
        .str.strip()
        .str.lower()
        .isin(lof_classes),
        "ProtDescription"
    ]
)

# 6. Identify tumors with at least one qualifying TP53 mutation
lof_tumors = [
    tumor
    for tumor, mutations in tumor_list
    if any(mutation in lof_mutations for mutation in mutations)
]

# 7. Count tumors
n_lof_tumors = len(set(lof_tumors))

#print("Tumors with TP53 loss-of-function mutations:", n_lof_tumors)

# 8. Identify tumors without any p53 mutations
# All tumors represented in the mutation dataset
all_tumors = df["Tumor_Sample_Barcode"].dropna().unique()

# Tumors with at least one TP53 LOF mutation
lof_tumor_set = set(lof_tumors)

no_lof_tumors = [
    tumor for tumor in all_tumors
    if tumor not in lof_tumor_set
]
n_no_lof_tumors = len(set(no_lof_tumors))

# 8.1. Clean all tumors with mutation but their TransActivation class is NA
#turns out there are none


#8.2. It turns out not all cohort mutations are in database
cohort_mutations = {
    mutation
    for _, mutations in tumor_list
    for mutation in mutations
}

annotated_mutations = set(
    tp53_small["ProtDescription"].dropna()
)

#mutation ends up in unannotated if it is missing from TP53 database or its transactivation class is not stated
unannotated_mutations = cohort_mutations - annotated_mutations

#Turnication is a predicted loss of function, I add them to lof list
truncating = [
    mutation for mutation in unannotated_mutations
    if "*" in mutation or "fs" in mutation
]
truncating_tumors = {
    tumor
    for tumor, mutations in tumor_list
    if any(mutation in truncating for mutation in mutations)
}

lof_tumors = list(set(lof_tumors) | truncating_tumors)

#for the sake of simpler analysis I will exclude others missing from no_lof_tumors
#due to difficulty in identifying each individual case

other_missing = [
    mutation for mutation in unannotated_mutations
    if "*" not in mutation and "fs" not in mutation
]

other_missing_tumors = {
    tumor
    for tumor, mutations in tumor_list
    if any(mutation in other_missing for mutation in mutations)
}

no_lof_tumors = [
    tumor for tumor in all_tumors
    if tumor not in set(lof_tumors)
    and tumor not in other_missing_tumors
]

# 8.3. a bit of checking
num_of_lof_tumors = len(set(lof_tumors))
num_of_no_lof_tumors = len(set(no_lof_tumors))
#print(num_of_lof_tumors, num_of_no_lof_tumors, len(set(truncating_tumors)))
#print(set(lof_tumors)&set(no_lof_tumors))


#9. KRas

kras = df[df["Hugo_Symbol"] == "KRAS"]

tumor_list_kras = [
    (tumor, group["HGVSp_Short"].dropna().unique().tolist())
    for tumor, group in kras.groupby("Tumor_Sample_Barcode")
]
#print(tumor_list_kras)

#9.1 I need a file with a list of all unique KRAS mutations so that I can sort through them manually(
# Identify all unique KRAS mutations
kras = df[df["Hugo_Symbol"] == "KRAS"]

unique_kras_mutations = sorted(
    kras["HGVSp_Short"].dropna().unique()
)

#print("Number of unique KRAS mutations:", len(unique_kras_mutations))

#for mutation in unique_kras_mutations:
#    print(mutation)
# Save unique KRAS mutations to a file
#pd.DataFrame({
#    "HGVSp_Short": unique_kras_mutations
#}).to_csv("unique_KRAS_mutations.csv", index=False)


#9. KRas

kras = df[df["Hugo_Symbol"] == "KRAS"]

tumor_list_kras = [
    (tumor, group["HGVSp_Short"].dropna().unique().tolist())
    for tumor, group in kras.groupby("Tumor_Sample_Barcode")
]

#9.1 I need a file with a list of all unique KRAS mutations so that I can sort through them manually(
# Identify all unique KRAS mutations
unique_kras_mutations = sorted(
    kras["HGVSp_Short"].dropna().unique()
)

#print("Number of unique KRAS mutations:", len(unique_kras_mutations))

#for mutation in unique_kras_mutations:
#    print(mutation)
# Save unique KRAS mutations to a file
#pd.DataFrame({
#    "HGVSp_Short": unique_kras_mutations
#}).to_csv("unique_KRAS_mutations.csv", index=False)

# 9.2 Load manually prepared KRAS annotations
kras_data = pd.read_csv(
    "unique_KRAS_mutations.csv",
    sep=";",
    keep_default_na=False
)

# 9.3 Select GOF, no-GOF, and uncertain classifications
GOF_classes = {"Oncogenic", "Likely Oncogenic"}

GOF_mutations = set(
    kras_data.loc[
        kras_data["Oncogenicity"].str.strip().isin(GOF_classes),
        "HGVSp_Short"
    ].str.strip()
)

no_GOF_mutations = set(
    kras_data.loc[
        kras_data["Oncogenicity"].str.strip() == "no_GOF",
        "HGVSp_Short"
    ].str.strip()
)

NA_mutations = set(
    kras_data.loc[
        kras_data["Oncogenicity"].str.strip() == "NA",
        "HGVSp_Short"
    ].str.strip()
)

# 9.4 Identify tumors with at least one KRAS GOF mutation
GOF_tumors = {
    tumor
    for tumor, mutations in tumor_list_kras
    if any(mutation in GOF_mutations for mutation in mutations)
}

# 9.5 Identify tumors with at least one KRAS mutation annotated as NA
NA_tumors = {
    tumor
    for tumor, mutations in tumor_list_kras
    if any(mutation in NA_mutations for mutation in mutations)
}

# 9.6 Identify tumors with at least one KRAS mutation annotated as no_GOF
no_GOF_mutation_tumors = {
    tumor
    for tumor, mutations in tumor_list_kras
    if any(mutation in no_GOF_mutations for mutation in mutations)
}

# Include tumors without any KRAS mutation
all_tumors = set(df["Tumor_Sample_Barcode"].dropna())
kras_tumors = set(kras["Tumor_Sample_Barcode"].dropna())
no_KRAS_tumors = all_tumors - kras_tumors

# Combine no_GOF mutation tumors and tumors without KRAS mutations
not_gof_tumors = no_GOF_mutation_tumors | no_KRAS_tumors

# Exclude tumors with NA mutations, unless they have a known GOF mutation
not_gof_tumors = not_gof_tumors - NA_tumors - GOF_tumors


# 10 Preparation for Fisher analysis
# 10.1
TP53_LOF = set(lof_tumors)
TP53_not_LOF = set(no_lof_tumors)
TP53_known = TP53_LOF | TP53_not_LOF

# Use tumors with known status for both genes
eligible_tumors = TP53_known & (GOF_tumors | not_gof_tumors)

a = len(eligible_tumors & TP53_LOF & GOF_tumors)
b = len(eligible_tumors & TP53_LOF & not_gof_tumors)
c = len(eligible_tumors & TP53_not_LOF & GOF_tumors)
d = len(eligible_tumors & TP53_not_LOF & not_gof_tumors)

table = [[a, b],
         [c, d]]

print("Contingency table:", table)

odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")

print("Odds ratio:", odds_ratio)
print("P-value:", p_value)
print("Eligible tumors:", len(eligible_tumors))
