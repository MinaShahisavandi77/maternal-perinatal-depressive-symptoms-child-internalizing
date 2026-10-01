# Maternal perinatal depressive symptoms impacts on child internalizing difficulties across development : How family based genomic risks design can clarify this relationship 

This repository contains the instructions and scripts that participating cohorts need to calculate polygenic risk scores (PRSs), build the analysis dataset, run the cohort-specific analyses, and share their results.

> **Before you start:** please check the latest version of the analysis plan and make sure your genotype data has gone through the same QC steps described there.

---

## Table of contents

- [Step 0 – Calculate polygenic risk scores](#step-0--calculate-polygenic-risk-scores)
  - [Weight files](#weight-files)
  - [File formats](#file-formats)
  - [Install PLINK](#install-plink)
  - [Run the PRS calculation](#run-the-prs-calculation)
  - [After calculating the PRSs](#after-calculating-the-prss)
- [Step 1 – Create the dataset](#step-1--create-the-dataset)
- [Step 2 – Run the analysis scripts](#step-2--run-the-analysis-scripts)
- [Step 3 – Share your results](#step-3--share-your-results)
- [Contact](#contact)

---

## Step 0 – Calculate polygenic risk scores

### Weight files

You will receive an email from **SURFfilesender** with the PRS weight files. Weights were derived with **LDpred2-auto**, and one set with **SBayesRC**.

| File | Phenotype | Method | Reference |
|---|---|---|---|
| `PPD_adjusted_weights_hapmap3_ldpred2.txt` | Postpartum depression (PPD) | LDpred2-auto | PMID:37849304 |
| `MDD_adjusted_weights_hapmap3_ldpred2.txt` | Major depressive disorder (MDD) | LDpred2-auto | PMID:39814019 |
| `MDD_female_adjusted_weights_hapmap3_ldpred2.txt` | MDD, females only | LDpred2-auto | PMID:40858613 |
| `PPD_unique_adjusted_weights_hapmap3_sbrc.txt` | PPD-specific (unique) component | SBayesRC | [See this repository](https://github.com/MinaShahisavandi77/Perinatal-Depression-as-a-Distinct-Genetic-Entity) |

### File formats

The weight files **do not contain a header**. The columns are as follows.

**LDpred2-auto files (6 columns)**

| Column | Name | Description |
|---|---|---|
| 1 | `CHR` | Chromosome |
| 2 | `POS` | Base-pair position |
| 3 | `SNPID` | Variant ID |
| 4 | `A1` | Effect allele |
| 5 | `A2` | Other allele |
| 6 | `BETA` | Adjusted weight (LDpred2-auto) |

**SBayesRC file**

| Column | Name | Description |
|---|---|---|
| 1 | `SNPID` | Variant ID |
| 2 | `A1` | Effect allele |
| 3 | `BETA` | Weight (SBayesRC) |
| 4 | `SE` | Standard error of the weight |
| 5 | `PIP` | Posterior inclusion probability |

### Install PLINK

PRSs are calculated with [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/).

- **Windows / macOS:** download PLINK from <https://www.cog-genomics.org/plink/>. Run it from a terminal, with the executable (`plink.exe` or `plink`) in the same working directory as the weight files.
- **Linux:** install it from the terminal:

  ```bash
  sudo apt install plink1.9
  ```

### Run the PRS calculation

Your genotype data must be in PLINK binary format (`.bed`, `.bim`, `.fam`), after your own QC and imputation. If your data is in VCF format, convert it first:

```bash
plink1.9 --vcf my_genetic_data.vcf.gz --keep-allele-order --make-bed --out my_genetic_data
```

Then calculate the scores:

```bash
# LDpred2-auto weights
plink1.9 --bfile my_genetic_data --score my_scores.txt 3 4 6 sum --out output_filename

# SBayesRC weights
plink1.9 --bfile my_genetic_data --score my_scores.txt 1 2 3 sum --out output_filename
```

| Argument | Meaning |
|---|---|
| `my_genetic_data` | Name of your genotype files, **without** extension |
| `my_scores.txt` | Name of the weight file, **with** extension |
| `output_filename` | Name of the output files, **without** extension |
| `3 4 6` | Columns for variant ID, effect allele, and weight in the LDpred2 files |
| `1 2 3` | Columns for variant ID, effect allele, and weight in the SBayesRC file |
| `sum` | Calculate the PRS as the **sum** of per-allele scores (PLINK uses the mean by default) |

The full script is available in [`PRS_calculation_file`](https://github.com/MinaShahisavandi77/maternal-perinatal-depressive-symptoms-child-internalizing/tree/main/scripts_PRS_calculation).

For more information on PLINK commands and flags, see the [PLINK 1.9 documentation](https://www.cog-genomics.org/plink/1.9/index).

### After calculating the PRSs

Standardize all scores to **mean = 0** and **SD = 1** before running any models. This step is also included in the shared R scripts.

---

## Step 1 – Create the dataset

### Inclusion criteria

Include participants who meet **all** of the following:

1. **European ancestry**
2. **Genotype data available for the child, mother, and father** (complete trios)
3. **At least one assessment of maternal depression:**
   - prenatal
   - postpartum, around 6 months
   - postpartum, around 36 months
4. **At least one assessment of child psychopathology:**
   - age 1–5 years
   - age 6–12 years
   - age 12–16 years

### Siblings

Keep only one child per family. Exclude siblings at random, or with another method that suits your cohort.

### Variables

Build the dataset following the variable [demo Excel file](https://github.com/MinaShahisavandi77/maternal-perinatal-depressive-symptoms-child-internalizing/blob/main/demo_variables_gene_mediation_project.xlsx), using **exactly the same variable names**.

---

## Step 2 – Run the analysis scripts
Before starting with scripts you need to install the R package designed specifically for this project [genmediation](https://github.com/MinaShahisavandi77/genmediation).
Each cohort has its own folder in this repository. Open the folder for your cohort and run the scripts in it.

> **Don't forget** to run the descriptive statistics scripts as well.

---

## Step 3 – Share your results

We only need the **`.csv` result files** produced by the scripts. 

The preferred way to send them is **SURFfilesender**. If you don't have access to it, contact us and we'll find another way.

---

## Contact

For questions or problems with the PRS calculation, please contact:

**M. Shahisavandi** – [m.shahisavandi@erasmusmc.nl](mailto:m.shahisavandi@erasmusmc.nl)
