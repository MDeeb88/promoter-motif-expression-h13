# Computational Analysis of Transcription Factor Motif Landscapes and Gene Expression Prediction in Human Promoters and Intergenic Regions

## Project Overview

This project investigates how transcription factor (TF) motif composition is related to gene expression in human promoter regions and whether promoter-derived regulatory patterns can be used to predict transcriptional activity in intergenic regions.

The project has three main parts:

1. Motif enrichment analysis using Cohen's d effect size.
2. Gene expression prediction from promoter motif features using Elastic Net regression.
3. Prediction and characterization of intergenic regions with promoter-like transcriptional potential.

The central biological question is:

> Can motif architecture explain promoter activity, and can a promoter-trained model identify intergenic regions with promoter-like regulatory potential?

---

## Repository Structure

```text
data/
├── processed/
│   └── H13_promoters_651motifs_with_TPM.csv

figures/

notebooks/
├── 01_inspect_H13_new_matrices.ipynb
├── 02_H13_ElasticNet_analysis.ipynb
└── 03_H13_intergenic_prediction.ipynb

results/
├── H13_ElasticNet_coefficients.csv
├── H13_cohens_d_enrichment.csv
├── H13_motif_enrichment_with_ElasticNet_weights.csv
├── H13_intergenic_predicted_TPM.csv
└── H13_promoter_intergenic_PCA_predictions.csv
```

Large raw matrix files are excluded from Git using `.gitignore`.

---

## Input Data

The analysis uses updated H13 motif matrices:

- Human promoter motif matrix
- Human intergenic motif matrix

Each matrix contains:

- 651 TF motif features
- one region identifier column

The promoter matrix was merged with expression values from the previous project using gene symbols extracted from promoter identifiers.

Expression values are represented as:

```text
log2(TPM + 1)
```

After merging:

```text
Total new promoter rows: 20,090
TPM-matched promoter rows: 20,007
Missing TPM rows: 83
```

The 83 promoters without TPM values were excluded from model training.

---

## Notebook 1: Matrix Inspection

Notebook:

```text
01_inspect_H13_new_matrices.ipynb
```

This notebook checks the structure and quality of the new H13 motif matrices.

Main checks:

- matrix dimensions
- column consistency between promoter and intergenic matrices
- missing values
- negative values
- duplicate identifiers
- motif score distributions
- motif burden per region

Key finding:

The promoter and intergenic matrices have identical motif columns, allowing the same trained promoter model to be applied directly to intergenic regions.

---

## Notebook 2: Elastic Net Regression and Motif Enrichment

Notebook:

```text
02_H13_ElasticNet_analysis.ipynb
```

This notebook trains a regression model to predict promoter expression from motif composition.

### Elastic Net Model

Elastic Net regression was used because it combines two forms of regularization:

- **L1 regularization**, which can shrink some motif coefficients exactly to zero.
- **L2 regularization**, which stabilizes coefficients when motif features are correlated.

This is useful for motif-based genomics because TF motifs are often highly correlated or belong to related TF families.

The model pipeline was:

```text
Motif matrix
    ↓
StandardScaler
    ↓
Elastic Net regression
```

Selected model parameters:

```text
alpha = 0.006551
l1_ratio = 0.9
```

Approximate model performance:

```text
Pearson correlation ≈ 0.62
R² ≈ 0.38
RMSE ≈ 2.0
```

These results were consistent with the original motif-expression analysis.

---

## Cohen's d Motif Enrichment Analysis

Motif enrichment was measured using **Cohen's d**.

### What is Cohen's d?

Cohen's d is an effect size statistic that measures how strongly two groups differ relative to their variation.

For each motif, Cohen's d was calculated as:

```text
Cohen's d =
(mean motif score in high-expression promoters - mean motif score in low-expression promoters)
/
pooled standard deviation
```

In this project:

- A **positive Cohen's d** means the motif is more enriched in high-expression promoters.
- A **negative Cohen's d** means the motif is more enriched in low-expression promoters.
- A larger absolute value means a stronger difference between groups.

### Why Cohen's d works here

Cohen's d is appropriate here because the goal is not only to test whether a motif differs statistically between high- and low-expression promoters, but to measure how large and biologically meaningful that difference is.

This is important because the dataset contains thousands of promoters. With large datasets, very small differences can become statistically significant even when the biological effect is weak. Cohen's d focuses on the **magnitude of the effect**, not only on statistical significance.

Cohen's d also standardizes motif differences by their pooled standard deviation. This makes it useful for comparing motifs with different score distributions or different score ranges.

In practical terms, Cohen's d answers:

> Which motifs show the strongest standardized difference between highly expressed and weakly expressed promoters?

This makes it useful for ranking motifs by biological effect size.

### Enrichment Strategy

Promoters were ranked by expression.

The main enrichment comparison used:

- top 10% high-expression promoters
- bottom 10% low-expression promoters

Cohen's d was calculated for all 651 motifs.

A second expression-only enrichment analysis was also considered, where promoters with TPM = 0 were excluded before defining high- and low-expression groups. This was useful because many genes had exactly zero expression, which can make the low-expression group larger when using quantile thresholds.

The two enrichment views answer related but different biological questions:

1. Including TPM = 0 asks:

   > Which motifs distinguish active promoters from inactive or weak promoters?

2. Excluding TPM = 0 asks:

   > Among expressed promoters, which motifs distinguish high expression from low expression?

Both are biologically valid, but they should be reported separately because they ask different questions.

### Main High-Expression Motifs

Strongly enriched motifs in high-expression promoters included:

- SP140
- MBD1
- SP100
- MYCN
- E2F4
- ATF1
- GMEB1
- MYPOP
- KAISO

Several motifs from the original project remained highly ranked, including MBD1, SP100, MYPOP, KAISO, and ATF1.

Two highly ranked motifs were newly introduced in the expanded H13 motif matrix:

- SP140
- MYCN

This suggests that the new motif matrix preserved the original regulatory signal while adding additional informative motifs.

---

## Combining Cohen's d with Elastic Net Weights

Cohen's d and Elastic Net coefficients provide complementary information.

Cohen's d shows whether a motif differs between high- and low-expression promoters.

Elastic Net coefficients show whether the motif contributes to expression prediction in a multivariate model.

The strongest biological candidates are motifs that have:

- large positive Cohen's d,
- positive Elastic Net coefficient,
- non-zero Elastic Net weight.

This intersection helps identify motifs that are both enriched in highly expressed promoters and useful for predicting expression.

Output:

```text
results/H13_motif_enrichment_with_ElasticNet_weights.csv
```

---

## Notebook 3: Intergenic TPM Prediction

Notebook:

```text
03_H13_intergenic_prediction.ipynb
```

This notebook applies the promoter-trained Elastic Net model to the intergenic motif matrix.

Workflow:

```text
Train Elastic Net on promoter motifs + measured TPM
        ↓
Apply trained model to intergenic motif matrix
        ↓
Predict log2(TPM + 1) for intergenic regions
        ↓
Convert predictions to TPM scale
```

Output:

```text
results/H13_intergenic_predicted_TPM.csv
```

The prediction output includes:

- intergenic region ID
- predicted log2(TPM + 1)
- predicted TPM

---

## Intergenic Prediction Results

Most intergenic regions were predicted to have low transcriptional activity.

However, a subset of intergenic regions showed promoter-like predicted expression.

Summary of predicted intergenic log2(TPM + 1):

```text
count = 20,000
mean ≈ 0.62
median ≈ 0.49
max ≈ 10.83
```

The corresponding predicted TPM distribution was highly skewed, with a small subset of intergenic regions receiving very high predicted TPM values.

This suggests that most intergenic regions are transcriptionally weak, but some have motif architectures that resemble active promoters.

---

## PCA of Promoter and Intergenic Motif Landscapes

PCA was performed using only the 651 motif features.

Expression values were not used to compute the PCA.

### PCA by Region Type

Promoter and intergenic regions were projected into the same motif feature space.

Main observations:

- Promoters form a broad diagonal cloud in PCA space.
- Most intergenic regions overlap with the promoter motif landscape.
- A subset of intergenic regions forms distinct clusters, especially along PC2.

The first two principal components explained approximately:

```text
PC1 ≈ 33.6%
PC2 ≈ 13.2%
Total ≈ 46.9%
```

This indicates that a large portion of motif variation is captured in the first two components.

### PCA Colored by Expression

A second PCA plot colored points by expression:

- Promoters were colored by measured log2(TPM + 1).
- Intergenic regions were colored by predicted log2(TPM + 1).

This showed that high predicted intergenic expression was not randomly distributed. Instead, high-scoring intergenic regions tended to localize to specific areas of motif space.

### Top Predicted Intergenic Regions

The top 1% predicted intergenic regions were highlighted in PCA space.

These high-scoring regions formed localized clusters rather than being randomly scattered, suggesting that their high predictions are driven by specific motif architectures.

---

## Promoter vs Intergenic Prediction Distribution

Predicted expression distributions were compared between promoters and intergenic regions.

Main observations:

- Promoters had higher predicted expression overall.
- Intergenic regions had a lower central distribution.
- Intergenic predictions contained a right tail of promoter-like regions.

A Kolmogorov-Smirnov test comparing promoter and intergenic predicted expression distributions showed:

```text
KS statistic ≈ 0.546
p-value ≈ 0.0
```

The p-value appears as 0.0 due to floating-point underflow, meaning the difference is extremely significant.

This confirms that promoter and intergenic prediction distributions are statistically distinct.

---

## Biological Interpretation

The project supports several conclusions:

1. TF motif composition contains meaningful information about promoter expression.
2. Elastic Net regression can predict promoter activity from motif features with performance consistent with the original analysis.
3. Cohen's d identifies motifs strongly associated with high or low promoter expression.
4. The expanded H13 motif matrix preserves many motifs found in the original analysis while introducing new high-ranking motifs.
5. Most intergenic regions are predicted to have low expression.
6. A subset of intergenic regions shows promoter-like predicted activity.
7. PCA shows that high predicted intergenic regions occupy specific motif-space clusters, suggesting non-random regulatory architecture.

---

## Important Methodological Notes

### TPM = 0 promoters

The original promoter matrix contained many promoters with:

```text
log2(TPM + 1) = 0
```

These were retained in the main analysis for consistency with the original project and because zero-expression promoters are biologically meaningful.

However, zero-expression values create many ties in the bottom expression group. Therefore, an optional expression-only enrichment analysis can be performed by excluding TPM = 0 before ranking promoters. This analysis asks a different question and should be reported separately if used.

### Why Elastic Net was used instead of a black-box model

Black-box models may improve prediction, but this project requires interpretable motif weights.

Elastic Net provides:

- predictive performance,
- coefficient-based interpretation,
- feature selection through zero coefficients,
- compatibility with enrichment analysis.

This makes it suitable for connecting machine learning results to biological motif interpretation.

---

## Future Work

Possible next steps:

- Motif enrichment analysis of highly predicted intergenic regions.
- Overlap with enhancer annotations.
- Integration with ATAC-seq data.
- Integration with H3K27ac or H3K4me3 ChIP-seq.
- Comparison with CAGE peaks or eRNA annotations.
- Validation of top predicted intergenic regions as candidate regulatory elements.
- Functional clustering of TF motif families enriched in high-scoring intergenic regions.

---

## Software

Main Python libraries:

- pandas
- numpy
- matplotlib
- scikit-learn
- scipy

---

## Author

Mohammed Deeb

Master's Project  
Bioinformatics / Computational Genomics
