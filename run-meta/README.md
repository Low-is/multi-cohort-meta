# 📊 Bioinformatics Multi-Cohort Meta-Analysis Toolkit

This toolkit provides a framework for computing **effect sizes**, pooling **meta-analytic estimates**, and performing statistical analysis on **gene expression datasets**.

It helps identify genes or gene sets with **consistent differential expression across multiple studies**, supporting **biomarker discovery** and **predictive modeling**.

## Tools Used
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Git](https://img.shields.io/badge/Git-F05032?style=for-the-badge&logo=git&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-4EAA25?style=for-the-badge&logo=gnubash&logoColor=white)

---

# 📚 Statistical Concepts

## Standard Deviation (SD)

Measures how **spread out individual data points are around the mean**.

* A **large SD** → data points are widely spread.
* A **small SD** → data points are clustered close to the mean.

---

## Standard Error (SE)

Measures how much the **sample mean would vary if the experiment were repeated many times**.

* Represents the **standard deviation of the distribution of sample means**.
* Indicates the **precision of the sample mean estimate**.

---

# 📏 Effect Size

## What is Effect Size?

Effect size quantifies the **magnitude of a relationship or difference between groups**.

Unlike statistical significance, it indicates **practical significance**.

### Example

If:

* Treatment group mean weight loss = **10 kg**
* Control group mean weight loss = **5 kg**

Then:

```
Cohen's d ≈ 0.67
```

This represents a **medium effect size**, meaning the treatment had a noticeable effect.

---

# 📊 Cohen's d

Used to **compare the difference between two group means**.

### Formula

```
d = (M2 - M1) / SD_pooled
```

Where:

* **M1** = Mean of group 1
* **M2** = Mean of group 2
* **SD_pooled** = pooled standard deviation

### Interpretation

| Cohen's d | Effect Size |
| --------- | ----------- |
| 0.2       | Small       |
| 0.5       | Medium      |
| 0.8       | Large       |

---

# 📉 Hedges' g (Bias Correction)

Hedges' g corrects Cohen's d for **small sample sizes**.

### Correction Factor

```
cf = (1 - 3) / (4 * (n1 + n2) - 9)
```

### Hedges' g

```
g = cf * (mean difference) / (pooled standard deviation)
```

Where:

* **n1, n2** = sample sizes of the two groups

---

# 🔗 Pearson's r

Measures the **strength and direction of a linear relationship between two variables**.

### Range

```
-1  → Perfect negative relationship
0   → No relationship
+1  → Perfect positive relationship
```

---

# 🧮 Pooled Standard Deviation

A **weighted average of standard deviations from multiple groups**.

Used when assuming **equal variance between groups**.

### Formula

```
Pooled SD = sqrt((SD1² + SD2²) / 2)
```

Where:

* **SD1** = standard deviation of group 1
* **SD2** = standard deviation of group 2

---

# 📦 Summary Effect Size

A **combined effect size calculated across multiple studies**.

Commonly used in **meta-analysis** to produce a **single overall estimate of an effect**.

This allows researchers to:

* Integrate results from multiple datasets
* Increase statistical power
* Identify consistent biological signals

---

# 🧬 Applications in Gene Expression Analysis

This toolkit supports:

* Differential gene expression meta-analysis
* Cross-study biomarker discovery
* Robust effect size estimation
* Predictive modeling using gene expression profiles

---

⭐ Designed for **bioinformatics researchers**, **computational biologists**, and **data scientists working with genomic datasets**.


# Key Features 🔑:
[Add content here]


## Applications 🖥️:
  - Identification of robust biomarkers across independent studies.
  - Construction of predictive models for early diagnosis or prognosis.
  - Insight into consistent molecular responses to conditions or treatments despite study-specific variability.
  - Combinging studies improves statistical power in selecting gene signatures.


---

## Below is an example of the output:

<!-- Forest plots side by side -->
<table>
  <tr>
    <td><img src="docs/forest_plots/BOLA1.jpg" alt="BOLA1" width="250"></td>
    <td><img src="docs/forest_plots/CYP4F3.jpg" alt="CYP4F3" width="250"></td>
    <td><img src="docs/forest_plots/VEGFA.jpg" alt="VEGFA" width="250"></td>
  </tr>
</table>

---


---

## Future Directions 🚀:
  - Include mulitple feature selection methods, deafult setting is selecting shared genes across platforms before performing meta-analysis.
  - Additional machine learning algorithms optomized for diverse data types.
  - User-defined parameters will allow filtering by effect size thresholds and pooled effect size cutoffs.
  - Adapt to other data types such as clinical data, single-cell RNA-Seq, spatial, proteomics, metabolomics, microbiome, GWAS, etc. 
  - Planned support for more data soures, including ArrayExpress, TCGA, cBioPortal, etc.
  - Eventually expand for exploratory analysis (LDA clustering), regression and classification task using ML. 



# Citations: 
Zheng, H., Rao, A.M., Ganesan, A. et al. Multi-cohort analysis identifies a blood-based immune transcriptomic signature for early lung cancer detection. npj Precis. Onc. 9, 246 (2025).[https://doi.org/10.1038/s41698-025-01043-z](https://doi.org/10.1038/s41698-025-01043-z)

Walsh CJ, Batt J, Herridge MS, Mathur S, Bader GD, Hu P, Khatri P, Dos Santos CC. Comprehensive multi-cohort transcriptional meta-analysis of muscle diseases identifies a signature of disease severity. Sci Rep. 2022 Jul 4;12(1):11260.[https://www.nature.com/articles/s41598-022-15003-1](https://www.nature.com/articles/s41598-022-15003-1)

Rashid NU, Li Q, Yeh JJ, Ibrahim JG. Modeling Between-Study Heterogeneity for Improved Replicability in Gene Signature Selection and Clinical Prediction. J Am Stat Assoc. 2020;115(531):1125-1138.[https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/)

Pollard, Katherine S.; Dudoit, Sandrine; and van der Laan, Mark J., "Multiple Testing Procedures: R multtest Package and Applications to Genomics" (December 2004). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 164.[https://biostats.bepress.com/ucbbiostat/paper164](https://biostats.bepress.com/ucbbiostat/paper164)

Sweeney TE, Shidham A, Wong HR, Khatri P. A comprehensive time-course-based multicohort analysis of sepsis and sterile inflammation reveals a robust diagnostic gene set. Sci Transl Med. 2015 May 13;7(287):287ra71.[https://www.science.org/doi/10.1126/scitranslmed.aaa5993]([https://doi.org/10.1126/scitranslmed.aaa5993](https://www.science.org/doi/10.1126/scitranslmed.aaa5993))
