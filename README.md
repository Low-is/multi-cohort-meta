# Multi-Cohort Meta-Analysis

## Clone Repo
```
# Bash  

git clone https://github.com/Low-is/multi-cohort-meta.git
cd multi-cohort-meta
```

## Create Python venv
```
# Bash

python -m venv venv
source venv/Scripts/activate # Git Bash command

# Install dependencies
pip install -r miner/requirements.txt
```

## Install R packages
```
# Bash

Rscript meta/setup.R # only need to run this once! 
```

## Run Data Processing Pipeline
```
# Bash

cd miner
python main.py

cd ..
# Process expression matrices
RScript meta/expr_pipeline.R # Takes ~15 mintutes to run

# Process pData
RScript meta/extract_pData.R

# Add 'condition' column to pData
RScript meta/add_condition_column.R
```

## Run Meta-Analysis Pipeline (under development)
```
# Bash

RScript run-meta/run_meta.R 
```


# Citations: 
Zheng, H., Rao, A.M., Ganesan, A. et al. Multi-cohort analysis identifies a blood-based immune transcriptomic signature for early lung cancer detection. npj Precis. Onc. 9, 246 (2025).[https://doi.org/10.1038/s41698-025-01043-z](https://doi.org/10.1038/s41698-025-01043-z)

Walsh CJ, Batt J, Herridge MS, Mathur S, Bader GD, Hu P, Khatri P, Dos Santos CC. Comprehensive multi-cohort transcriptional meta-analysis of muscle diseases identifies a signature of disease severity. Sci Rep. 2022 Jul 4;12(1):11260.[https://www.nature.com/articles/s41598-022-15003-1](https://www.nature.com/articles/s41598-022-15003-1)

Rashid NU, Li Q, Yeh JJ, Ibrahim JG. Modeling Between-Study Heterogeneity for Improved Replicability in Gene Signature Selection and Clinical Prediction. J Am Stat Assoc. 2020;115(531):1125-1138.[https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/)

Pollard, Katherine S.; Dudoit, Sandrine; and van der Laan, Mark J., "Multiple Testing Procedures: R multtest Package and Applications to Genomics" (December 2004). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 164.[https://biostats.bepress.com/ucbbiostat/paper164](https://biostats.bepress.com/ucbbiostat/paper164)

Sweeney TE, Shidham A, Wong HR, Khatri P. A comprehensive time-course-based multicohort analysis of sepsis and sterile inflammation reveals a robust diagnostic gene set. Sci Transl Med. 2015 May 13;7(287):287ra71.[https://www.science.org/doi/10.1126/scitranslmed.aaa5993]([https://doi.org/10.1126/scitranslmed.aaa5993](https://www.science.org/doi/10.1126/scitranslmed.aaa5993))
