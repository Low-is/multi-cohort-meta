#!/bin/bash

set -e

echo "🔎 Running GEO mining pipeline..."
python biomedical-literature-miner/main.py

echo "🧬 Running R meta-analysis..."
Rscript multi-cohort-meta/meta_analysis.R

echo "✅ Full pipeline complete"
