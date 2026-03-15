#!/bin/bash

#Exit immediately if any command fails
set -e

#Configuring variables
REPO_DIR='/Users/lnash1/Documents/Birds of Norfolk/Web Dev'
SCRIPT="API Data Extraction.R"
DATA="All Observations.csv"
COMMIT_MSG="Automated data collection update $(date '+%d/%m/%Y')"

#Run script
cd "$REPO_DIR"
Rscript "$SCRIPT"

#Commit to GitHub
git add "$DATA"

if git diff --cached --quiet; then
    echo "No changes to commit."
else
    git commit -m "$COMMIT_MSG"
    git push origin main
fi
