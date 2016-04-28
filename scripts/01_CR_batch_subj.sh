#!/bin/bash -eu
# first arg: number of permutation runs; second is n bootstrap runs
# e.g., to execute from the command line: 
# scripts/DP_CR_analysis_overall_by_subj.sh 1000 1000

R_EXEC="R --no-save --no-restore"

# by subject
$R_EXEC --args $1 $2 < scripts/01_CR_subj_group.R
$R_EXEC --args $1 $2 < scripts/01_CR_subj_competition.R

$R_EXEC --args $1 $2 < scripts/01_CR_subj_consolidation.R
$R_EXEC --args $1 $2 < scripts/01_CR_consolidation_collapsed.R
