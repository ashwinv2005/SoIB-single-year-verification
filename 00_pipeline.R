library(tidyverse)
library(glue)
library(tictoc)
# for parallel iterations
library(furrr)
library(parallel)

source("00_scripts/00_functions.R")


# PART 2 (subsample) ------------------------------------------------------------------

# Preparing data for trends analysis


# STEP 1: Subsample data for locations (create set of randomly selected GROUP.IDs)
# a file with random GROUP.IDs is first created so that the more time consuming step (creating the data files)
# can be repeated without sampling a different set of GROUP.IDs each time
# Run:
# - every time "sub_samp_locs.csv" is updated
# Requires:
# - tidyverse, parallel, foreach, doParallel, tictoc
# - data files:
#   - "sub_samp_locs.csv" for whole country and individual mask versions
# Outputs:
# - "randomgroupids.RData" for whole country and individual mask versions

load("00_data/analyses_metadata.RData")

# not functionising because parallelisation doesn't work inside functions
cur_mask <- "none"
tic("generated random group IDs for full country")
source("00_scripts/create_random_groupids.R")
toc() # 86 min

# STEP 2: Create subsampled data files using subsampled GROUP.IDs
# Run:
# - after above step (P2, S1)
# Requires:
# - tidyverse, tictoc
# - data files:
#   - "dataforanalyses.RData" for whole country and individual mask versions
#   - "randomgroupids.RData" for whole country and individual mask versions
# Outputs:
# - "dataforsim/dataX.csv" for whole country and individual mask versions

load("00_data/analyses_metadata.RData")

cur_mask <- "none"
my_assignment <- 1:50 # CHANGE FOR YOUR SUBSET
tic(glue("Generated subsampled data for full country (# {min(my_assignment)}:{max(my_assignment)})"))
source("00_scripts/create_random_datafiles.R")
toc() # 462 min (~ 8 h)






# PART 3 (run) ------------------------------------------------------------------

# STEP 1: Run trends models for all selected species
# Run:
# - after above step (P2, S2)
# Requires:
# - tidyverse, tictoc, lme4, VGAM, parallel, foreach, doParallel
# - data files:
#   - "dataforsim/dataX.csv" for whole country and individual mask versions
#   - "specieslists.RData" for whole country and individual mask versions
# Outputs:
# - "trends/trendsX.csv" for whole country and individual mask versions

load("00_data/analyses_metadata.RData")

# not functionising because parallelisation doesn't work inside functions
cur_mask <- "none"
my_assignment <- 1:2 # CHANGE FOR YOUR SUBSET
tic(glue("Species trends for full country (sims {min(my_assignment)}--{max(my_assignment)})"))
source("00_scripts/run_species_trends.R")
toc() # 102 hours