# sc_virus
Regulatory motif analysis workflow in a single-cell setup

## Set up environment
1. To set up the environment for this project, you'll need to have conda installed.
2. Create the conda environment running the following:
   ```bash
   conda env create -f environment.yaml
3. Activate the environment and check the installation
   ```bash
   conda activate sc_virus
   conda list

## Running miReact
Specify according paths in `code/02run_mireact.R` and run `code/02run_mireact.sh` in a cluster. 

You can see more on the tutorial on how to run miReact here: https://github.com/muhligs/miReact

## Running the workflow
Main workflow for `sc_virus` can be found in the snakefile. Scripts can be found in `code/flow`

Running the workflow in the cluster:
   ```bash
   snakemake -j 10 --slurm --default-resources slurm_account=your_account


add `-n` flag for a dry-run
