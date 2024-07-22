#!/usr/bin/bash

cat ../sample_submit_script.sh script_randomise.sh > submit.sh

sbatch submit.sh
