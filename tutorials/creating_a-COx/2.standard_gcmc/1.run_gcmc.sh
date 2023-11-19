#!/usr/bin/bash

cat ../sample_submit_script.sh script_standard_gcmc.sh > submit.sh

sbatch submit.sh
