#!/usr/bin/bash

rm -f submit.sh

cat ../sample_submit_script.sh script_quench_randomised.sh > submit.sh

sbatch submit.sh
