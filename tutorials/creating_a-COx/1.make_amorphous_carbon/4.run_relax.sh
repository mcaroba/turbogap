#!/usr/bin/bash

rm -f submit.sh

cat ../sample_submit_script.sh script_box_relax.sh > submit.sh

sbatch submit.sh
