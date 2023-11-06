#!/usr/bin/bash

rm -f submit.sh

cat ../sample_submit_script.sh script_hamiltonian.sh > submit.sh

bash submit.sh
