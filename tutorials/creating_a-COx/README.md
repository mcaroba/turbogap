# Running this tutorial
Please see the Turbogap wiki [Creating oxygenated amorphous carbon](https://turbogap.fi/wiki/index.php/Creating_oxygenated_amorphous_carbon) for the detailed tutorial.

## If you are using this on a cluster
- Copy the directory `turbogap/tutorials/creating_a-COx` to where you
  want to run.
- Change `sample_submit_script.sh` to reflect the type of job
  scheduler you use (here, it's slurm), the modules you've loaded for
  `turbogap` and python, and change the `PATH` environment variable to
  where you've installed `turbogap/bin`.
- Make sure the project account is correct!
- Change the `srun turbogap` commands in the `script_*.sh` files to
  the standard for your cluster (e.g. `mpirun -np $N turbogap`).
- In each of the directories enumerated with "1.,2., etc", run the
  scripts in order after the preceding job has finished. They are
  enumerated with "1.,2., etc" with bash, e.g. `bash
  1.run_randomise.sh`.


## If you are using this locally
- Copy the directory `turbogap/tutorials/creating_a-COx` to where you
  want to run.
- Change `sample_submit_script.sh` by removing the `module load`
  lines, and change the `PATH` environment variable to where you've
  installed `turbogap/bin`.
- Change the `srun turbogap` commands in the `script_*.sh` files to
  how you run mpi commands (e.g. `mpirun -np $N turbogap`).
- Make sure your local python installation/venv has `numpy` and ASE
  installed `pip install ase`.
- In each of the directories enumerated with "1.,2., etc", run the
  scripts in order after the preceding job has finished. They are
  enumerated with "1.,2., etc" with bash, e.g. `bash
  1.run_randomise.sh`.
