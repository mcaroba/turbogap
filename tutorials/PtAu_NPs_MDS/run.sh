# This value scales the annealing temperature
f=1.1

mkdir -p trajs/
mkdir -p logs/
mkdir -p trajs/trajs_$f
mkdir -p logs/logs_$f

for i in $(seq 1 1 660); do

SECONDS=0

n_Au=$(awk '{if($1=="Au"){n+=1} else {n+=0}} END {print n}' init_xyz/$i.xyz)
n_Pt=$(awk '{if($1=="Pt"){n+=1} else {n+=0}} END {print n}' init_xyz/$i.xyz)

T=$(echo $n_Au $n_Pt | awk -v f=$f '{print f*($1/($1+$2)*750.+$2/($1+$2)*1150.)}')

echo "Doing $i/660..."

##########################################
# Anneal at $T for 1 ps
cat>input<<eof
atoms_file = 'init_xyz/$i.xyz'
pot_file = 'gap_files/PtAuH.gap'

n_species = 3
species = Pt Au H

md_nsteps = 250
md_step = 4.

optimize = "vv"
thermostat = "bussi"
t_beg = $T
t_end = $T
tau_t = 10.
eof

mpirun -np 4 turbogap md &> /dev/null
n=$(head -1 trajectory_out.xyz | awk '{print $1+2}')
tail -$n trajectory_out.xyz > atoms.xyz
mv thermo.log logs/logs_$f/anneal_$i.log
##########################################



##########################################
# Quench to 100K over 1 ps
cat>input<<eof
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/PtAuH.gap'

n_species = 3
species = Pt Au H

md_nsteps = 250
md_step = 4.

optimize = "vv"
thermostat = "bussi"
t_beg = $T
t_end = 100.
tau_t = 100.

write_xyz = 250
eof

mpirun -np 4 turbogap md &> /dev/null
n=$(head -1 trajectory_out.xyz | awk '{print $1+2}')
tail -$n trajectory_out.xyz > atoms.xyz
mv trajectory_out.xyz trajs/trajs_$f/$i.xyz
mv thermo.log logs/logs_$f/quench_$i.log
##########################################



##########################################
# Relax
cat>input<<eof
atoms_file = 'atoms.xyz'
pot_file = 'gap_files/PtAuH.gap'

n_species = 3
species = Pt Au H

md_nsteps = 1000

optimize = "gd"
eof

mpirun -np 4 turbogap md &> /dev/null
n=$(head -1 trajectory_out.xyz | awk '{print $1+2}')
tail -$n trajectory_out.xyz >> trajs/trajs_$f/$i.xyz
mv thermo.log logs/logs_$f/relax_$i.log
rm trajectory_out.xyz
rm input
##########################################


echo "    ... in $SECONDS s"

done
