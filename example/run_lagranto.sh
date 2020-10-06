#!/bin/bash

module load dyn_tools
MODEL=mpas; export MODEL
LAGRANTO=/home/steidani/phd/mpas_lagranto/code/prog
which caltra
which seltra
module list
echo ${LAGRANTO}

#=======
# define paths
#pdir=/net/rossby/lhome/sprenger/lagranto.mpas/josh.alland    
pdir=/net/thermo/atmosdyn/steidani/mpas/data/MPAS_lagranto

#=======
# create link to files
echo " *********************** get files ***********************"	
ln -sf ${pdir}/* .

#create tracevars file
cat > tracevars << EOF
temperature 1. P
EOF

# calculate trajectories
echo " *********************** calc traj ***********************"
${LAGRANTO}/caltra startf_z -72 mpas.backward.1 -ref 19881231_1800 -o 360
