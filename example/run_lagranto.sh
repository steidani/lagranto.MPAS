#!/bin/bash

module load dyn_tools
module list
LAGRANTO=/home/steidani/phd/mpas_lagranto/code/prog
echo ${LAGRANTO}

#=======
# define paths  
pdir=/net/thermo/atmosdyn/steidani/mpas/data/cdf

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
${LAGRANTO}/caltra startf_z 72 lsl_forward.1 -ref 19881228_1800 -o 360
