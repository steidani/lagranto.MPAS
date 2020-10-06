#!/bin/bash

module load dyn_tools
module list
#Currently Loaded Modulefiles:
#  1) legacy               4) pgi/13.9             7) dyn_tools
#  2) ncl/6.1.2            5) gcc/system
#  3) hdf5/1.8.13          6) netcdf/4.3.2-pgf90

LAGRANTO=/home/steidani/phd/mpas_lagranto/code/prog
echo ${LAGRANTO}

# =======
# define paths  
pdir=/net/thermo/atmosdyn/steidani/mpas/data/cdf		# location of mpas files
startdir=/net/thermo/atmosdyn/steidani/mpas/data/startf	# location of startfiles
tradir=/net/thermo/atmosdyn/steidani/mpas/data/tra	# working directory / location to store trajectories
if [ ! -d "$tradir" ]; then
    mkdir -p /net/thermo/atmosdyn/steidani/mpas/data/tra
fi

# =======						 
#create tracevars file
cat > $tradir/tracevars << EOF
temperature 1. P
EOF

# =======        
#select startfiles (startf_YYYMMDD_HH_FLAG)
cd ${startdir}
startfiles=`ls -v startf*`

# =======
# go to working directory
cd ${tradir}

# =======
# loop through startfiles (sorted with -v), select date, FLAG
for startf in ${startfiles}
do
    echo " *** $year *** "
    echo " *** $startf *** "

    # link startf 
    ln -sf ${startdir}/${startf} .

    # get date (YYYYMMDD_HH)
    date=`echo $startf | cut -c 8-18`
    
    # get FLAG
    bid=`echo $startf | cut -c 20-24`

    # =======
    # get MPAS files (7 days)
    for timestep in 0 24 48 72 96 120 144 168
    do
	newdate=`newtime $date -$timestep`	# remove "-" for forward trajectories
        
        year=`echo $newdate | cut -c 1-4`
        month=`echo $newdate | cut -c 5-6`
        day=`echo $newdate | cut -c 7-8`
        file=latlon.current.${year}-${month}-${day}_00:00:00.nc
	
	# create link
        if [ ! -f "$tradir/$file" ]; then
            ln -sf ${pdir}/${file} .
        fi
    done

    # =======
    # calculate trajectories
    echo " *********************** calc traj ***********************"
	
    ${LAGRANTO}/caltra ${startf} -168 lsl_${date}_${bid} -i 360 -o 360 -ref ${date}
    # 7-d (-168) backward trajectories from <startf>. remove "-" for forward trajectories
    # Temporal resolution of input is 6 h (360 min), und trajectories are written out every 6 h (360 min).

    # =======
    # clean up
    gzip ${tradir}/lsl_${date}_${bid}
    
    # remove old files
    rm ${tradir}/${startf}
    newdate=`newtime $date -192`	# -96 for 3d backward trajectories
    year=`echo $newdate | cut -c 1-4`
    month=`echo $newdate | cut -c 5-6`
    day=`echo $newdate | cut -c 7-8`
    file=latlon.current.${year}-${month}-${day}_00:00:00.nc
    if [ -f "$tradir/$file" ]; then
        rm ${tradir}/${file}
    fi
done   
