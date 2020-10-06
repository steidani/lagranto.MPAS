#!/bin/csh

#module load dyn_tools
module load netcdf/4.3.2-pgf90

# -----------------------------------------------------------------------------
# Set some general parameters
# -----------------------------------------------------------------------------

# Usage
if ( $#argv == 0 ) then
  echo "install.sh [lib|core|all|clean|test]"
  exit 0
endif

# Set the mode
set mode = $1

# Set path for development
set path_devel = "${PWD}/"

# Init Fortran compiler and set netCDF acccordingly
setenv FORTRAN pgf90
setenv FLAGS '-O'

# Set Lagranto path
setenv LAGRANTO ${path_devel}

# -----------------------------------------------------------------------------
# Set internal parameters and detailed installation mode
# -----------------------------------------------------------------------------

# Set netCDF paths
setenv NETCDF_LIB `nc-config --flibs`
setenv NETCDF_INC `nc-config --fflags`

# Set list of core programs
set core  = "caltra newtime" 

# Set list of libraries
set libs  = "iotra ioinp inter times rotate"

# Core programs
foreach prog ( $core )
   if ( "${prog}" == "${mode}" ) then
      set core  = ${prog}
      set mode  = "core"
   endif
end

# Libraries
foreach lib ( $libs )
   if ( "${lib}" == "${mode}" ) then
      set libs  = ${lib}
      set mode  = "lib"
   endif
end

# Check that the mode is ok 
if ( "${mode}" == "all"     ) goto modeok
if ( "${mode}" == "lib"     ) goto modeok
if ( "${mode}" == "core"    ) goto modeok
if ( "${mode}" == "clean"   ) goto modeok
if ( "${mode}" == "test"   ) goto modeok
echo "Unsupported mode ${mode} ... Stop"
exit 1

modeok:

# -----------------------------------------------------------------------------
# Make clean 
# -----------------------------------------------------------------------------
if ( "${mode}" == "clean" ) then

    cd ${LAGRANTO}/
    # clean programs
    foreach prog ( $core )
       \rm -f ${prog}/${prog} ${prog}/${prog}.o
    end

    # clean libraries
    \rm -f lib/*.a lib/*.o

    # clean test
    \rm -f test/*.log test/test.txt test/*.out
    exit 0

endif

# -----------------------------------------------------------------------------
# Install libraries 
# -----------------------------------------------------------------------------

if ( ("${mode}" == "lib") | ("${mode}" == "all" ) ) then

    echo "-----------------------------------------------------------------"
    echo "Installing libraries"
    echo "-----------------------------------------------------------------"
    echo
 
    # Change to library directory
    cd ${LAGRANTO}/lib

    # Loop over all libraries - compile and make library
    foreach lib ( $libs )
    
        \rm -f ${lib}.a
        \rm -f ${lib}.o
        if ( -f ${lib}.f ) then
          echo ${FORTRAN} -c ${FLAGS} ${lib}.f
          ${FORTRAN} -c ${FLAGS} ${NETCDF_INC} ${lib}.f
        endif
        if ( -f ${lib}.f90 ) then
          echo ${FORTRAN} -c ${FLAGS} ${lib}.f90
          ${FORTRAN} -c ${FLAGS} ${NETCDF_INC} ${lib}.f90
        endif
        ar r ${lib}.a ${lib}.o
        \rm -f ${lib}.l ${lib}.o
        ranlib ${lib}.a
        if ( ! -f ${lib}.a ) then
          echo "Problem in compiling ${lib} ... Stop"
          exit 1
        endif

    end

endif

if ( "${mode}" == "lib" ) exit 0

# -----------------------------------------------------------------------------
# Check that libraries are ok
# -----------------------------------------------------------------------------

echo
echo "-----------------------------------------------------------------"
echo "Check that all libraries are available"
echo "-----------------------------------------------------------------"
echo

# Change to library directory
cd ${LAGRANTO}/lib

# Check whether all libraries are available
foreach lib ( $libs )

    if ( ! -f ${lib}.a ) then
      echo "Library ${lib} missing... Stop"
      exit 1
    else
      ls -l ${lib}.a
    endif

end

# Exit if only libraries shoudl be installed
if ( "${mode}" == "lib" ) exit 0

# -----------------------------------------------------------------------------
# Compile Lagrango core programs
# -----------------------------------------------------------------------------

if ( ("${mode}" == "all" ) | ("${mode}" == "core" ) ) then

    echo
    echo "-----------------------------------------------------------------"
    echo "Installing Lagranto core programs"
    echo "-----------------------------------------------------------------"

    foreach prog ( $core )

        echo
        echo "----- $prog"
        echo
        cd ${LAGRANTO}/prog
        \rm -f ${prog}.o
        \rm -f ${prog}
        if ( "${prog}" == "trace"  ) \rm calvar.o

        if ( -f ${prog}.make ) then
           make -f ${prog}.make
        else if ( -f ${prog}.install ) then
           ./${prog}.install
        endif

        if ( ! -f ${prog} ) then
          echo "Problem in compiling ${prog} ... Stop"
          exit 1
        endif

    end

endif

exit 0
