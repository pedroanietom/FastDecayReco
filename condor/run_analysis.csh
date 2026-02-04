#!/bin/csh

set baseDir=`pwd`

setenv HOME /sphenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end

source /opt/sphenix/core/bin/sphenix_setup.csh -n ana.532

setenv MYINSTALL /sphenix/u/pnietomar/install
setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
setenv ROOT_INCLUDE_PATH $MYINSTALL/include:$ROOT_INCLUDE_PATH
source /opt/sphenix/core/bin/setup_local.csh $MYINSTALL

# -----------------------------
# Condor arguments
# -----------------------------
set runnumber  = $1
set nevents    = $2
set segment    = $3
set CALOFILE   = $4
set CLUSTERFILE = $5
set CLUSTERDIR = $6
set OUTDIR     = $7

# -----------------------------
# Output file
# -----------------------------
mkdir -p ${OUTDIR}
set OUTPUTFILE = "${OUTDIR}/Fun4All_TrackSeeding_${runnumber}_${segment}.root"

echo "===================================="
echo "Initializing Fun4All_TrackSeeding"
echo "Run number      : ${runnumber}"
echo "Segment         : ${segment}"
echo "Events          : ${nevents}"
echo "Calo file list  : ${CALOFILE}"
echo "Cluster file    : ${CLUSTERFILE}"
echo "Output dir      : ${CLUSTERDIR}"
echo "Output file     : ${OUTPUTFILE}"
echo "===================================="

cd ${baseDir}
pwd

# -----------------------------
# Run ROOT macro
# -----------------------------
root -b -q 'Fun4All_TrackSeeding.C('${runnumber}', '${segment}', '${nevents}', "'${CALOFILE}'", "'${CLUSTERFILE}'", "'${CLUSTERDIR}'", "'${OUTPUTFILE}'")'

