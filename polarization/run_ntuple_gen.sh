#!/bin/bash -xe

sleep $(( ( RANDOM % 200 ) + 1 ))

JOBNUM=$1
OUTPUT_DIR=$2

WORKDIR=`pwd`

source /cvmfs/cms.cern.ch/cmsset_default.sh

############ Start DNNTuples ############
export SCRAM_ARCH=el8_amd64_gcc12

CMSSWDIR="."
scram p CMSSW CMSSW_15_0_18
cd CMSSW_15_0_18/src
eval `scram runtime -sh`

# clone this repo into "DeepNTuples" directory
git clone https://github.com/yiyangzha/Onia2MuMu.git Analyzers/MuMu

scram b -j8

cd $WORKDIR

function retry {
  local n=1
  local max=5
  local delay=5
  while true; do
    "$@" && break || {
      if [[ $n -lt $max ]]; then
        ((n++))
        echo "Command failed. Attempt $n/$max:"
        sleep $delay;
      else
        echo "The command has failed after $n attempts."
        return 1
      fi
    }
  done
}

seed=$(((${JOBNUM}) * 100))
seed=$((seed + ($(date +%s) % 10000)))

### process file
retry cmsRun ParticleGun-Upsilon2MM.py SEED=$seed
retry cmsRun $CMSSWDIR/CMSSW_15_0_18/src/Analyzers/MuMu/test/run_upsilon_condor_gen.py inputFiles=file:step.root
### end processing file

if ! [ -z "$OUTPUT_DIR" ]; then
  retry xrdcp --silent -p -f rootuple.root $OUTPUT_DIR
fi

touch ${WORKDIR}/dummy.cc
