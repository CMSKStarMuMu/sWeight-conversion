#!/bin/bash

year=2016
subs=${1}
nsubs=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/sWeight-conversion
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

# if [ ! -r $SAMPLEDIR/${year}_data_beforsel.root ]; then
#     echo $SAMPLEDIR/${year}_data_beforsel.root not found
#     exit 1
# fi
if [ ! -r $HOME/merge_convert_sWeights ]; then
    echo $HOME/merge_convert_sWeights not found
    exit 1
fi

cp $HOME/merge_convert_sWeights .

echo ./merge_convert_sWeights ${subs} ${nsubs} ${year}
./merge_convert_sWeights ${subs} ${nsubs} ${year}
