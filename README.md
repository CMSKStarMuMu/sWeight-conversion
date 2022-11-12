# sWeight-conversion
Code to convert the data weights from a sPlot to be positive

## Quick run

### Set up working area

```sh
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/ && cmsenv && cd ../..
git clone git@github.com:CMSKStarMuMu/sWeight-conversion.git
cd sWeight-conversion
```

Adjust the paths to input and output files in (convert_sWeights.cc), (merge_convert_sWeights.cc), (merge2_convert_sWeights.cc), (run_convert_sWeights.sh), (run_merge_convert_sWeights.sh)

### Submit first iteration

```sh
source sub_convert_sWeights.sh
```

### Submit second iteration

Adjust the `nsubs` variable in (sub_merge_convert_sWeights.sh), to specify the number of parallel jobs to be used in the splitting

```sh
source sub_merge_convert_sWeights.sh
```

### Merge output

```sh
make merge2_convert_sWeights
./merge2_convert_sWeights <NJOBS>
```

where `<NJOBS>` is the value set for the `nsubs` variable in the previous step
