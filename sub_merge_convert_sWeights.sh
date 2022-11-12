#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

for nsubs in {1..6}
do

    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_convert_sWeights.sub
Executable  = run_merge_convert_sWeights.sh
subs        = \$(ProcId)
Arguments   = \$INT(subs) ${nsubs}
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/merge_convert_sWeights_\$INT(subs)_${nsubs}.out
Error       = logs_parSub/merge_convert_sWeights_\$INT(subs)_${nsubs}.err
transfer_output_files = ""
+JobFlavour = "nextweek"
EOF

    echo "Queue ${nsubs}">>temp_sub_convert_sWeights.sub

    # Compilation, submission and file removal
    if make merge_convert_sWeights
    then condor_submit temp_sub_convert_sWeights.sub
    fi
    rm temp_sub_convert_sWeights.sub

done
