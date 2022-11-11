#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_convert_sWeights.sub
Executable  = run_convert_sWeights.sh
subs        = \$(ProcId)
Arguments   = \$INT(subs)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/convert_sWeights_\$INT(subs).out
Error       = logs_parSub/convert_sWeights_\$INT(subs).err
transfer_output_files = ""
+JobFlavour = "nextweek"
EOF

echo 'Queue 8192'>>temp_sub_convert_sWeights.sub

# Compilation, submission and file removal
if make convert_sWeights
then condor_submit temp_sub_convert_sWeights.sub
fi
rm temp_sub_convert_sWeights.sub
