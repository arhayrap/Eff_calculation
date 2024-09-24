#!/bin/bash

version="UL"
current_dir=${PWD}

certificate="x509up_u123840"

output_folder="/eos/user/a/arhayrap/JPsi_${version}/"
executable="run_${version}.csh"
submitfile="Run_${version}.sub"

mkdir -p $output_folder
rm $executable
rm $submitfile

printf "#!/bin/tcsh \n\n" >> $executable

printf "setenv current_DIr `pwd` \n" >> $executable
printf "setenv X509_USER_PROXY ~/public/${certificate} \n" >> $executable
printf "# source /afs/cern.ch/cms/cmsset_default.csh \n" >> $executable
printf "source /cvmfs/cms.cern.ch/cmsset_default.csh \n" >> $executable
printf "cmsrel CMSSW_14_0_5 \n" >> $executable
printf "cd ./CMSSW_14_0_5/src \n" >> $executable
printf "cp ~/public/${certificate} /tmp/ \n" >> $executable
printf "cmsenv \n" >> $executable
printf "set newid=`expr \$1 + 1` \n" >> $executable
printf "cp -r -f ${current_dir}/../../PreSelection . \n" >> $executable
printf "cp -r -f ${current_dir}/Samples_JPsi_${version} ./sample_names \n" >> $executable
printf "cp       ${current_dir}/../../Trig_ScaleFactors_OnMC.py . \n" >> $executable
printf "scram b \n" >> $executable
printf "set filelist=( ./sample_names/* ) \n" >> $executable
printf "echo \${filelist[\$newid]} \n" >> $executable
printf "mv \${filelist[\$newid]} BToKPhi.txt \n" >> $executable
printf "set filenames = \"\" \n" >> $executable
printf "foreach line ( \"\`cat ./BToKPhi.txt\`\" ) \n" >> $executable
printf "    set filenames = \"\$filenames \$line\" \n" >> $executable
printf "end \n" >> $executable
printf "cmsRun Trig_ScaleFactors_OnMC.py \$filenames \n" >> $executable
printf "cp Trigger_ScaleFactors_OnMC.root ${output_folder}/output_file_\$1.root \n" >> $executable
printf "rm -r -f * \n" >> $executable

n_jobs=`ls ${current_dir}/Samples_JPsi_${version} | wc -l`
echo $n_jobs

printf "executable       = $executable \n" >> $submitfile
printf "arguments        = \$(ProcId) \n" >> $submitfile
printf "# log              = log/log_\$(ClusterId).\$(ProcId).log \n" >> $submitfile
printf "# output           = output/output_\$(ClusterId).\$(ProcId).out \n" >> $submitfile
printf "# error            = error/error_\$(ClusterId).\$(ProcId).err \n" >> $submitfile
printf "universe         = vanilla \n" >> $submitfile
printf "should_transfer_files   = Yes \n" >> $submitfile
printf "+JobFlavour      = \"longlunch\" \n" >> $submitfile
printf "queue $n_jobs \n" >> $submitfile

chmod 777 $executable
