#!/bin/tcsh

setenv current_DIr `pwd` > output_log.txt
setenv X509_USER_PROXY /afs/cern.ch/user/a/arhayrap/public/x509up_u123840
#cd /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/
#source /afs/cern.ch/cms/cmsset_default.csh

cmsrel CMSSW_14_0_5
cd ./CMSSW_14_0_5/src
cp ~/public/x509up_u123840 /tmp/

cmsenv

# cd $current_DIr

set newid=`expr $1 + 1`

cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/PreSelection .
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/Condor/BParking/Samples_A_1_UL ./sample_names
cp /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/Trig_ScaleFactors_OnBParking_1.py .

scram b

set filelist=( ./sample_names/* )
echo ${filelist[$newid]}

mv ${filelist[$newid]} BToKPhi.txt

set filenames = ""

foreach line ( "`cat ./BToKPhi.txt`" )
    set filenames = "$filenames $line"
end

cmsRun Trig_ScaleFactors_OnBParking_1.py $filenames

cp Trigger_ScaleFactors_OnBParking.root /eos/user/a/arhayrap/BParking_skimmed_A_1_UL_111/output_file_$1.root

# sleep 20m

rm -r *
