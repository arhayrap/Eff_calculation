#!/bin/tcsh 

setenv current_DIr /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/JPsi 
setenv X509_USER_PROXY ~/public/x509up_u123840 
# source /afs/cern.ch/cms/cmsset_default.csh 
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cmsrel CMSSW_14_0_5 
cd ./CMSSW_14_0_5/src 
cp ~/public/x509up_u123840 /tmp/ 
cmsenv 
set newid=1 
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/JPsi/../../PreSelection . 
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/JPsi/Samples_JPsi_PL ./sample_names 
cp       /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/JPsi/../../Trig_ScaleFactors_OnMC.py . 
scram b 
set filelist=( ./sample_names/* ) 
echo ${filelist[$newid]} 
mv ${filelist[$newid]} BToKPhi.txt 
set filenames = "" 
foreach line ( "`cat ./BToKPhi.txt`" ) 
    set filenames = "$filenames $line" 
end 
cmsRun Trig_ScaleFactors_OnMC.py $filenames 
cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_PL//output_file_$1.root 
rm -r -f * 
