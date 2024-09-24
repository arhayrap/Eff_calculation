#!/bin/tcsh 

setenv current_DIr /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/BParking 
setenv X509_USER_PROXY ~/public/x509up_u123840 
# source /afs/cern.ch/cms/cmsset_default.csh 
source /cvmfs/cms.cern.ch/cmsset_default.csh 
cmsrel CMSSW_14_0_5 
cd ./CMSSW_14_0_5/src 
cp ~/public/x509up_u123840 /tmp/ 
cmsenv 
set newid=1 
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/BParking/../../PreSelection . 
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/BParking/Samples_A_2_PL ./sample_names 
cp       /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_5/src/Eff_calculation/Condor/BParking/../../Trig_ScaleFactors_OnBParking.py . 
scram b 
set filelist=( ./sample_names/* ) 
echo  
mv  BToKPhi.txt 
set filenames = "" 
foreach line ( "`cat ./BToKPhi.txt`" ) 
    set filenames = "$filenames $line" 
end 
cmsRun Trig_ScaleFactors_OnBParking.py $filenames 
cp Trigger_ScaleFactors_OnBParking.root /eos/user/a/arhayrap/BParking_skimmed_A_2_PL//output_file_$1.root 
rm -r -f * 
