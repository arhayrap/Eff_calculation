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
# cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/Condor/JPsi/Samples_JPsi_2018_Autumn ./sample_names
cp -r -f /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/Condor/JPsi/Samples_JPsi ./sample_names
cp /afs/cern.ch/user/a/arhayrap/scratch0/Eff_Calculation/CMSSW_14_0_5/src/Trig_ScaleFactors_OnMC.py .

scram b

set filelist=( ./sample_names/* )
echo ${filelist[$newid]}

mv ${filelist[$newid]} BToKPhi.txt

foreach line ( "`cat ./BToKPhi.txt`" )

    echo "xrdcp -r -f " $line " ."> download.sh

end

chmod 777 download.sh

./download.sh

set old_file_path = `find ./*.root -type f | head -n 1`

mv $old_file_path ./input_sample.root

echo $current_DIr/sample.root > sample.txt
cmsRun Trig_ScaleFactors_OnMC.py > output_log.txt

# cp output_log.txt /eos/user/a/arhayrap/JPsi_output_data/output_log_$1.txt
# cp -r * /eos/user/a/arhayrap/ALL_CONDOR_FILES
# cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_output_data_UL_Tight_Loose_ID_test_old_version_222/Trigger_ScaleFactors_OnMC_$1.root
# cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_output_data_PL_Tight_Loose_ID_test_old_version_444/Trigger_ScaleFactors_OnMC_$1.root
cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_output_data_UL_Tight_Loose_ID_test_old_version_555/Trigger_ScaleFactors_OnMC_$1.root
# cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_output_data_UL_Tight_Loose_ID_test_old_version_444/Trigger_ScaleFactors_OnMC_$1.root
# cp Trigger_ScaleFactors_OnMC.root /eos/user/a/arhayrap/JPsi_output_data_2018Autumn_Tight_Loose_ID_1/Trigger_ScaleFactors_OnMC_$1.root

# sleep 20m

rm -r *

