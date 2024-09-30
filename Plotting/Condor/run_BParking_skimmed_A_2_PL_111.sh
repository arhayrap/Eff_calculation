#!/bin/bash 
cp /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_13/src/Eff_calculation/Plotting/Condor/../Grid.C . 
cp /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_13/src/Eff_calculation/Plotting/Condor/../BParking_skimmed_A_2_PL_111.txt . 
index=$1 
hadd input.root `head -$((index*5)) BParking_skimmed_A_2_PL_111.txt | tail -5` 
root -b -q Grid.C 
cp output.root /eos/user/a/arhayrap/BParking/TriggerScaleFactors/BParking_skimmed_A_2_PL_111//output_$1.root 
rm -r -f * 
