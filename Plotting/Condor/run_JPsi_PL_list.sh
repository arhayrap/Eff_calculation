#!/bin/bash 
cp /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_13/src/Eff_calculation/Plotting/Condor/../Grid.C . 
cp /afs/cern.ch/user/a/arhayrap/scratch0/CMSSW_14_0_13/src/Eff_calculation/Plotting/Condor/../JPsi_PL_list.txt . 
index=$1 
hadd input.root `head -$((index*5)) JPsi_PL_list.txt | tail -5` 
root -b -q Grid.C 
cp output.root /eos/user/a/arhayrap/BParking/TriggerScaleFactors/JPsi_PL_list//output_$1.root 
rm -r -f * 
