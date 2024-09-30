#!/bin/bash

list=$1
base="${PWD}"

if [[ $list == "" ]];
then
    echo "Give the project name."
    exit
fi

if ! test -f $base/../$list;
then
    echo "A list of data files with the name of '$list' does not exist."
    exit
fi

nFiles_per_job=5
nFiles=`cat ../$list | wc -l`
nJobs=$((nFiles/nFiles_per_job))
project_name=`echo $list | grep -o '^[^\.]*'`
username=`whoami`
first_letter=`echo $username | head -c 1 | tail -c 1`
echo $project_name

rm run_$project_name.sh Run_$project_name.sub hadd_$project_name.sh

output="/eos/user/$first_letter/$username/BParking/TriggerScaleFactors/$project_name/"
echo $output
mkdir -p $output

printf "Executable = run_$project_name.sh \n" >> Run_$project_name.sub
printf "universe = vanilla \n" >> Run_$project_name.sub
printf "Should_Transfer_Files = YES \n" >> Run_$project_name.sub
printf "arguments = \$(ProcId) \n" >> Run_$project_name.sub
# printf "Output = logs/message_\$(Cluster)_\$(Process).out \n" >> Run_$project_name.sub
printf "Error = logs/message_\$(Cluster)_\$(Process).err \n" >> Run_$project_name.sub
# printf "Log = logs/message_\$(Cluster)_\$(Process).log \n" >> Run_$project_name.sub
printf "queue $nJobs \n" >> Run_$project_name.sub

printf "#!/bin/bash \n" >> run_$project_name.sh
printf "cp $base/../Grid.C . \n" >> run_$project_name.sh
printf "cp $base/../$list . \n" >> run_$project_name.sh
printf "index=\$1 \n" >> run_$project_name.sh
printf "hadd input.root \`head -\$((index*$nFiles_per_job)) $list | tail -$nFiles_per_job\` \n" >> run_$project_name.sh
printf "root -b -q Grid.C \n" >> run_$project_name.sh
printf "cp output.root $output/output_\$1.root \n" >> run_$project_name.sh
printf "rm -r -f * \n" >> run_$project_name.sh

printf "hadd $project_name.root $output/*.root \n" >> hadd_$project_name.sh

chmod +777 run_$project_name.sh
chmod +777 hadd_$project_name.sh

mkdir -p logs
