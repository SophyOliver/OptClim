#!/bin/bash
# set the number of nodes and processes per node
#SBATCH --nodes=6

# set the number of tasks (processes) per node.
#SBATCH --ntasks-per-node=16

# set max wallclock time
#SBATCH --time=40:00:00

module load python/anaconda2/5.0.1
source activate myenv

finished="False"
iRun=1 # change this if optimisation process crashed half way last time. At least don't have as 1 if you have previous runs saved, as will restart the whole thing and erase previous runs.
iMax=30 # max number of iterations that can fit in the set wall clock time.

# Loop through optimisation iterations while we're below the set max iterations and the algorithm is not finished
while [ $iRun -le $iMax ] && [ "$finished" == "False" ]; do

	if [ $iRun -eq 1 ]; then
		# Need to restart
		python runOptimise_oneJob.py --restart -dir testloop -v ../Configurations/OxfordMOPS_Configs/TWIN_Configs/OxfordMOPS_dfols_twin_6p.json &>prog_testloop.txt
	
	else
		# Can carry on from previous iteration
		python runOptimise_oneJob.py -dir testloop -v ../Configurations/OxfordMOPS_Configs/TWIN_Configs/OxfordMOPS_dfols_twin_6p.json &>prog_testloop.txt

	fi
	
	# Update counter
	let iRun++
	# Get last line from progress text, which should be either True or False, determining whether the optimisation process is finished or not.
	finished=`tail -1 prog_testloop.txt`
	
	echo The iRun is $iRun
	echo "$finished"
	
done