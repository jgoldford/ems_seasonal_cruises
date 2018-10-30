import subprocess
import glob


# make a function that 
def make_run_script(pipelinePath,run_dir,sample_dir,sample_file,learn_file,working_dir,fTrunc):
	fname = run_dir + '/run_script.sh';
	base = "Rscript " + pipelinePath;
	line = base + ' -w ' + working_dir + " -i " + sample_dir + " -s " + sample_file+  " -l " + learn_file + " -f " + str(fTrunc);
	with open(fname,'w') as f:
		f.write("#!/bin/bash -l")
		f.write("\n")
		f.write("module load R/3.3.0")
		f.write("\n")
		f.write(line)
		f.write("\n")
	return fname

def run_qsub(fname,time):
	call = 'chmod 777 ' + fname;
	subprocess.call(call,shell=True);

	call = 'qsub -l h_rt=' + time + ' -P bioinfor ' + fname;
	subprocess.call(call,shell=True);


# define samples directory
sample_dir = "/projectnb/microcrm/marine/flash_final";
base_dir = "/projectnb/microcrm/marine/dada2_merged";


# reconstruct sample and learning files:
sample_file = '/projectnb/microcrm/marine/dada2_merged/samples.samples';
learn_file = '/projectnb/microcrm/marine/dada2_merged/samples.learn';

# the working directories (where all filtered fastq files go, are named by the sample ids)
working_dir =  '/projectnb/microcrm/marine/dada2_merged/batch1';
# make the dirs
subprocess.call('mkdir ' + working_dir,shell=True);

# the run directories 
run_dir =  working_dir + '/run';
# make directories for runs on each sub_sample
subprocess.call('mkdir ' + run_dir,shell=True);

# path to the R dada2 pipeline (one R script)
pipelinePath = "/projectnb/microcrm/marine/dada2_merged/dada2pipelineForMergedSequences.R";

# set parameters for 16S region (paired end reads are ~ 409 nucleoties long, so there are )
run_script = make_run_script(pipelinePath,run_dir,sample_dir,sample_file,learn_file,working_dir,400);
run_qsub(run_script,"48:00:00");
