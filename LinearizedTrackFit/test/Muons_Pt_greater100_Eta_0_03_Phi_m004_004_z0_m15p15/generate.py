import os

max_events = 1000
total_jobs = 1000

os.system("mkdir -p JobFiles")

# Starts from 1 to match the SLURM array index
for i in range(1, total_jobs+1):
    template_file = open("oneStep_wideRange_template.py", "r")
    output_file = open("JobFiles/oneStep_wideRange_"+str(i)+".py", "w")
    for line in template_file:
        line = line.replace("GENERATORSEED", str(i))
        line = line.replace("VTXSMEAREDSEED", str(i+1))
        line = line.replace("G4SIMHITSSEED", str(i+2))
        line = line.replace("MIXSEED", str(i+3))
        line = line.replace("MAXEVENTS", str(max_events))
        line = line.replace("INDEX", str(i))
        output_file.write(line)

    output_file.close()

slurm_template_file = open("jobFile_template.slrm", "r")
slurm_output_file = open("JobFiles/jobFile.slrm", "w")
for line in slurm_template_file:
    line = line.replace("ARRAYMAX", str(total_jobs))
    line = line.replace("JOBSDIR", os.getcwd()+"/JobFiles")
    slurm_output_file.write(line)
slurm_output_file.close()
