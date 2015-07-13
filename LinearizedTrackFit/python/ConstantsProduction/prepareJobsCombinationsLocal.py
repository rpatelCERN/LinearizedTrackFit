import os
import subprocess
import time

__author__ = 'demattia'

# The executable used to run training and testing
executable = "/Users/demattia/RemoteProjects/RemoteProjects"

# queue
jobs_queue = []


# Class to store all the needed information for a combination
class Combination:
    def __init__(self, name, layers_list, radius_cut_min, radius_cut_max, input_file_name):
        self.layers = layers_list
        self.name = name
        self.radius_cut_min = radius_cut_min
        self.radius_cut_max = radius_cut_max
        self.input_file_name = input_file_name


# Function to generate all combinations from a given set of 6 layers or disks
# No repetition of the same layer is allowed in the input layers_list
def generate_layers_combinations(base_name, layers_list, radius_cut_min, radius_cut_max, input_file_name):
    combinations = [Combination(base_name+"_All", layers_list, radius_cut_min, radius_cut_max, input_file_name)]
    for layer in layers_list:
        name = base_name+"_Removed_" + str(layer)
        layers_combination = [x for x in layers_list if x != layer]
        combinations.append(Combination(name, layers_combination, radius_cut_min, radius_cut_max, input_file_name))
    return combinations


def array_to_string(a, separator):
    s_a = ""
    for x in a:
        s_a += str(x)+separator
    return s_a.rstrip(separator)


# Function to generate the configuration file for a given combination
def prepare_job(c, fit_type, oneOverPt_min, oneOverPt_max, stub_coordinates, track_parameters, pre_estimate):
    print "preparing files for", c.name
    job_dir = "Combinations_"+fit_type+"/"+c.name+"/"
    pre_estimate_file_transverse = ""
    pre_estimate_file_longitudinal = ""
    if pre_estimate:
        job_dir = "PreEstimate_"+fit_type+"/"+c.name+"/"
    else:
        current_dir = os.getcwd()
        pre_estimate_file_transverse = current_dir+"/PreEstimate_Transverse/"+c.name+"/matrixVD_0_pre_chargeOverPt.txt"
        pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal/"+c.name+"/matrixVD_0_pre_cotTheta.txt"
        if fit_type.find("Transverse_Pz") != -1:
            pre_estimate_file_transverse = current_dir+"/PreEstimate_Transverse_Pz/"+c.name+"/matrixVD_0_pre_chargeOverPz.txt"
        if fit_type.find("Longitudinal_Rz") != -1:
            pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal_Rz/"+c.name+"/matrixVD_0_pre_cotTheta.txt"

    os.system("mkdir -p "+job_dir)

    configuration_file = open(job_dir+"conf_train.txt", "w")

    trainStart = 0.
    trainEnd = 0.2
    testStart = 0.8
    testEnd = 1.

    configuration_file.write(str(0)+"\n")
    configuration_file.write(c.input_file_name+"\n")
    configuration_file.write(array_to_string(c.layers, " ")+"\n")
    configuration_file.write(str(trainStart)+"\n")
    configuration_file.write(str(trainEnd)+"\n")
    configuration_file.write(str(testStart)+"\n")
    configuration_file.write(str(testEnd)+"\n")
    configuration_file.write(array_to_string(stub_coordinates, " ")+"\n")
    configuration_file.write(array_to_string(track_parameters, " ")+"\n")
    configuration_file.write(array_to_string(c.radius_cut_min, " ")+"\n")
    configuration_file.write(array_to_string(c.radius_cut_max, " ")+"\n")
    configuration_file.write(pre_estimate_file_transverse+"\n")
    configuration_file.write(pre_estimate_file_longitudinal+"\n")
    configuration_file.write(str(oneOverPt_min)+"\n")
    configuration_file.write(str(oneOverPt_max)+"\n")
    configuration_file.write(str(1)+"\n")
    configuration_file.write(str(0.)+"\n")
    configuration_file.write(str(0.8)+"\n")
    configuration_file.write(str(1)+"\n")
    configuration_file.write(str(-3.)+"\n")
    configuration_file.write(str(3.)+"\n")
    configuration_file.write(str(1)+"\n")
    configuration_file.write(str(-15.)+"\n")
    configuration_file.write(str(15.)+"\n")
    configuration_file.write(str(1)+"\n")
    configuration_file.write(str(1)+"\n")
    configuration_file.write(str(0)+"\n")

    configuration_file.close()

    # Generate the testing file where the only difference is the first line
    os.system("sed '1s/.*/1/' "+os.getcwd()+"/"+job_dir+"conf_train.txt > "+os.getcwd()+"/"+job_dir+"conf_test.txt")


    job_file = open(job_dir+"jobFile.sh", "w")
    job_file.write("#!/bin/sh\n")
    job_file.write("cd "+os.getcwd()+"/"+job_dir+"\n")
    job_file.write(executable+" conf_train.txt &> log_train.txt\n")
    job_file.write(executable+" conf_test.txt &> log_test.txt\n")
    job_file.close()

    # print "Submitting job for", c.name
    # print "cd "+job_dir+"; sbatch jobFile.slrm; cd -"
    # jobs_queue.append("cd "+job_dir+"; sbatch jobFile.slrm; cd -")
    # os.system("cd "+job_dir+"; sbatch jobFile.slrm; cd -")
    jobs_queue.append(os.getcwd()+"/"+job_dir+"jobFile.sh")
    os.system("chmod +x "+os.getcwd()+"/"+job_dir+"jobFile.sh")


# Function to make the combinations unique. If there is a duplicate take the one with the biggest region index since
# this ensures that radius cuts are applied.
# We split this step for now. We instead run all the combinations and use the duplicates to check that their results
# are the same.


# It has to be able to run with a pre- or post-fix to derive the pre-estimate coefficients. When running
# the second PCA it will look for the pre-estimates in the directories with the same name of the current
# ones but with the added pre- or post-fix.
# The layers and disks must be given in order from the smallest index to the biggest index.


files_dir = "/Users/demattia/RemoteProjects/"
file_barrel = files_dir+"extracted_prompt_extrapolated.root"
file_hybrid = files_dir+"extracted_hybrid.root"
file_endcaps = files_dir+"extracted_endcaps.root"

combinations = []

name = "Region_1"
layers_list = [5, 6, 7, 8, 9, 10]
# combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_barrel))
combinations.append(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_barrel))

name = "Region_2"
layers_list = [5, 6, 7, 8, 9, 11]
# combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))
combinations.append(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))

name = "Region_3"
layers_list = [5, 6, 7, 8, 11, 12]
# combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))
combinations.append(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))

name = "Region_4"
layers_list = [5, 6, 7, 11, 12, 13]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[11-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_hybrid))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_hybrid))

name = "Region_5"
layers_list = [5, 6, 11, 12, 13, 14]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[13-5] = 61.
radius_cut_min[12-5] = 61.
radius_cut_max[11-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_6"
# layers_list = [5, 6, 11, 13, 14, 15]
layers_list = [5, 6, 11, 12, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[14-5] = 61.
radius_cut_min[13-5] = 61.
radius_cut_max[12-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_7"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[15-5] = 61.
radius_cut_min[14-5] = 61.
radius_cut_max[13-5] = 61.
radius_cut_max[12-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_8"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[15-5] = 61.
radius_cut_max[14-5] = 61.
radius_cut_max[13-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_9"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_max[14-5] = 61.
radius_cut_max[15-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))


class JobType:
    def __init__(self, fit_type, stub_coordinates, track_parameters):
        self.fit_type = fit_type
        self.stub_coordinates = stub_coordinates
        self.track_parameters = track_parameters


def prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 0.):
    oneOverPtMin = 0.
    oneOverPtMax = 1./pt_min
    if pt_max != 0.:
        oneOverPtMin = 1/pt_max
    backup_dir = "Backup/Backup"
    if pre_estimate:
        backup_dir += "Pre"
    for j in job_types:
        os.system("rm -rf "+backup_dir+"_"+j.fit_type+"_Old; mv "+backup_dir+"_"+j.fit_type+" "+backup_dir+"_"+j.fit_type+"_Old")
        job_dir = "Combinations_"+j.fit_type+"/"
        if pre_estimate:
            job_dir = "PreEstimate_"+j.fit_type+"/"
        os.system("mkdir -p "+backup_dir+"_"+j.fit_type+"; mv "+job_dir+" "+backup_dir)
        os.system("mkdir -p "+job_dir)
        # Use the combinations to modify the cfg, generate directories with the given name, and submit the jobs.
        for region in combinations:
            print ""
            for c in region:
                print c.name
                prepare_job(c, j.fit_type, oneOverPtMin, oneOverPtMax, j.stub_coordinates, j.track_parameters, pre_estimate)


# # Prepare and submit jobs for the pre-estimates of pt (or pz) and cot(theta)
pre_estimate = True
job_types_pre = []
# job_types_pre.append(JobType("Transverse", ["phi"], ["charge/pt"]))
job_types_pre.append(JobType("Transverse_Pz", ["phi"], ["chargeOverPz"]))
job_types_pre.append(JobType("Longitudinal", ["z"], ["cotTheta"]))
# job_types_pre.append(JobType("Longitudinal_Rz", ["R", "z"], ["cotTheta"]))
# Full pT range
# prepare_all_jobs(job_types_pre, combinations, pre_estimate)
# Low pT only
prepare_all_jobs(job_types_pre, combinations, pre_estimate, pt_min = 2., pt_max = 10.)


# Prepare and submit jobs for the second PCA
pre_estimate = False
job_types = []
job_types.append(JobType("Transverse", ["CorrectedPhi"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
# job_types.append(JobType("Longitudinal", ["CorrectedZ"], ["cotTheta", "z0"]))
job_types.append(JobType("Longitudinal_Rz", ["CorrectedZ"], ["cotTheta", "z0"]))
job_types.append(JobType("Longitudinal_SecondOrder", ["CorrectedZSecondOrder"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_Rz_SecondOrder", ["CorrectedZSecondOrder"], ["cotTheta", "z0"]))
prepare_all_jobs(job_types, combinations, pre_estimate)

# Transverse plane only and splitting low and high pT

# 2 < pT < 10 GeV/c
job_types = []
job_types.append(JobType("Transverse_2_10", ["CorrectedPhi"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_2_10", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_2_10", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 10.)

# pT > 10 GeV/c
job_types = []
job_types.append(JobType("Transverse_10_more", ["CorrectedPhi"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_10_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_10_more", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 10.)

# # 2 < pT < 15 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_2_15", ["CorrectedPhi"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_2_15", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_2_15", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 15.)

# # pT > 15 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_15_more", ["CorrectedPhi"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_15_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_15_more", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 15.)


def check_pid(pid):
    """ Check For the existence of a unix pid. """
    try:
        # Kill with value 0 does not kill the process
        os.kill(pid, 0)
        return True
    except OSError:
        return False


# Now submit the jobs to keep the total number of jobs running in parallel < N
print jobs_queue[0]
maximum_parallel_jobs = 8

failed_processes = []
running_processes = []

while len(jobs_queue) > 0:
    print "jobs left", len(jobs_queue)

    # Get the list of all running process ids
    failed_processes.extend([p for p in running_processes if p[1].poll() != None and p[1].poll() != 0])
    running_processes = [p for p in running_processes if p[1].poll() == None]
    print "running jobs:", [p[0] for p in running_processes]

    for i in range(min(maximum_parallel_jobs - len(running_processes), len(jobs_queue))):
        job_command = jobs_queue.pop()
        process = subprocess.Popen(job_command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        running_processes.append([job_command.lstrip(os.getcwd()).rstrip("jobFile.sh"), process])

    if len(jobs_queue) == 0:
        break
    # Every 10 seconds
    time.sleep(10)


print "All jobs finished"
if len(failed_processes) == 0:
    print "No failed jobs"
else:
    print "Failed jobs:"
    for p in failed_processes:
        print p[0]