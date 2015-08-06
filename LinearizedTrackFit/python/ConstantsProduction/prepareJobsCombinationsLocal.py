import os
import subprocess
import time

__author__ = 'demattia'

# The executable used to run training and testing
executable = "/Users/demattia/RemoteProjects/RemoteProjects"

# queue
jobs_queue = []


def array_to_string(a, separator):
    s_a = ""
    for x in a:
        s_a += str(x)+separator
    return s_a.rstrip(separator)


# Function to generate the configuration file for a given combination
def prepare_job(layers, input_file_name, radius_cut_min, radius_cut_max,
                fit_type, oneOverPt_min, oneOverPt_max,
                stub_coordinates, track_parameters, pre_estimate):
    job_dir = "Combinations_"+fit_type+"/"
    pre_estimate_file_transverse = ""
    pre_estimate_file_longitudinal = ""
    if pre_estimate:
        job_dir = "PreEstimate_"+fit_type+"/"
    else:
        current_dir = os.getcwd()
        pre_estimate_file_transverse = current_dir+"/PreEstimate_Transverse/"
        pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal_Rz/"
        if fit_type.find("Longitudinal") != -1 and fit_type.find("Longitudinal_Rz") == -1:
            pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal/"
        # if fit_type.find("Longitudinal_Rz") != -1 or fit_type.find("ExtrapolatedR") != -1 or fit_type.find("_phiRz") != -1:
        #     pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal_Rz/"
        if fit_type.find("Transverse_Pz") != -1:
            pre_estimate_file_transverse = current_dir+"/PreEstimate_Transverse_Pz/"
            pre_estimate_file_longitudinal = current_dir+"/PreEstimate_Longitudinal_Rz/"

    os.system("mkdir -p "+job_dir)

    configuration_file = open(job_dir+"conf_train.txt", "w")

    train_start = 0.
    train_end = 0.5
    test_start = 0.8
    test_end = 1.

    configuration_file.write(str(0)+"\n")
    configuration_file.write(input_file_name+"\n")
    configuration_file.write(array_to_string(layers, " ")+"\n")
    configuration_file.write(str(train_start)+"\n")
    configuration_file.write(str(train_end)+"\n")
    configuration_file.write(str(test_start)+"\n")
    configuration_file.write(str(test_end)+"\n")
    configuration_file.write(array_to_string(stub_coordinates, " ")+"\n")
    configuration_file.write(array_to_string(track_parameters, " ")+"\n")
    configuration_file.write(array_to_string(radius_cut_min, " ")+"\n")
    configuration_file.write(array_to_string(radius_cut_max, " ")+"\n")
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

    # print os.getcwd()+"/"+job_dir+"jobFile.sh"
    # print "chmod +x "+os.getcwd()+"/"+job_dir+"jobFile.sh"
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
file_full = files_dir+"extracted_fullTracker_bigger.root"

layers_list = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
radius_cut_min = [0.]*len(layers_list)
radius_cut_max = [1000.]*len(layers_list)


class JobType:
    def __init__(self, fit_type, stub_coordinates, track_parameters):
        self.fit_type = fit_type
        self.stub_coordinates = stub_coordinates
        self.track_parameters = track_parameters


def prepare_all_jobs(job_types, pre_estimate, pt_min=2., pt_max=0.):
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
        prepare_job(layers_list, file_full, radius_cut_min, radius_cut_max,
                    j.fit_type, oneOverPtMin, oneOverPtMax, j.stub_coordinates, j.track_parameters, pre_estimate)


# # Prepare and submit jobs for the pre-estimates of pt (or pz) and cot(theta)
# pre_estimate = True
# job_types_pre = []
# # job_types_pre.append(JobType("Transverse", ["phi"], ["charge/pt"]))
# # job_types_pre.append(JobType("Transverse_PhiR", ["phi", "R"], ["charge/pt"]))
# job_types_pre.append(JobType("Transverse_PhiZ", ["phi", "z"], ["charge/pt"]))
# # job_types_pre.append(JobType("Transverse_Pz", ["phi"], ["chargeOverPz"]))
# job_types_pre.append(JobType("Longitudinal", ["z"], ["cotTheta"]))
# job_types_pre.append(JobType("Longitudinal_Rz", ["R", "z"], ["cotTheta"]))
# job_types_pre.append(JobType("Longitudinal_Rz_tgTheta", ["R", "z"], ["tgTheta"]))
# Full pT range
# prepare_all_jobs(job_types_pre, pre_estimate)
# # Low pT only
# # prepare_all_jobs(job_types_pre, pre_estimate, pt_min=2., pt_max=10.)


# # # Prepare and submit jobs for the second PCA
pre_estimate = False
# job_types = []
# # job_types.append(JobType("Transverse", ["CorrectedPhiFirstOrder"], ["charge/pt", "phi"]))
# # job_types.append(JobType("Transverse_SecondOrder", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# # # job_types.append(JobType("Transverse_Pz", ["CorrectedPhiFirstOrderPz"], ["charge/pt", "phi"]))
# job_types.append(JobType("Longitudinal", ["CorrectedZFirstOrder"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_Rz", ["CorrectedZFirstOrder"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_SecondOrder", ["CorrectedZSecondOrder"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_Rz_SecondOrder", ["CorrectedZSecondOrder"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_SecondOrder_GEN", ["CorrectedZSecondOrderGen"], ["cotTheta", "z0"]))
# job_types.append(JobType("Longitudinal_Exact_GEN", ["CorrectedZExactGen"], ["cotTheta", "z0"]))
# prepare_all_jobs(job_types, pre_estimate)

# Transverse plane only and splitting low and high pT
# 2 < pT < 10 GeV/c
job_types = []
# job_types.append(JobType("Transverse_2_10", ["CorrectedPhiFirstOrder"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_2_10", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_2_10", ["CorrectedPhiFirstOrderPz"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_2_10", ["CorrectedPhiSecondOrderGen"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Exact_GEN_2_10", ["CorrectedPhiExactGen"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Exact_GEN_ExactR_2_10", ["CorrectedPhiExactGenExactR"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_DeltaZ_2_10", ["CorrectedPhiSecondOrderGenDeltaZ"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_ExactR_2_10", ["CorrectedPhiSecondOrderGenExactR"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedR_2_10", ["CorrectedPhiSecondOrderExtrapolatedR"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_phiRz_2_10",
                         ["CorrectedPhiSecondOrder", "R", "CorrectedZSecondOrder"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrder_2_10",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrder"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_GenExactRSecondOrderNonRadialStripCorrection_2_10",
                         ["CorrectedPhiSecondOrderGenExactRNonRadialStripCorrection"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrection_2_10",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection"],
                         ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrection_phiExtrapolatedRz_2_10",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection",
                          "ExtrapolatedRSecondOrder",
                          "CorrectedZSecondOrder"],
                         ["charge/pt", "phi"]))
prepare_all_jobs(job_types, pre_estimate, pt_min=2., pt_max=10.)

# pT > 10 GeV/c
job_types = []
# job_types.append(JobType("Transverse_10_more", ["CorrectedPhiFirstOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_10_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Pz_10_more", ["CorrectedPhiFirstOrderPz"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_10_more", ["CorrectedPhiSecondOrderGen"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Exact_GEN_10_more", ["CorrectedPhiExactGen"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_Exact_GEN_ExactR_10_more", ["CorrectedPhiExactGenExactR"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_DeltaZ_10_more", ["CorrectedPhiSecondOrderGenDeltaZ"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_GEN_ExactR_10_more", ["CorrectedPhiSecondOrderGenExactR"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedR_10_more", ["CorrectedPhiSecondOrderExtrapolatedR"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_phiRz_10_more",
#                          ["CorrectedPhiSecondOrder", "R", "CorrectedZSecondOrder"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrder_10_more",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrder"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_GenExactRSecondOrderNonRadialStripCorrection_10_more",
                         ["CorrectedPhiSecondOrderGenExactRNonRadialStripCorrection"], ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrection_10_more",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection"],
                         ["charge/pt", "phi"]))
job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrection_phiExtrapolatedRz_10_more",
                         ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection",
                          "ExtrapolatedRSecondOrder",
                          "CorrectedZSecondOrder"],
                         ["charge/pt", "phi"]))
prepare_all_jobs(job_types, pre_estimate, pt_min=10.)


# job_types = []
# job_types.append(JobType("Transverse_SecondOrder_10_12", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedR_10_12", ["CorrectedPhiSecondOrderExtrapolatedR"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, pre_estimate, pt_min=10., pt_max=12)

# job_types = []
# job_types.append(JobType("Transverse_SecondOrder_15_20", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedR_15_20", ["CorrectedPhiSecondOrderExtrapolatedR"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, pre_estimate, pt_min=15., pt_max=20)

# job_types = []
# job_types.append(JobType("Transverse_SecondOrder_20_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# job_types.append(JobType("Transverse_SecondOrder_ExtrapolatedR_20_more", ["CorrectedPhiSecondOrderExtrapolatedR"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, pre_estimate, pt_min=20.)


# job_types = []
# job_types.append(JobType("Transverse_SecondOrder_phiRz_20_more", ["CorrectedPhiSecondOrder", "R", "CorrectedZSecondOrder"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, pre_estimate, pt_min=20.)
#
# job_types = []
# job_types.append(JobType("Transverse_SecondOrder_phiRz_15_20", ["CorrectedPhiSecondOrder", "R", "CorrectedZSecondOrder"], ["charge/pt", "phi"]))
# prepare_all_jobs(job_types, pre_estimate, pt_min=15., pt_max=20.)


# # # 2 < pT < 15 GeV/c
# # job_types = []
# # job_types.append(JobType("Transverse_2_15", ["CorrectedPhi"], ["charge/pt", "phi"]))
# # job_types.append(JobType("Transverse_SecondOrder_2_15", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# # job_types.append(JobType("Transverse_Pz_2_15", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
# # prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 15.)
#
# # # pT > 15 GeV/c
# # job_types = []
# # job_types.append(JobType("Transverse_15_more", ["CorrectedPhi"], ["charge/pt", "phi"]))
# # job_types.append(JobType("Transverse_SecondOrder_15_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"]))
# # job_types.append(JobType("Transverse_Pz_15_more", ["CorrectedPhiPz"], ["charge/pt", "phi"]))
# # prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 15.)


def check_pid(pid):
    """ Check For the existence of a unix pid. """
    try:
        # Kill with value 0 does not kill the process
        os.kill(pid, 0)
        return True
    except OSError:
        return False


# Now submit the jobs to keep the total number of jobs running in parallel < N
print jobs_queue
maximum_parallel_jobs = 2

failed_processes = []
running_processes = []

while len(jobs_queue) > 0 or len(running_processes) > 0:
    print "jobs left", len(jobs_queue)

    # Get the list of all running process ids
    failed_processes.extend([p for p in running_processes if p[1].poll() is not None and p[1].poll() != 0])
    running_processes = [p for p in running_processes if p[1].poll() is None]
    print "running jobs:", [p[0] for p in running_processes]

    for i in range(min(maximum_parallel_jobs - len(running_processes), len(jobs_queue))):
        job_command = jobs_queue.pop()
        process = subprocess.Popen(job_command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        running_processes.append([job_command.lstrip(os.getcwd()).rstrip("jobFile.sh"), process])

    if len(jobs_queue) == 0 and len(running_processes) == 0:
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
