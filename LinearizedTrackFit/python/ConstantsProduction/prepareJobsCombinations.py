import os
import subprocess
import time

__author__ = 'demattia'


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


def array_to_string(a):
    s_a = ""
    for x in a:
        s_a += str(x)+","
    return s_a.rstrip(",")


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

    configuration_file = open(job_dir+"BuildLinearizedTrackFitMatrix_cfg.py", "w")
    for line in open("BuildLinearizedTrackFitMatrix_template_cfg.py"):
        line = line.replace("INPUT_FILE_NAME", c.input_file_name)
        line = line.replace("ONE_OVER_PT_MIN", str(oneOverPt_min))
        line = line.replace("ONE_OVER_PT_MAX", str(oneOverPt_max))
        line = line.replace("STUB_COORDINATES", stub_coordinates)
        line = line.replace("TRACK_PARAMETERS", track_parameters)
        line = line.replace("CHARGE_REGIONS", "1")
        line = line.replace("LAYERS_LIST", array_to_string(c.layers))
        line = line.replace("RADIUS_CUTS_MIN", array_to_string(c.radius_cut_min))
        line = line.replace("RADIUS_CUTS_MAX", array_to_string(c.radius_cut_max))
        line = line.replace("FIRST_ORDER_CHARGEOVERPT_COEFFICIENTS_FILE_NAME", pre_estimate_file_transverse)
        line = line.replace("FIRST_ORDER_COTTHETA_COEFFICIENTS_FILE_NAME", pre_estimate_file_longitudinal)
        configuration_file.write(line)
    configuration_file.close()
    slurm_file = open(job_dir+"jobFile.slrm", "w")
    job_name = c.name.split("_")[1]
    if len(c.name.split("_")) > 3:
        job_name += "_"+c.name.split("_")[3]
    for line in open("jobFile_template.slrm"):
        line = line.replace("WORKING_DIR", os.getcwd()+"/"+job_dir)
        line = line.replace("PCA", job_name)
        slurm_file.write(line)
    slurm_file.close()
    # print "Submitting job for", c.name
    # print "cd "+job_dir+"; sbatch jobFile.slrm; cd -"
    jobs_queue.append("cd "+job_dir+"; sbatch jobFile.slrm; cd -")
    # os.system("cd "+job_dir+"; sbatch jobFile.slrm; cd -")


# Function to make the combinations unique. If there is a duplicate take the one with the biggest region index since
# this ensures that radius cuts are applied.
# We split this step for now. We instead run all the combinations and use the duplicates to check that their results
# are the same.


# It has to be able to run with a pre- or post-fix to derive the pre-estimate coefficients. When running
# the second PCA it will look for the pre-estimates in the directories with the same name of the current
# ones but with the added pre- or post-fix.
# The layers and disks must be given in order from the smallest index to the biggest index.


files_dir = "/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/"
file_barrel = files_dir+"Muons_Pt_2_more_Eta_m01p09_Phi_m05p1_3_z0_m15p15/JobFiles/extracted.root"
file_hybrid = files_dir+"Muons_Pt_2_more_Eta_0p7_1p7_Phi_m05p1_3_z0_m15p15_HYBRID/JobFiles/extracted.root"
# file_endcaps = files_dir+"Muons_Pt_2_more_Eta_1p5_2p5_Phi_m05p1_3_z0_m15p15_ENDCAPS/JobFiles/extracted.root"
file_endcaps = files_dir+"Endcaps/A/extracted.root"

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
radius_cut_min[12-5] = 61.
radius_cut_max[11-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_6"
# layers_list = [5, 6, 11, 13, 14, 15]
layers_list = [5, 6, 11, 12, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[13-5] = 61.
radius_cut_max[12-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_7"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[14-5] = 61.
radius_cut_max[13-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_8"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[15-5] = 61.
radius_cut_max[14-5] = 61.
# combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))
combinations.append(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_9"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
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
            # Generate the jobFile for the full region and the configuration files for each region
            slurm_file_name = "jobFile_"+region[0].name+".slrm"
            slurm_file = open(job_dir+slurm_file_name, "w")
            for line in open("jobFile_region_template.slrm"):
                line = line.replace("REGION", region[0].name)
                slurm_file.write(line)
            slurm_file.write("\n")
            for c in region:
                print c.name
                prepare_job(c, j.fit_type, oneOverPtMin, oneOverPtMax, j.stub_coordinates, j.track_parameters, pre_estimate)
                slurm_file.write("cd "+os.getcwd()+"/"+job_dir+"/"+c.name+"/ \n")
                slurm_file.write("cmsRun BuildLinearizedTrackFitMatrix_cfg.py \n")
            slurm_file.write("\n")
            slurm_file.write('echo "Ended at `date` on `hostname`" \n')
            slurm_file.write("exit 0 \n")
            slurm_file.close()
            # sbatch_command = "cd "+job_dir+"; sbatch "+slurm_file_name+"; cd -"
            # jobs_queue.append(sbatch_command)
            # os.system(sbatch_command)


# # Prepare and submit jobs for the pre-estimates of pt (or pz) and cot(theta)
# pre_estimate = True
# job_types_pre = [JobType("Transverse", '"phi"', '"charge/pt"')]
# job_types_pre.append(JobType("Transverse_Pz", '"phi"', '"chargeOverPz"'))
# job_types_pre.append(JobType("Longitudinal", '"z"', '"cotTheta"'))
# job_types_pre.append(JobType("Longitudinal_Rz", '"R", "z"', '"cotTheta"'))
# prepare_all_jobs(job_types_pre, combinations, pre_estimate)

# Prepare and submit jobs for the second PCA
pre_estimate = False
job_types = []
# job_types.append(JobType("Transverse", '"CorrectedPhi"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_SecondOrder", '"CorrectedPhiSecondOrder"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_Pz", '"CorrectedPhiPz"', '"charge/pt", "phi"'))
job_types.append(JobType("Longitudinal", '"CorrectedZ"', '"cotTheta", "z0"'))
# job_types.append(JobType("Longitudinal_Rz", '"CorrectedZ"', '"cotTheta", "z0"'))
# job_types.append(JobType("Longitudinal_SecondOrder", '"CorrectedZSecondOrder"', '"cotTheta", "z0"'))
job_types.append(JobType("Longitudinal_Rz_SecondOrder", '"CorrectedZSecondOrder"', '"cotTheta", "z0"'))
prepare_all_jobs(job_types, combinations, pre_estimate)

# # Transverse plane only and splitting low and high pT
# 
# # 2 < pT < 10 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_2_10", '"CorrectedPhi"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_SecondOrder_2_10", '"CorrectedPhiSecondOrder"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_Pz_2_10", '"CorrectedPhiPz"', '"charge/pt", "phi"'))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 10.)
# 
# # pT > 10 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_10_more", '"CorrectedPhi"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_SecondOrder_10_more", '"CorrectedPhiSecondOrder"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_Pz_10_more", '"CorrectedPhiPz"', '"charge/pt", "phi"'))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 10.)
# 
# # 2 < pT < 15 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_2_15", '"CorrectedPhi"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_SecondOrder_2_15", '"CorrectedPhiSecondOrder"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_Pz_2_15", '"CorrectedPhiPz"', '"charge/pt", "phi"'))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 2., pt_max = 15.)
# 
# # pT > 15 GeV/c
# job_types = []
# job_types.append(JobType("Transverse_15_more", '"CorrectedPhi"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_SecondOrder_15_more", '"CorrectedPhiSecondOrder"', '"charge/pt", "phi"'))
# job_types.append(JobType("Transverse_Pz_15_more", '"CorrectedPhiPz"', '"charge/pt", "phi"'))
# prepare_all_jobs(job_types, combinations, pre_estimate, pt_min = 15.)



# Now submit the jobs to keep the total number of jobs running in parallel < N
# print jobs_queue
maximum_parallel_jobs = 30

while len(jobs_queue) > 0:
    output = subprocess.check_output("squeue -u demattia", shell=True)
    number_of_jobs = 0
    for line in output.split("\n"):
        splitted_line = line.split()
        if len(splitted_line) > 0 and splitted_line[0] != "JOBID":
            # if splitted_line[2] == "PCA":
            number_of_jobs += 1

    for i in range(min(maximum_parallel_jobs - number_of_jobs, len(jobs_queue))):
        job_command = jobs_queue.pop()
        print "submitting job:"
        print job_command
        os.system(job_command)
    if len(jobs_queue) == 0:
        break
    # Every 10 seconds
    time.sleep(10)
