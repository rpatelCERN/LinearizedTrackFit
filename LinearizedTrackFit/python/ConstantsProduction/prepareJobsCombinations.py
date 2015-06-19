import os

__author__ = 'demattia'


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
    job_dir = "Combinations/"+c.name+"/"
    pre_estimate_file_transverse = ""
    pre_estimate_file_longitudinal = ""
    if pre_estimate:
        job_dir = "PreEstimate_"+fit_type+"/"+c.name+"/"
    else:
        pre_estimate_file_transverse = "PreEstimate_Transverse/"+c.name+"/file.txt"
        pre_estimate_file_longitudinal = "PreEstimate_Longitudinal/"+c.name+"/file.txt"
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
    for line in open("jobFile_template.slrm"):
        line = line.replace("WORKING_DIR", os.getcwd()+"/"+job_dir)
        slurm_file.write(line)
    slurm_file.close()
    print "Submitting job for", c.name
    print "cd "+job_dir+"; sbatch jobFile.slrm; cd -"
    os.system("cd "+job_dir+"; sbatch jobFile.slrm; cd -")


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
file_endcaps = files_dir+"Muons_Pt_2_more_Eta_1p5_2p5_Phi_m05p1_3_z0_m15p15_ENDCAPS/JobFiles/extracted.root"


combinations = []

name = "Region_1"
layers_list = [5, 6, 7, 8, 9, 10]
combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_barrel))

name = "Region_2"
layers_list = [5, 6, 7, 8, 9, 11]
combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))

name = "Region_3"
layers_list = [5, 6, 7, 8, 11, 12]
combinations.extend(generate_layers_combinations(name, layers_list, [0]*16, [1000.]*16, file_hybrid))

name = "Region_4"
layers_list = [5, 6, 7, 11, 12, 13]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[11-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_hybrid))

name = "Region_5"
layers_list = [5, 6, 11, 12, 13, 14]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[12-5] = 61.
radius_cut_max[11-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_6"
layers_list = [5, 6, 11, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[13-5] = 61.
radius_cut_max[12-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_7"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[14-5] = 61.
radius_cut_max[13-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_8"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_min[15-5] = 61.
radius_cut_max[14-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))

name = "Region_9"
layers_list = [5, 11, 12, 13, 14, 15]
radius_cut_min = [0]*16
radius_cut_max = [1000.]*16
radius_cut_max[15-5] = 61.
combinations.extend(generate_layers_combinations(name, layers_list, radius_cut_min, radius_cut_max, file_endcaps))


class JobType:
    def __init__(self, fit_type, stub_coordinates, track_parameters):
        self.fit_type = fit_type
        self.stub_coordinates = stub_coordinates
        self.track_parameters = track_parameters


def prepare_pre_estimate_jobs(job_types, combinations):
    for j in job_types:
        os.system("rm -rf Backup/BackupPre_"+j.fit_type+"_Old; mv Backup/BackupPre_"+j.fit_type+" Backup/BackupPre_"+j.fit_type+"_Old")
        os.system("mkdir -p Backup/BackupPre_"+j.fit_type+"; mv PreEstimate_"+j.fit_type+" Backup/BackupPre")
        # Use the combinations to modify the cfg, generate directories with the given name, and submit the jobs.
        for c in combinations:
            # pT > 2 GeV/c
            prepare_job(c, j.fit_type, 0., 1./2., j.stub_coordinates, j.track_parameters, True)


job_types = [JobType("Transverse", '"phi"', '"charge/pt"')]
job_types.append(JobType("Transverse_Pz", '"phi"', '"chargeOverPz"'))
job_types.append(JobType("Longitudinal", '"z"', '"cotTheta"'))
job_types.append(JobType("Longitudinal_Rz", '"R", "z"', '"cotTheta"'))

prepare_pre_estimate_jobs(job_types, combinations)



os.system("rm -rf Backup_Old; mv Backup Backup_Old; mkdir Backup; mv Combinations Backup")




# prepare_job(c, "Transverse_lowPt", 1./10., 1./2., '"phi"', '"charge/pt"', True)
# prepare_job(c, "Transverse_highPt", 0., 1./10., '"phi"', '"charge/pt"', True)
# prepare_job(c, "Transverse_lowPt_15", 1./15., 1./2., '"phi"', '"charge/pt"', True)
# prepare_job(c, "Transverse_highPt_15", 0., 1./15., '"phi"', '"charge/pt"', True)
