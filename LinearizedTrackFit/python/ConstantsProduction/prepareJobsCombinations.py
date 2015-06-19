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
        # self.set_radius_cuts()

    # def set_radius_cuts(self):
    #     first = -1
    #     # Count the number of barrel layers
    #     barrel_layers = 0
    #     for l in self.layers:
    #         if l < 11:
    #             barrel_layers += 1
    #     if barrel_layers >= 4:
    #         return
    #     # Add radius cuts based on the configuration
    #     for i in range(len(self.layers)):
    #         if 10 < self.layers[i] < 16:
    #             if first == -1:
    #                 self.radius_cut_min[self.layers[i]-5] = 61.
    #                 # If we have 3 barrel layers we do not need to add extra cuts
    #                 if barrel_layers >= 3:
    #                     break
    #                 first = i
    #             else:
    #                 self.radius_cut_min[self.layers[first]-5] = 0.
    #                 self.radius_cut_max[self.layers[first]-5] = 61.
    #                 self.radius_cut_min[self.layers[i]-5] = 61.
    #                 break


# Function to generate all combinations from a given set of 6 layers or disks
# No repetition of the same layer is allowed in the input layers_list
def generate_layers_combinations(base_name, layers_list, radius_cut_min, radius_cut_max, input_file_name):
    combinations = [Combination(base_name+"_All", layers_list, radius_cut_min, radius_cut_max, input_file_name)]
    for layer in layers_list:
        name = base_name+"_Removed_" + str(layer)
        layers_combination = [x for x in layers_list if x != layer]
        combinations.append(Combination(name, layers_combination, radius_cut_min, radius_cut_max, input_file_name))

    return combinations


# Function to generate the configuration file for a given combination
def prepare_job(c, pt_min, pt_max, stub_coordinates, track_parameters, pre_estimate):
    print "preparing files for", c.name
    job_dir = "Combinations/"+c.name
    pre_estimate_file_transverse = ""
    pre_estimate_file_longitudinal = ""
    if pre_estimate:
        job_dir = "PreEstimate/"+c.name
    else:
        pre_estimate_file_transverse = "PreEstimate/"+c.name+"/file.txt"
        pre_estimate_file_longitudinal = "PreEstimate/"+c.name+"/file.txt"
    os.system("mkdir -p "+job_dir)
    configuration_file = open("BuildLinearizedTrackFitMatrix_cfg.py", "w")
    for line in open("BuildLinearizedTrackFitMatrix_template_cfg.py"):
        line = line.replace("INPUT_FILE_NAME", c.input_file_name)
        line = line.replace("ONE_OVER_PT_MIN", 1./pt_max)
        line = line.replace("ONE_OVER_PT_MAX", 1./pt_min)
        line = line.replace("STUB_COORDINATES", stub_coordinates)
        line = line.replace("TRACK_PARAMETERS", track_parameters)
        line = line.replace("CHARGE_REGIONS", 1)
        line = line.replace("LAYERS_LIST", c.layers)
        line = line.replace("FIRST_ORDER_CHARGEOVERPT_COEFFICIENTS_FILE_NAME", pre_estimate_file_transverse)
        line = line.replace("FIRST_ORDER_COTTHETA_COEFFICIENTS_FILE_NAME", pre_estimate_file_longitudinal)
        configuration_file.write(line)
    configuration_file.close()
    slurm_file = open(job_dir+"jobFile.slrm", "w")
    for line in open("jobFile_template.slrm"):
        line = line.replace("WORKING_DIR", os.getcwd()+"/"+job_dir)
        slurm_file.write(line)
    slurm_file.close()

# Function to make the combinations unique. If there is a duplicate take the one with the biggest region index since
# this ensures that radius cuts are applied.
# We split this step for now. We instead run all the combinations and use the duplicates to check that their results
# are the same.


# It has to be able to run with a pre- or post-fix to derive the pre-estimate coefficients. When running
# the second PCA it will look for the pre-estimates in the directories with the same name of the current
# ones but with the added pre- or post-fix.
# The layers and disks must be given in order from the smallest index to the biggest index.


files_dir = "/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/"
file_barrel = files_dir+"Muons_Pt_2_more_Eta_m01p09_Phi_m05p1_3_z0_m15p15"
file_hybrid = files_dir+"Muons_Pt_2_more_Eta_0p7_1p7_Phi_m05p1_3_z0_m15p15_HYBRID"
file_endcaps = files_dir+"Muons_Pt_2_more_Eta_1p5_2p5_Phi_m05p1_3_z0_m15p15_ENDCAPS"


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


# cleaning
os.system("rm -rf BackupPre_Transverse_Old; mv BackupPre_Transverse BackupPre_Transverse_Old")
os.system("mkdir BackupPre_Transverse; mv PreEstimate_Transverse Backup")
# Use the combinations to modify the cfg, generate directories with the given name, and submit the jobs.
for c in combinations:
    prepare_job(c, '"phi"', '"charge/pt"')

os.system("rm -rf BackupPre_Longitudinal_Old; mv BackupPre_Longitudinal BackupPre_Longitudinal_Old")
os.system("mkdir BackupPre_Longitudinal; mv PreEstimate_Longitudinal Backup")
# Use the combinations to modify the cfg, generate directories with the given name, and submit the jobs.
for c in combinations:
    print c.name.split("_")[1]
    prepare_job(c, '"phi"', '"charge/pt"')



os.system("rm -rf Backup_Old; mv Backup Backup_Old; mkdir Backup; mv Combinations Backup")
