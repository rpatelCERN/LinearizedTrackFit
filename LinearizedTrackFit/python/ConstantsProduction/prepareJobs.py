import os

__author__ = 'demattia'


# The executable used to run training and testing
executable = "/Users/demattia/RemoteProjects/RemoteProjects"
files_dir = "/Users/demattia/RemoteProjects/"
layers_list = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
radius_cut_min = [0.]*len(layers_list)
radius_cut_max = [1000.]*len(layers_list)
# queue
jobs_queue = []


def array_to_string(a, separator):
    s_a = ""
    for x in a:
        s_a += str(x)+separator
    return s_a.rstrip(separator)


def prepare_cfg(name, job_dir, job_type, input_file_name, train_start, train_end, test_start, test_end,
                stub_coordinates, track_parameters, layers,
                pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                linear_fit_low_pt_file_name, linear_fit_high_pt_file_name, linear_fit_longitudinal_file_name,
                oneOverPt_min, oneOverPt_max, input_regions_number):
    configuration_file = open(job_dir+"conf_"+name+".txt", "w")

    configuration_file.write(str(job_type)+"\n")
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
    configuration_file.write(linear_fit_low_pt_file_name+"\n")
    configuration_file.write(linear_fit_high_pt_file_name+"\n")
    configuration_file.write(linear_fit_longitudinal_file_name+"\n")
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
    configuration_file.write(str(input_regions_number)+"\n")

    configuration_file.close()


def prepare_job(layers, input_file_name, base_dir, radius_cut_min, radius_cut_max,
                fit_type, oneOverPt_min, oneOverPt_max,
                stub_coordinates, track_parameters,
                pre_estimate_transverse, pre_estimate_longitudinal, input_regions_number):
    """
    Function to generate the configuration file for a given combination.
    """
    current_dir = os.getcwd()
    print "input_regions_number = ", input_regions_number
    if input_regions_number == 9:
        current_dir += "/NineRegions"
    elif input_regions_number == 14:
        current_dir += "/FourteenRegions"
    else:
        raise ValueError("Error: regions number can only be 9 or 14. Value requested = "+str(input_regions_number))
    job_dir = current_dir+"/"+base_dir+"/Combinations_"+fit_type+"/"
    pre_estimate_file_transverse = ""
    pre_estimate_file_longitudinal = ""
    if pre_estimate_transverse == "" and pre_estimate_longitudinal == "":
        job_dir = current_dir+"/"+base_dir+"/PreEstimate_"+fit_type+"/"
    else:
        pre_estimate_file_transverse = current_dir+"/"+pre_estimate_transverse
        pre_estimate_file_longitudinal = current_dir+"/"+pre_estimate_longitudinal
    os.system("mkdir -p "+job_dir)
    train_start = 0.
    train_end = 0.5
    test_start = 0.8
    test_end = 1.
    # This is the job type for training
    job_type = 0
    prepare_cfg("train", job_dir, job_type, input_file_name, train_start, train_end, test_start, test_end,
                stub_coordinates, track_parameters, layers,
                pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                "", "", "",
                oneOverPt_min, oneOverPt_max, input_regions_number)
    # This is the job type for testing
    job_type = 1
    prepare_cfg("test", job_dir, job_type, input_file_name, train_start, train_end, test_start, test_end,
                stub_coordinates, track_parameters, layers,
                pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                "", "", "",
                oneOverPt_min, oneOverPt_max, input_regions_number)
    job_file = open(job_dir+"jobFile.sh", "w")
    job_file.write("#!/bin/sh\n")
    job_file.write("cd "+job_dir+"\n")
    job_file.write(executable+" conf_train.txt &> log_train.txt\n")
    job_file.write(executable+" conf_test.txt &> log_test.txt\n")
    job_file.close()
    jobs_queue.append(job_dir+"jobFile.sh")
    os.system("chmod +x "+job_dir+"jobFile.sh")


def prepare_job_ltf(layers, input_file_name, base_dir, radius_cut_min, radius_cut_max,
                    fit_type, oneOverPt_min, oneOverPt_max,
                    stub_coordinates, track_parameters,
                    pre_estimate_transverse, pre_estimate_longitudinal,
                    linear_fit_low_pt, linear_fit_high_pt, linear_fit_longitudinal,
                    input_regions_number):
    """
    Function to generate the configuration file for the linearized track fitter.
    """
    current_dir = os.getcwd()
    if input_regions_number == 9:
        current_dir += "/NineRegions"
    elif input_regions_number == 14:
        current_dir += "/FourteenRegions"
    else:
        raise ValueError("Error: regions number can only be 9 or 14. Value requested = "+str(input_regions_number))
    job_dir = current_dir+"/"+base_dir+"/LTF_"+fit_type+"/"
    pre_estimate_file_transverse = current_dir+"/"+pre_estimate_transverse
    pre_estimate_file_longitudinal = current_dir+"/"+pre_estimate_longitudinal
    linear_fit_low_pt_file_name = current_dir+"/"+linear_fit_low_pt
    linear_fit_high_pt_file_name = current_dir+"/"+linear_fit_high_pt
    linear_fit_longitudinal_file_name = current_dir+"/"+linear_fit_longitudinal
    os.system("mkdir -p "+job_dir)
    train_start = 0.
    train_end = 0.5
    test_start = 0.5
    test_end = 1.
    # This is the type of the linearized track fitter
    job_type = 2
    prepare_cfg("ltf", job_dir, job_type, input_file_name, train_start, train_end, test_start, test_end,
                stub_coordinates, track_parameters, layers,
                pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                linear_fit_low_pt_file_name, linear_fit_high_pt_file_name, linear_fit_longitudinal_file_name,
                oneOverPt_min, oneOverPt_max, input_regions_number)
    job_file = open(job_dir+"jobFile.sh", "w")
    job_file.write("#!/bin/sh\n")
    job_file.write("cd "+job_dir+"\n")
    job_file.write(executable+" conf_ltf.txt &> log_ltf.txt\n")
    job_file.close()
    jobs_queue.append(job_dir+"jobFile.sh")
    os.system("chmod +x "+job_dir+"jobFile.sh")


class JobType:
    def __init__(self, fit_type, stub_coordinates, track_parameters,
                 linear_fit_low_pt_file_name="",
                 linear_fit_high_pt_file_name="",
                 linear_fit_longitudinal_file_name=""):
        self.fit_type = fit_type
        self.stub_coordinates = stub_coordinates
        self.track_parameters = track_parameters
        self.linear_fit_low_pt_file_name = linear_fit_low_pt_file_name
        self.linear_fit_high_pt_file_name = linear_fit_high_pt_file_name
        self.linear_fit_longitudinal_file_name = linear_fit_longitudinal_file_name


def prepare_all_jobs(job_types, pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                     file_name, base_dir, input_regions_number, pt_min=2., pt_max=0., ltf=False):
    oneOverPtMin = 0.
    oneOverPtMax = 1./pt_min
    if pt_max != 0.:
        oneOverPtMin = 1/pt_max
    backup_dir = "Backup/Backup"
    prefix = "Combinations_"
    if not ltf:
        pre_estimate = pre_estimate_file_transverse != "" and pre_estimate_file_longitudinal != ""
        if pre_estimate:
            backup_dir += "Pre"
            prefix = "Pre_estimate_"
    else:
        backup_dir += "LTF"
        prefix = "LTF_"
    for j in job_types:
        os.system("rm -rf "+backup_dir+"_"+j.fit_type+"_Old; mv "+backup_dir+"_"+j.fit_type+" "+backup_dir+"_"+j.fit_type+"_Old")
        job_dir = prefix+j.fit_type+"/"
        os.system("mkdir -p "+backup_dir+"_"+j.fit_type+"; mv "+job_dir+" "+backup_dir)
        os.system("mkdir -p "+job_dir)
        # Use the combinations to modify the cfg, generate directories with the given name, and submit the jobs.
        if not ltf:
            prepare_job(layers_list, file_name, base_dir, radius_cut_min, radius_cut_max,
                        j.fit_type, oneOverPtMin, oneOverPtMax, j.stub_coordinates,
                        j.track_parameters, pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                        input_regions_number)
        else:
            prepare_job_ltf(layers_list, file_name, base_dir, radius_cut_min, radius_cut_max,
                            j.fit_type, oneOverPtMin, oneOverPtMax, j.stub_coordinates,
                            j.track_parameters, pre_estimate_file_transverse, pre_estimate_file_longitudinal,
                            j.linear_fit_low_pt_file_name, j.linear_fit_high_pt_file_name, j.linear_fit_longitudinal_file_name,
                            input_regions_number)

