from prepareJobs import *

__author__ = 'demattia'


def prepare_jobs_pre_estimates(input_file_full, base_dir, input_regions_number, fit_range="", pt_min=2., pt_max=0.):
    """
    Prepare and submit jobs for the pre-estimates of pt (or pz), cot(theta) and tan(theta).
    """
    job_types_pre = []
    job_types_pre.append(JobType("Transverse"+fit_range, ["phi"], ["charge/pt"]))
    # job_types_pre.append(JobType("Transverse_Pz", ["phi"], ["chargeOverPz"]))
    job_types_pre.append(JobType("Longitudinal_Rz"+fit_range, ["R", "z"], ["cotTheta", "tgTheta"]))
    prepare_all_jobs(job_types_pre, "", "", input_file_full, base_dir, input_regions_number,
                     pt_min=pt_min, pt_max=pt_max)


def prepare_pre_estimates():
    regions_numbers = [9, 14]

    for regions_number in regions_numbers:
        # Flat pT
        # -------
        file_full = files_dir+"extracted_flatPt.root"
        prepare_jobs_pre_estimates(file_full, "/FlatPt", regions_number)
        prepare_jobs_pre_estimates(file_full, "/FlatPt", regions_number, fit_range="_2_10", pt_min=2., pt_max=10.)
        prepare_jobs_pre_estimates(file_full, "/FlatPt", regions_number, fit_range="_10_more", pt_min=10.)

        # Flat 1/pT
        # ---------
        file_full = files_dir+"extracted_fullTracker_bigger.root"
        prepare_jobs_pre_estimates(file_full, "/FlatOneOverPt", regions_number)
        prepare_jobs_pre_estimates(file_full, "/FlatOneOverPt", regions_number, fit_range="_2_10", pt_min=2., pt_max=10.)
        prepare_jobs_pre_estimates(file_full, "/FlatOneOverPt", regions_number, fit_range="_10_more", pt_min=10.)
