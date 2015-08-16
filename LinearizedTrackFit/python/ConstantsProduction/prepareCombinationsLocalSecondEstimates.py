from prepareJobs import *

__author__ = 'demattia'


file_full_list = {"/FlatOneOverPt": files_dir+"extracted_fullTracker_bigger.root",
                  "/FlatPt": files_dir+"extracted_flatPt.root"}


constants_list_second_estimates = [
    ["FullCorrections",
     "",
     "FlatOneOverPt/PreEstimate_Transverse",
     "FlatOneOverPt/PreEstimate_Longitudinal_Rz"],
    ["FullCorrections",
     "_flatPtTgTheta",
     "FlatOneOverPt/PreEstimate_Transverse",
     "FlatPt/PreEstimate_Longitudinal_Rz"],
    ["FullCorrections",
     "_flatPtTgThetaAndPtPreEstimate",
     "FlatPt/PreEstimate_Transverse",
     "FlatPt/PreEstimate_Longitudinal_Rz"]
]



def prepare_second_estimates():
    """
    Prepare and submit jobs for the second PCA.
    """
    regions_number_names = {9: "NineRegions", 14: "FourteenRegions"}

    for regions_number in regions_number_names:
        regions_name = regions_number_names[regions_number]

        for pt_distribution_name in file_full_list:
            file_full = file_full_list[pt_distribution_name]

            # job_types = [JobType("Longitudinal_Rz", ["CorrectedZSecondOrder"], ["cotTheta", "z0"])]
            # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
            #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
            #                  "/FlatOneOverPt", regions_number)

            for c in constants_list_second_estimates:

                job_types = [JobType(c[0]+"_2_10"+c[1],
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
                prepare_all_jobs(job_types, regions_name+"/"+c[2],
                                 regions_name+"/"+c[3], file_full,
                                 regions_name, regions_number, pt_min=2., pt_max=10.)

                job_types = [JobType(c[0]+"_10_more"+c[1],
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
                prepare_all_jobs(job_types, regions_name+"/"+c[2],
                                 regions_name+"/"+c[3], file_full,
                                 regions_name, regions_number, pt_min=10.)
