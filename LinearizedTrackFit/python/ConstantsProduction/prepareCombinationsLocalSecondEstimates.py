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

    # # Old baseline
    # job_types = [JobType("OldBaseline_2_10", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"])]
    # prepare_all_jobs(job_types, "FlatOneOverPt/PreEstimate_Transverse", "FlatOneOverPt/PreEstimate_Longitudinal_Rz",
    #                  files_dir+"extracted_fullTracker_bigger.root",
    #                  "/FlatOneOverPt", 9, pt_min=2., pt_max=10.)
    #
    # job_types = [JobType("OldBaseline_10_more", ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"])]
    # prepare_all_jobs(job_types, "FlatOneOverPt/PreEstimate_Transverse", "FlatOneOverPt/PreEstimate_Longitudinal_Rz",
    #                  files_dir+"extracted_fullTracker_bigger.root",
    #                  "/FlatOneOverPt", 9, pt_min=10.)
    #

    # # Old baseline everything flat pT
    # job_types = [JobType("OldBaseline_2_10_flatPtTgThetaAndPtPreEstimate",
    #                      ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"])]
    # prepare_all_jobs(job_types, "FlatPt/PreEstimate_Transverse", "FlatPt/PreEstimate_Longitudinal_Rz",
    #                  files_dir+"extracted_fullTracker_bigger.root",
    #                  "/FlatPt", 9, pt_min=2., pt_max=10.)
    #
    # job_types = [JobType("OldBaseline_10_more_flatPtTgThetaAndPtPreEstimate",
    #                      ["CorrectedPhiSecondOrder"], ["charge/pt", "phi"])]
    # prepare_all_jobs(job_types, "FlatPt/PreEstimate_Transverse", "FlatPt/PreEstimate_Longitudinal_Rz",
    #                  files_dir+"extracted_fullTracker_bigger.root",
    #                  "/FlatPt", 9, pt_min=10.)

    regions_number_names = {9: "NineRegions", 14: "FourteenRegions"}

    for regions_number in regions_number_names:
        # regions_name = regions_number_names[regions_number]

        for pt_distribution_name in file_full_list:
            file_full = file_full_list[pt_distribution_name]

            job_types = [JobType("FullCorrections_GEN_2_10",
                                 ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN"],
                                 ["charge/pt", "phi"])]
            prepare_all_jobs(job_types, "/FlatOneOverPt/PreEstimate_Transverse",
                             "/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
                             pt_distribution_name, regions_number, pt_min=2., pt_max=10.)

            job_types = [JobType("FullCorrections_GEN_10_more",
                                 ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN"],
                                 ["charge/pt", "phi"])]
            prepare_all_jobs(job_types, "/FlatOneOverPt/PreEstimate_Transverse",
                             "/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
                             pt_distribution_name, regions_number, pt_min=10.)
    #
    #         for c in constants_list_second_estimates:
    #
    #             job_types = [JobType("Longitudinal_Rz"+c[1],
    #                                  ["CorrectedZSecondOrder"],
    #                                  ["cotTheta", "z0"])]
    #             prepare_all_jobs(job_types, c[2],
    #                              c[3], file_full,
    #                              pt_distribution_name+"/", regions_number)
    #
    #             job_types = [JobType(c[0]+"_2_10"+c[1],
    #                          ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
    #                          ["charge/pt", "phi"])]
    #             prepare_all_jobs(job_types, c[2],
    #                              c[3], file_full,
    #                              pt_distribution_name+"/", regions_number, pt_min=2., pt_max=10.)
    #
    #             job_types = [JobType(c[0]+"_10_more"+c[1],
    #                          ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
    #                          ["charge/pt", "phi"])]
    #             prepare_all_jobs(job_types, c[2],
    #                              c[3], file_full,
    #                              pt_distribution_name+"/", regions_number, pt_min=10.)
