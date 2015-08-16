from prepareJobs import *

__author__ = 'demattia'


def prepare_second_estimates():
    """
    Prepare and submit jobs for the second PCA.
    """
    regions_number_names = {9: "NineRegions", 14: "FourteenRegions"}

    for regions_number in regions_number_names:
        regions_name = regions_number_names[regions_number]

        # Flat 1/pT
        # ---------
        file_full = files_dir+"extracted_fullTracker_bigger.root"

        # job_types = [JobType("Longitudinal_Rz", ["CorrectedZSecondOrder"], ["cotTheta", "z0"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatOneOverPt", regions_number)

        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatOneOverPt", regions_number, pt_min=10.)
        #
        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_alsoPreEstimates",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse_10_more",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz_10_more", file_full,
        #                  "/FlatOneOverPt", regions_number, pt_min=10.)
        #
        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatPtPreEstimates",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatOneOverPt", regions_number, pt_min=10.)
        #
        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_15_more",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatOneOverPt", regions_number, pt_min=15.)

        job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatTgTheta",
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
        prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
                         regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
                         "/FlatOneOverPt", regions_number, pt_min=10.)

        job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10",
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
        prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
                         regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
                         "/FlatOneOverPt", regions_number, pt_min=2., pt_max=10.)

        job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10_flatTgTheta",
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
        prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
                         regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
                         "/FlatOneOverPt", regions_number, pt_min=2., pt_max=10.)

        # Flat pT
        # -------
        file_full = files_dir+"extracted_flatPt.root"

        # job_types = [JobType("Longitudinal_Rz", ["CorrectedZSecondOrder"], ["cotTheta", "z0"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatPt", regions_number)

        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatPt", regions_number, pt_min=10.)

        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatOneOverPtPreEstimates",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatPt", regions_number, pt_min=10.)
        # job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_15_more",
        #                      ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
        #                      ["charge/pt", "phi"])]
        # prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
        #                  regions_name+"/FlatOneOverPt/PreEstimate_Longitudinal_Rz", file_full,
        #                  "/FlatPt", regions_number, pt_min=15.)

        job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatOneOverPt_PtPreEstimate",
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
        prepare_all_jobs(job_types, regions_name+"/FlatOneOverPt/PreEstimate_Transverse",
                         regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
                         "/FlatPt", regions_number, pt_min=2., pt_max=10.)

        job_types = [JobType("Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10_flatOneOverPt",
                             ["CorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup"],
                             ["charge/pt", "phi"])]
        prepare_all_jobs(job_types, regions_name+"/FlatPt/PreEstimate_Transverse",
                         regions_name+"/FlatPt/PreEstimate_Longitudinal_Rz", file_full,
                         "/FlatPt", regions_number, pt_min=2., pt_max=10.)
