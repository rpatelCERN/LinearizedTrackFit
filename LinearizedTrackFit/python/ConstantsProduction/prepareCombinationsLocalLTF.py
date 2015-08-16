from prepareJobs import *

__author__ = 'demattia'


ltf_stub_coordinates = ["phi", "R", "z"]
ltf_track_parameters = ["charge/pt", "phi", "cotTheta", "z0"]

file_full_list = {"/FlatOneOverPt": files_dir+"extracted_fullTracker_bigger.root",
                  "/FlatPt": files_dir+"extracted_flatPt.root"}


constants_list = [
    ["FullCorrections",
     "FlatOneOverPt/PreEstimate_Transverse",
     "FlatOneOverPt/PreEstimate_Longitudinal_Rz",
     "FlatOneOverPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10",
     "FlatOneOverPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatOneOverPtPreEstimates",
     "FlatOneOverPt/Combinations_Longitudinal_Rz"],
    ["FullCorrections_flatPtHighPt",
     "FlatOneOverPt/PreEstimate_Transverse",
     "FlatOneOverPt/PreEstimate_Longitudinal_Rz",
     "FlatOneOverPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10",
     "FlatPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatOneOverPtPreEstimates",
     "FlatOneOverPt/Combinations_Longitudinal_Rz"],
    ["FullCorrections_flatPtHighPtAndTgTheta",
     "FlatOneOverPt/PreEstimate_Transverse",
     "FlatPt/PreEstimate_Longitudinal_Rz",
     "FlatOneOverPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10_flatTgTheta",
     "FlatPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more_flatTgTheta",
     "FlatOneOverPt/Combinations_Longitudinal_Rz"],
    ["FullCorrections_flatPtHighPtAndTgThetaAndPtPreEstimate",
     "FlatPt/PreEstimate_Transverse",
     "FlatPt/PreEstimate_Longitudinal_Rz",
     "FlatOneOverPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_2_10_flatOneOverPt",
     "FlatPt/Combinations_Transverse_SecondOrder_ExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_10_more",
     "FlatOneOverPt/Combinations_Longitudinal_Rz"]
]


def prepare_linearized_track_fit():
    """
    Prepare and submit jobs for the LTF.
    """
    regions_number_names = {9: "NineRegions", 14: "FourteenRegions"}

    for regions_number in regions_number_names:
        regions_name = regions_number_names[regions_number]

        for pt_distribution_name in file_full_list:
            file_full = file_full_list[pt_distribution_name]

            for c in constants_list:

                job_types = [JobType(c[0]+"_2_5", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=2., pt_max=5., ltf=True)

                job_types = [JobType(c[0]+"_5_15", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=5., pt_max=15., ltf=True)

                job_types = [JobType(c[0]+"_15_more", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=15., ltf=True)

                job_types = [JobType(c[0]+"_15_25", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=15., pt_max=25., ltf=True)

                job_types = [JobType(c[0]+"_25_50", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=25., pt_max=50., ltf=True)

                job_types = [JobType(c[0]+"_50_100", ltf_stub_coordinates, ltf_track_parameters, c[3], c[4], c[5])]
                prepare_all_jobs(job_types, regions_name+c[1], regions_name+c[2], file_full,
                                 pt_distribution_name, regions_number, pt_min=50., pt_max=100., ltf=True)
