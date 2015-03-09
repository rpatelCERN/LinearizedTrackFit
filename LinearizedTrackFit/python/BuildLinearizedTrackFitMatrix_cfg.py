import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.source = cms.Source("PoolSource",
#     # replace 'myfile.root' with the source file you want to use
#     fileNames = cms.untracked.vstring(
#         'file:myfile.root'
#     )
# )

# process.TFileService = cms.Service("TFileService", fileName = cms.string("matrixHists.root") )

process.demo = cms.EDAnalyzer('BuildLinearizedTrackFitMatrix',
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_1p5_5_Eta_0_03_Phi_m004_004_z0_m15p15/extracted.root"),
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_10_20_Eta_01_012_phi_004_01_z0_m1p1/extracted.root"),

                              # Bigger range to cover a set of single modules
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_01_03_phi_0_02_z0_m15p15/extracted.root"),

                              # Full trigger tower 2 < pt < 10 GeV/c
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_m04p06_Phi_m05p1_3_z0_m15p15/extracted.root"),

                              # Full trigger tower and full pt range
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_2_more_Eta_m04p06_Phi_m05p1_3_z0_m15p15/JobFiles/extracted_doubleDigits.root"),

                              # Full trigger tower 10 < pt < 30 GeV/c
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_10_30_Eta_m04p06_Phi_m05p1_3_z0_m15p15/JobFiles/extracted.root"),

                              # Smaller range to cover only a set of single modules
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_0_04_phi_m025_025_z0_m15p15/extracted.root"),

                              # OfficialStub
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_2_10_Eta_0_04_phi_m025_025_z0_m1p1/extracted.root"),

                              # Templated input file
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_2_more_Eta_m04p06_Phi_m05p1_3_z0_m15p15/JobFiles/extracted_doubleDigits.root"),


                              # Specify to select a single sector (made of single modules)
                              SingleModules = cms.bool(False),

                              # Select layers to use for each stub coordinate
                              LayersPhi = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersPhiOverR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersChargeSignedPhi = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersGenChargeSignedPhi = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersOneOverR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersChargeSignedR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersChargeCorrectedR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersChargeOverPtCorrectedR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersChargeOverPtCorrectedRCube = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersZ = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersDeltaS = cms.vint32(5, 6, 7, 8, 9, 10),

                              # Fix all phi means to 0 if True
                              FixMeansPhi = cms.bool(False),

                              # To remove outlier stubs. Each number is, for each layer, the maximum allowed distance in cm in
                              # the transverse plane between the stub and the extrapolated generator-level track associated to it.
                              DistanceCutsTransverse = cms.vdouble(1., 2., 4., 5., 8., 10.),
                              DistanceCutsLongitudinal = cms.vdouble(0.3, 0.5, 0.8, 3., 3.5, 4.),

                              VariableNames = cms.vstring("phi", "R"),
                              TrackParameterNames = cms.vstring("charge/pt", "phi"),
                              # VariableNames = cms.vstring("z", "R"),
                              # TrackParameterNames = cms.vstring("cotTheta", "z0"),
                              # LayersZ = cms.vint32(5, 6, 7),

                              # Geometric cuts
                              OneOverPtMin = cms.double(0.2),
                              OneOverPtMax = cms.double(0.5),
                              OneOverPtRegions = cms.int32(1),
                              PhiMin = cms.double(0.),
                              PhiMax = cms.double(0.8),
                              PhiRegions = cms.int32(1),
                              EtaMin = cms.double(-0.3),
                              EtaMax = cms.double(0.5),
                              EtaRegions = cms.int32(1),
                              Z0Min = cms.double(-15.),
                              Z0Max = cms.double(15),
                              Z0Regions = cms.int32(1),
                              # Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
                              ChargeRegions = cms.int32(1),
                              
                              BuildMatrix = cms.bool(True),
                              # Fraction of the events in the input file to use. Use only 1/2 of them so that the testing phase can use the rest as a statistically independent sample.
                              EventsFractionStartBuild = cms.double(0.),
                              EventsFractionEndBuild = cms.double(0.5),
                              # Optional for the buildMatrix step
                              MapSectors = cms.bool(False),
                              ComputeDistances = cms.bool(False),
                              ComputeCorrelations = cms.bool(True),
                              PhiSymmetricFit = cms.bool(False),
                              UsePcs = cms.bool(True),

                              TestMatrix = cms.bool(True),
                              EventsFractionStartTest = cms.double(0.5),
                              EventsFractionEndTest = cms.double(1.)
)


process.p = cms.Path(process.demo)
