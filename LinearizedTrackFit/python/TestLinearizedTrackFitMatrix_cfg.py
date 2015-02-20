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

process.demo = cms.EDAnalyzer('TestLinearizedTrackFitMatrix',
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_1p5_5_Eta_0_03_Phi_m004_004_z0_m15p15/extracted.root"),
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_10_20_Eta_01_012_phi_004_01_z0_m1p1/extracted.root"),

                              # Bigger range to cover a set of single modules
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_01_03_phi_0_02_z0_m15p15/extracted.root"),

                              # Full trigger tower 2 < pt < 10 GeV/c
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_m04p06_Phi_m05p1_3_z0_m15p15/extracted.root"),

                              # Full trigger tower and full pt range
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_2_more_Eta_m04p06_Phi_m05p1_3_z0_m15p15/JobFiles/extracted_doubleDigits.root"),

                              # Full trigger tower 10 < pt < 30 GeV/c
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_10_30_Eta_m04p06_Phi_m05p1_3_z0_m15p15/extracted.root"),

                              # Smaller range to cover only a set of single modules
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_0_04_phi_m025_025_z0_m15p15/extracted.root"),

                              # OfficialStub
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Official/Muons_Pt_2_10_Eta_0_04_phi_m025_025_z0_m1p1/extracted.root"),

                              # Specify to select a single sector (made of single modules)
                              SingleModules = cms.bool(False),
                              UsePcs = cms.bool(True),

                              # Select layers to use for each stub coordinate
                              LayersPhi = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersR = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersZ = cms.vint32(5, 6, 7, 8, 9, 10),
                              LayersDeltaS = cms.vint32(5, 6, 7, 8, 9, 10),

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
                              OneOverPtMin = cms.double(1./5.),
                              OneOverPtMax = cms.double(1./2.),
                              PhiMin = cms.double(0.),
                              PhiMax = cms.double(0.8),
                              EtaMin = cms.double(-0.3),
                              EtaMax = cms.double(0.5),
                              Z0Min = cms.double(-15.),
                              Z0Max = cms.double(15),
                              
                              EventsFractionStart = cms.double(0.5),
                              EventsFractionEnd = cms.double(0.5001)
)


process.p = cms.Path(process.demo)
