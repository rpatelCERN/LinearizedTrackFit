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

process.TFileService = cms.Service("TFileService", fileName = cms.string("matrixHists.root") )

process.demo = cms.EDAnalyzer('BuildLinearizedTrackFitMatrix',
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_1p5_5_Eta_0_03_Phi_m004_004_z0_m15p15/extracted.root"),
                              # InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_10_20_Eta_01_012_phi_004_01_z0_m1p1/extracted.root"),
                              # Full trigger tower
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_2_10_Eta_m04p06_Phi_m05p1_3_z0_m15p15/extracted.root"),

                              # Specify to select a single sector (made of single modules)
                              SingleModules = cms.bool(True),

                              # Tracks are required to have at least one stub in this number of layers. More stubs (e.g. overlaps) are neglected, but the track is still used.
                              RequiredLayers = cms.uint32(6),


                              # Select layers to use (must be consistent with the required layers
                              LayersPhi = cms.vint32(5, 6, 7, 8, 9, 10);
                              LayersR = cms.vint32(5, 6, 7, 8, 9, 10);
                              LayersZ = cms.vint32(5, 6, 7, 8, 9, 10);

                              VariableNames = cms.vstring("phi"),
                              TrackParameterNames = cms.vstring("charge/pt", "phi"),
                              # VariableNames = cms.vstring("z"),
                              # TrackParameterNames = cms.vstring("cotTheta", "z0"),

                              # Geometric cuts
                              OneOverPtMin = cms.double(1./1000.),
                              OneOverPtMax = cms.double(1./2.),
                              OneOverPtRegions = cms.int32(1),
                              PhiMin = cms.double(-4.),
                              PhiMax = cms.double(4.),
                              PhiRegions = cms.int32(1),
                              EtaMin = cms.double(-4.),
                              EtaMax = cms.double(4.),
                              EtaRegions = cms.int32(1),
                              Z0Min = cms.double(-15.),
                              Z0Max = cms.double(15),
                              Z0Regions = cms.int32(1),
                              # Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
                              ChargeRegions = cms.int32(2),
                              
                              BuildMatrix = cms.bool(True),
                              # Fraction of the events in the input file to use. Use only 1/2 of them so that the testing phase can use the rest as a statistically independent sample.
                              EventsFractionStartBuild = cms.double(0.),
                              EventsFractionEndBuild = cms.double(0.5),

                              TestMatrix = cms.bool(True),
                              EventsFractionStartTest = cms.double(0.5),
                              EventsFractionEndTest = cms.double(1.)
)


process.p = cms.Path(process.demo)
