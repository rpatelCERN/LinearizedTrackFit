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
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_Pt_1p5_5_Eta_0_03_Phi_m004_004_z0_m15p15/extracted.root"),
                              # Fraction of the events in the input file to use. Use only 1/2 of them so that the testing phase can use the rest as a statistically independent sample.
                              EventsFraction = cms.double(0.5),
                              # Tracks are required to have at least one stub in this number of layers. More stubs (e.g. overlaps) are neglected, but the track is still used.
                              RequiredLayers = cms.uint32(6),
                              VariableNames = cms.vstring("z", "phi"),
                              TrackParameterNames = cms.vstring("1/pt", "phi", "cotTheta")
)


process.p = cms.Path(process.demo)
