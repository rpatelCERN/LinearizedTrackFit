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

process.demo = cms.EDAnalyzer('LinearizedTrackFit',
                              # InputFileName = cms.string("/afs/cern.ch/user/d/demattia/work/Upgrade/CMSSW_6_2_0_SLHC20_patch1/src/Fixed/extracted.root")
                              # InputFileName = cms.string("/afs/cern.ch/user/d/demattia/work/Upgrade/CMSSW_6_2_0_SLHC20_patch1/src/Fixed/extracted_highPt.root")
                              # InputFileName = cms.string("/afs/cern.ch/user/d/demattia/work/Upgrade/CMSSW_6_2_0_SLHC20_patch1/src/Wide/extracted.root")
                              # InputFileName = cms.string("/afs/cern.ch/user/d/demattia/work/Upgrade/CMSSW_6_2_0_SLHC20_patch1/src/EtaZWide/extracted.root")
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_wide/extracted.root")
)


process.p = cms.Path(process.demo)
