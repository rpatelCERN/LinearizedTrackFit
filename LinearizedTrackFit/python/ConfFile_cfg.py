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
                              HInvPtErrRangeMin = cms.double(-0.01),
                              HInvPtErrRangeMax = cms.double(0.01),
                              HZ0ErrRangeMin = cms.double(-1.),
                              HZ0ErrRangeMax = cms.double(1.),

                              # PhiMinCut = cms.double(0.),
                              # PhiMaxCut = cms.double(0.04),
                              # InvPtMinCut = cms.double(1./105.),
                              # InvPtMaxCut = cms.double(1./95.),
                              # Z0MinCut = cms.double(-15.),
                              # Z0MaxCut = cms.double(15.),
                              # EtaMinCut = cms.double(0.),
                              # EtaMaxCut = cms.double(0.3),
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/AllOpenSmallRanges/JobFiles/extracted.root")


                              # PhiMinCut = cms.double(-0.1),
                              # PhiMaxCut = cms.double(0.1),
                              # InvPtMinCut = cms.double(1./100.),
                              # InvPtMaxCut = cms.double(1./50.),
                              # Z0MinCut = cms.double(-0.0001),
                              # Z0MaxCut = cms.double(0.0001),
                              # EtaMinCut = cms.double(0.15),
                              # EtaMaxCut = cms.double(0.25),
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/InvPtRange/JobFiles/extracted.root")


                              # PhiMinCut = cms.double(-0.04),
                              # PhiMaxCut = cms.double(0.04),
                              # InvPtMinCut = cms.double(1./4.),
                              # InvPtMaxCut = cms.double(1./1.5),
                              # Z0MinCut = cms.double(-15.),
                              # Z0MaxCut = cms.double(15.),
                              # EtaMinCut = cms.double(0.25),
                              # EtaMaxCut = cms.double(0.35),
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/AllOpenVerySmallPt/JobFiles/extracted.root")


                              # PhiMinCut = cms.double(-0.00001),
                              # PhiMaxCut = cms.double(0.00001),
                              # InvPtMinCut = cms.double(1./100.-0.00001),
                              # InvPtMaxCut = cms.double(1./100.+0.00001),
                              # Z0MinCut = cms.double(-0.00001),
                              # Z0MaxCut = cms.double(0.00001),
                              # EtaMinCut = cms.double(0.2-0.00001),
                              # EtaMaxCut = cms.double(0.2+0.00001),
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/Fixed_Pt_100/JobFiles/extracted.root")


                              PhiMinCut = cms.double(-0.1),
                              PhiMaxCut = cms.double(0.1),
                              InvPtMinCut = cms.double(1./500.-0.00001),
                              InvPtMaxCut = cms.double(1./500.+0.00001),
                              Z0MinCut = cms.double(-15.),
                              Z0MaxCut = cms.double(15.),
                              EtaMinCut = cms.double(0.15),
                              EtaMaxCut = cms.double(0.25),
                              InputFileName = cms.string("/home/demattia/fdata_demattia/Upgrade/LinearizedTrackFit/Stubs/Muons_wide/extracted.root")


                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/InvPtRange/extracted.root")
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/AllOpenSmallerPt/JobFiles/extracted.root")
                              # InputFileName = cms.string("/home/demattia/CMSSW_6_2_0_SLHC20_patch1/src/LinearizedTrackFit/LinearizedTrackFit/test/Fixed_Pt_90_110/JobFiles/extracted.root")
)


process.p = cms.Path(process.demo)
