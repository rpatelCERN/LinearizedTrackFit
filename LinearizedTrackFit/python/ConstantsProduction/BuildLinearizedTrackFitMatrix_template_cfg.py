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

                              # Templated input file
                              InputFileName = cms.string("INPUT_FILE_NAME"),

                              # Specify to select a single sector (made of single modules)
                              SingleModules = cms.bool(False),

                              # Select layers to use for each stub coordinate
                              LayersAll = cms.vint32(LAYERS_LIST),

                              # Radius cuts
                              RadiusCutsMin = cms.vdouble(RADIUS_CUTS_MIN),
                              RadiusCutsMax = cms.vdouble(RADIUS_CUTS_MAX),

                              VariableNames = cms.vstring(STUB_COORDINATES),
                              TrackParameterNames = cms.vstring(TRACK_PARAMETERS),

                              # Geometric cuts
                              OneOverPtMin = cms.double(ONE_OVER_PT_MIN),
                              OneOverPtMax = cms.double(ONE_OVER_PT_MAX),
                              OneOverPtRegions = cms.int32(1),
                              PhiMin = cms.double(0.),
                              PhiMax = cms.double(0.8),
                              PhiRegions = cms.int32(1),
                              EtaMin = cms.double(-3.),
                              EtaMax = cms.double(3.),
                              EtaRegions = cms.int32(1),
                              Z0Min = cms.double(-15.),
                              Z0Max = cms.double(15),
                              Z0Regions = cms.int32(1),
                              # Specify 1 for no charge splitting and 2 for separating positive and negative charge in difference regions
                              ChargeRegions = cms.int32(CHARGE_REGIONS),
                              # Specify 0 for no endcap region splitting and != 0 for endcap regions splitting
                              EndcapRegions = cms.int32(0),
                              
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
                              EventsFractionEndTest = cms.double(1.),

                              FirstOrderChargeOverPtCoefficientsFileName = cms.string("FIRST_ORDER_CHARGEOVERPT_COEFFICIENTS_FILE_NAME"),
                              FirstOrderCotThetaCoefficientsFileName = cms.string("FIRST_ORDER_COTTHETA_COEFFICIENTS_FILE_NAME")
)


process.p = cms.Path(process.demo)
