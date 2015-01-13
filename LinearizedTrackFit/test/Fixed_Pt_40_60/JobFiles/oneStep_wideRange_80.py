#########################
#
# Configuration file for simple PGUN events
# production in tracker only geometry
#
# Instruction to run this script are provided on this page:
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP II
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 12/04/2013
# Maj. modif  : 17/06/2013 (adding the official stub producer)
# Maj. modif  : 10/01/2014 (going to new CMSSW release)
#
# Script tested with release CMSSW_6_2_0_SLHC20
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')

# Special geometry (Tracker only)
process.load('DataProduction.SkimGeometry.Sim_SKIM_cff')
process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')
process.load('DataProduction.SkimGeometry.mixNoPU_SKIM_cfi')
process.load('DataProduction.SkimGeometry.Digi_SKIM_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'DES23_62_V1::All', '')

# Random seeds
process.RandomNumberGeneratorService.generator.initialSeed      = 80
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = 81
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = 82
process.RandomNumberGeneratorService.mix.initialSeed            = 83

# Generate particle gun events

# Generate particle gun events
process.generator = cms.EDProducer("FlatRandomOneOverPtGunProducer",
    PGunParameters = cms.PSet(
	# 40 < pt < 60 GeV
	MinOneOverPt = cms.double(1./60.),
        MaxOneOverPt = cms.double(1./40.),
#	XFlatSpread = cms.double(1.5),  # In mm (requires an update 
#	YFlatSpread = cms.double(1.5),  # In mm  of the official 
#	ZFlatSpread = cms.double(150.), # In mm  PGUN code, see tutorial)
	towerID= cms.int32(-1),           # Trigger tower ID (put -1 for default params)
	XFlatSpread = cms.double(0.),  # In mm (requires an update 
	YFlatSpread = cms.double(0.),  # In mm  of the official 
	ZFlatSpread = cms.double(0.), # In mm  PGUN code, see tutorial)
        PartID = cms.vint32(-13),
	# Eta = 0.2
        MinEta = cms.double(0.2),
        MaxEta = cms.double(0.2),
	# Phi = 0
        MinPhi = cms.double(0.),
	MaxPhi = cms.double(0.)
    ),
    Verbosity = cms.untracked.int32(0),
    AddAntiParticle = cms.bool(False),
)

# Load the extractor
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

# Tune some options (see MIB_extractor_cfi.py for details)

process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doPixel          = True
process.MIBextraction.doL1TT           = True
process.MIBextraction.doMatch          = True
process.MIBextraction.extractedRootFile = "extracted_80.root"

process.MIBextraction.analysisSettings = cms.untracked.vstring(
    "matchedStubs 0",
    "posMatching  1",
    "zMatch  0",
    "maxClusWdth  4",
    "windowSize   10",
    "pdgSel -1",
    "verbose 0"
    )


process.generation_step      = cms.Path(process.pgen)
process.simulation_step      = cms.Path(process.psim)
process.genfiltersummary_step= cms.EndPath(process.genFilterSummary)
process.digitisation_step    = cms.Path(process.pdigi)
process.L1TrackTrigger_step  = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step  = cms.Path(process.TrackTriggerAssociatorClustersStubs)
process.MIBextraction_step = cms.Path(process.MIBextraction)
process.endjob_step          = cms.EndPath(process.endOfProcess)


process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1TrackTrigger_step,process.L1TTAssociator_step,process.endjob_step,process.MIBextraction_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq
	
# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
