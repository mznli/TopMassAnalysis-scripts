"""
Example analysis module: make a dimuon mass plot from a NanoAOD
"""
import logging
import re
import os.path

import bamboo.plots
import numpy as np
from functools import partial
from bamboo.analysismodules import (NanoAODHistoModule, NanoAODModule)
from bamboo import treefunctions as op
from bamboo.scalefactors import get_correction

logger = logging.getLogger(__name__)

def muonDef(era):
    
    def muonDefImpl(mu):
        if era == '2016ULpreVFP' or era == '2016ULpostVFP' or era == '2018UL':
            ptCut = 30.
        elif era == '2017UL':
            ptCut = 29.
        return op.AND(
            mu.pt > ptCut,
            op.abs(mu.eta) < 2.4,
            # tight ID
            mu.tightId,
            # tight deltabeta ISO R=0.4
            mu.pfRelIso04_all < 0.15,
            op.abs(mu.dxy) <= 0.05,
            op.abs(mu.dz) <= 0.1,
            mu.sip3d <= 8
        )
    return muonDefImpl


def vetoMuonDef(mu):
    return op.AND(
        mu.pt > 15., op.abs(mu.eta) < 2.4,
        # loose ID
        mu.looseId,
        # loose deltabeta ISO R=0.4
        mu.pfRelIso04_all < 0.25,
    )

def eleDef(era):
    def eleDefImpl(ele):
        absEtaSC = op.abs(ele.eta + ele.deltaEtaSC)
        ele_pt = ele.pt # / ele.eCorr # uncomment to use uncalibrated electron pt
        if era == '2016ULpreVFP' or era == '2016ULpostVFP':
            eraCut = op.AND(ele_pt > 29, op.abs(ele.eta) < 2.5)
        if era == '2017UL' or era == '2018UL':
            eraCut = op.OR(
                op.AND(ele_pt > 34., op.abs(ele.eta) < 2.5),
                op.AND(ele_pt > 30., op.abs(ele.eta) < 2.1),
            )
        return op.AND(
            eraCut,
            op.OR(absEtaSC < 1.4442, absEtaSC > 1.566),
            # "average" d0 and dz cuts, to be tuned?
            op.OR(
                # barrel
                op.AND(absEtaSC <= 1.479, op.abs(ele.dxy) < 0.05, op.abs(ele.dz) < 0.1),
                # endcap
                op.AND(absEtaSC > 1.479, op.abs(ele.dxy) < 0.1, op.abs(ele.dz) < 0.2),
            ),
            # tight cut-based ID
            ele.cutBased == 4,
        )
    return eleDefImpl

def vetoEleDef(ele):
    ele_pt = ele.pt # / ele.eCorr # uncomment to use uncalibrated electron pt
    return op.AND(
        ele_pt > 15., op.abs(ele.eta) < 2.5,
        # veto cut-based ID
        ele.cutBased >= 1,
    )

def jetDef(jet):
    return op.AND(
        jet.pt > 30., op.abs(jet.eta) < 2.4,
        # tight lepton veto jet ID
        jet.jetId & 4,
        # loose puID for jets with pt < 50
        op.OR(jet.pt > 50, jet.puId > 0)
    )

# Clean jets from leptons -> since we have veto'd extra loose leptons we
# don't have to use vetoElectrons/Muons at this point
# However we NEED rng_any since we don't know how many electrons/muons we have (could be 0 or 1)
# Also make sure the jets are sorted by Pt (not guaranteed since JER is applied)
def cleanJets(jets, muons, electrons, sort=True):
    jets = op.select(jets, lambda jet: op.AND(
            op.NOT(op.rng_any(electrons, lambda ele: op.deltaR(jet.p4, ele.p4) < 0.4)),
            op.NOT(op.rng_any(muons, lambda mu: op.deltaR(jet.p4, mu.p4) < 0.4))
        ))
    if sort:
        jets = op.sort(jets, lambda j: -j.pt)
    return jets


# scale factors

def pogEraFormat(era):
    return era.replace("UL", "") + "_UL"

def localizePOGSF(era, POG, fileName):
    return os.path.join("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration", "POG", POG, pogEraFormat(era), fileName)

def getYearFromEra(era):
    """ Go from '2017'/'2018UL'/'2016ULpreVFP' to '17'/'18'/'16' """
    return re.search(r"20([0-9]+).*", era).group(1)

def makePUWeight(tree, era, selection):
    from bamboo.analysisutils import makePileupWeight
    goldenJSON = f"Collisions{getYearFromEra(era)}_UltraLegacy_goldenJSON"
    puTuple = (localizePOGSF(era, "LUM", "puWeights.json.gz"), goldenJSON)
    return makePileupWeight(puTuple, tree.Pileup_nTrueInt, systName="pileup", sel=selection)

### definitions from: https://gitlab.cern.ch/swertz/ttbbRun2Bamboo/-/blob/master/python/definitions.py?ref_type=heads#L410
# maps name of systematic to name of correction inside of jsons
leptonSFLib = {
    "electron_ID": "UL-Electron-ID-SF",
    "electron_reco": "UL-Electron-ID-SF",
    "electron_trigger": "EleTriggerSF",
    "muon_reco": "NUM_TrackerMuons_DEN_genTracks",
    "muon_ID": "NUM_TightID_DEN_TrackerMuons",
    "muon_iso": "NUM_TightRelIso_DEN_TightIDandIPCut",
    "muon_trigger": {
        "2016ULpreVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        "2016ULpostVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        "2017UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
        "2018UL": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
    },
}

def getLeptonSF(era, systName):
    if systName == "muon_trigger":
        corrName = leptonSFLib[systName][era]
    else:
        corrName = leptonSFLib[systName]

    if "muon" in systName:
        path = localizePOGSF(era, "MUO", "muon_Z.json.gz")
        # uncomment the following to get ele trig sf:
#    elif systName == "electron_trigger":
#        path = os.path.join(
#            os.path.dirname(os.path.abspath(__file__)),
#            "..", "..", "scale-factors", "eleTrigSFs", era + "_EleTriggerSF_NanoAODv9_v0.json")
    elif "electron" in systName:
        path = localizePOGSF(era, "EGM", "electron.json.gz")

    return path, corrName

def getScaleFactor(era, noSel, systName, defineOnFirstUse=True):
    fileName, correction = getLeptonSF(era, systName)

    from bamboo.scalefactors import get_correction
    

    if "muon" in systName:
        etaParam = "eta"
        etaExpr = lambda mu: op.abs(mu.eta)
    elif "electron" in systName:
        etaParam = "eta"
        etaExpr = lambda el: el.eta + el.deltaEtaSC
    else:
        raise ValueError("Only muon or electron SFs are handled here!")

#    if systName == "muon_ID" or systName == "muon_iso" or systName == "muon_trigger" :
#        return get_correction(fileName, correction, params={"pt": lambda mu: mu.pt, etaParam: etaExpr, "year": pogEraFormat(era)},
#                              systParam="ValType", systNomName="sf",
#                              systVariations={f"{systName}up": "systup", f"{systName}down": "systdown"},
#                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    if systName == "muon_reco":  #"muon" in systName:
        return get_correction(fileName, correction, params={"pt": lambda mu: op.max(mu.pt, 40.), etaParam: etaExpr}, ## check if pt threshold was removed
                              systParam="scale_factors", systNomName="nominal",
                              systVariations={f"{systName}up": "systup", f"{systName}down": "systdown"},
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    elif "muon" in systName:
        return get_correction(fileName, correction, params={"pt": lambda mu: mu.pt, etaParam: etaExpr},
                              systParam="scale_factors", systNomName="nominal",
                              systVariations={f"{systName}up": "systup", f"{systName}down": "systdown"},
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    elif systName == "electron_trigger":
        return get_correction(fileName, correction, params={"pt": lambda el: el.pt, etaParam: etaExpr},
                              systParam="sf", systNomName="central",
                              systVariations=("up", "down"), systName=systName,
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)
    else:
        wp = "Tight" if systName == "electron_ID" else "RecoAbove20"
        return get_correction(fileName, correction, params={"pt": lambda el: el.pt, etaParam: etaExpr, "year": era.replace("UL", ""), "WorkingPoint": wp},
                              systParam="ValType", systNomName="sf",
                              systVariations={f"{systName}up": "sfup", f"{systName}down": "sfdown"},
                              defineOnFirstUse=defineOnFirstUse, sel=noSel)

def getL1PrefiringSystematic(tree):
    
    return op.systematic(tree.L1PreFiringWeight_Nom, name="L1prefire", up=tree.L1PreFiringWeight_Up, down=tree.L1PreFiringWeight_Dn)

def muonTriggerDef(HLT, sample, era, isMC):
    
    cuts = []
    if era == '2016ULpreVFP' or era == '2016ULpostVFP':
        cuts.append(op.OR(HLT.IsoMu24, HLT.IsoTkMu24))
    if era == '2017UL':
        cuts.append(HLT.IsoMu27)
    if era == '2018UL':
        cuts.append(HLT.IsoMu24)
    if not isMC:
        cuts.append(op.c_bool("SingleMuon" in sample))
    return cuts

def eleTriggerDef(TrigObj, HLT, ele, sample, era, isMC):
    cuts = []
    if era == '2016ULpreVFP' or era == '2016ULpostVFP':
        cuts.append(HLT.Ele27_WPTight_Gsf)
    if era == '2017UL':
        cuts.append(op.OR(
            op.AND(HLT.Ele32_WPTight_Gsf_L1DoubleEG, op.rng_any(TrigObj, lambda obj: op.AND(op.deltaR(obj.p4, ele.p4) < 0.1, obj.filterBits & 1024))),
            HLT.Ele28_eta2p1_WPTight_Gsf_HT150))
    if era == '2018UL':
        cuts.append(op.OR(HLT.Ele32_WPTight_Gsf, HLT.Ele28_eta2p1_WPTight_Gsf_HT150))
    if not isMC:
        if era == '2016ULpreVFP' or era == '2016ULpostVFP' or era == '2017UL':
            cuts.append(op.c_bool("SingleElectron" in sample))
        if era == '2018UL':
            cuts.append(op.c_bool("EGamma" in sample)) # only called EGamma in 2018
    return cuts


# Return exclusive muon selection (without and with applied trigger), with trigger, ID, and iso scale factors (for MC)
def buildMuonSelections(tree, noSel, muons, vetoMuons, electrons, vetoElectrons, sample, era, isMC):
    
    scaleFactors = []
    if isMC:
        muonRecoSF = getScaleFactor(era, noSel, systName="muon_reco")
        muonIDSF = getScaleFactor(era, noSel, systName="muon_ID")
        muonIsoSF = getScaleFactor(era, noSel, systName="muon_iso")
        scaleFactors = [ muonRecoSF(muons[0]), muonIDSF(muons[0]), muonIsoSF(muons[0]) ]

    oneMuSel = noSel.refine("muon",
                    cut=op.AND(
                        op.rng_len(muons) == 1,
                        #op.rng_len(vetoMuons) == 1,
                        op.rng_len(vetoElectrons) == 0
                    ),
                    weight=scaleFactors
                )
    
    triggerSFWeights = []
    if isMC:
        muonTriggerSF = getScaleFactor(era, oneMuSel, systName="muon_trigger")
        triggerSFWeights.append(muonTriggerSF(muons[0]))
        triggerSFWeights.append(getL1PrefiringSystematic(tree))
    oneMuTriggerSel = oneMuSel.refine("muonTrigger",
                                    cut=muonTriggerDef(tree.HLT, sample, era, isMC),
                                    weight=triggerSFWeights)

    return oneMuSel, oneMuTriggerSel

# Return exclusive electron selection (without and with applied trigger), with trigger, ID/iso, and reco scale factors (for MC)
def buildElectronSelections(tree, noSel, muons, vetoMuons, electrons, vetoElectrons, sample, era, isMC):
    scaleFactors = []
    if isMC:
        eleRecoSF = getScaleFactor(era, noSel, systName="electron_reco")
        eleIDSF = getScaleFactor(era, noSel, systName="electron_ID")
        scaleFactors = [ eleRecoSF(electrons[0]), eleIDSF(electrons[0]) ]

    oneEleSel = noSel.refine("electron",
                    cut=op.AND(
                        op.rng_len(vetoMuons) == 0,
                        op.rng_len(vetoElectrons) == 1,
                        op.rng_len(electrons) == 1
                    ),
                    weight=scaleFactors
                )
    triggerSFWeights = []
    if isMC:
        eleTriggerSF = getScaleFactor(era, oneEleSel, systName="electron_trigger")
        triggerSFWeights.append(eleTriggerSF(electrons[0]))
        triggerSFWeights.append(getL1PrefiringSystematic(tree))
        if era == "2017UL":
            # HLT Z_vtx correction
            triggerSFWeights.append(0.991)
    oneEleTriggerSel = oneEleSel.refine("electronTrigger",
                                    cut=eleTriggerDef(tree.TrigObj, tree.HLT, electrons[0], sample, era, isMC),
                                    weight=triggerSFWeights)

    return oneEleSel, oneEleTriggerSel


bTagWorkingPoints = {
    "2016ULpreVFP": {
        "btagDeepFlavB": {
            "L": 0.0508,
            "M": 0.2598,
            "T": 0.6502
        },
        "btagDeepB": { # DeepCSV
            "L": 0.2027,
            "M": 0.6001,
            "T": 0.8819
        },
    },
    "2016ULpostVFP": {
        "btagDeepFlavB": {
            "L": 0.0480,
            "M": 0.2489,
            "T": 0.6377
        },
        "btagDeepB": { # DeepCSV
            "L": 0.1918,
            "M": 0.5847,
            "T": 0.8767
        },
    },
    "2017UL": {
        "btagDeepFlavB": {
            "L": 0.0532,
            "M": 0.3040,
            "T": 0.7476
        },
        "btagDeepB": {
            "L": 0.1355,
            "M": 0.4506,
            "T": 0.7738
        },
    },
    "2018UL": {
        "btagDeepFlavB": {
            "L": 0.0490,
            "M": 0.2783,
            "T": 0.7100
        },
        "btagDeepB": {
            "L": 0.1208,
            "M": 0.4168,
            "T": 0.7665
        },
    },
}

def bTagDef(jets, era, wp="M", tagger="btagDeepFlavB"):
    
    return op.select(jets, lambda jet: getattr(jet, tagger) >= bTagWorkingPoints[era][tagger][wp])

def METFilter(flags, era, isMC):
    # from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    #
    # from  https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/blob/master/Legacy/event_selection.md
    #   primary vertex filter ("Flag_goodVertices"): data and MC 2016, 2017, 2018.
    #   beam halo filter ("Flag_globalSuperTightHalo2016Filter"): data and MC 2016, 2017, 2018.
    #   HBHE noise filter ("Flag_HBHENoiseFilter"): data and MC 2016, 2017, 2018.
    #   HBHEiso noise filter ("Flag_HBHENoiseIsoFilter"): data and MC 2016, 2017, 2018.
    #   ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter"): data and MC 2016, 2017, 2018.
    #   Bad PF Muon Filter ("Flag_BadPFMuonFilter"): data and MC 2016, 2017, 2018.
    #   ECAL bad calibration filter ("Flag_ecalBadCalibReducedMINIAODFilter"): data and MC 2017, 2018.
    #   ECAL endcap bad SC noise filter ("Flag_eeBadScFilter"): data 2016, 2017, 2018.
    cuts = [
           flags.goodVertices,
           flags.globalSuperTightHalo2016Filter,
           flags.HBHENoiseFilter,
           flags.HBHENoiseIsoFilter,
           flags.EcalDeadCellTriggerPrimitiveFilter,
           flags.BadPFMuonFilter
            ]
#    if "2017" in era or "2018" in era:
#        cuts.append(flags.ecalBadCalibReducedMINIAODFilter) # Only 2017-2018 : both MC and data
        # Could not find in NanoAOD : TODO 
    if not isMC:
        cuts.append(flags.eeBadScFilter) # Only data 

    return cuts

def isBHadron(pId):
        return op.OR((pId // 100) % 10 == 5, op.AND(pId >= 5000, pId <= 5999)) 

class RunEvtSel(NanoAODHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend", type=str, default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import NanoAODDescription
        return super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v5", year="2016", isMC=self.isMC(sample)), backend=self.args.backend)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import CutFlowReport, SummedPlot, Skim
        from bamboo.plots import EquidistantBinning as EqB
        

        plots = []

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight, autoSyst=False)

        era = sampleCfg.get("era") if sampleCfg else None
        cfr = CutFlowReport("yields", recursive=True)
        plots.append(cfr)
        cfr.add(noSel, "Initial")

        plots.append(Plot.make1D(
            "nMus", op.rng_len(t.Muon), noSel,
            EqB(10, -0.5, 9.5), title="Number of all muons (no-selection)"))

        ##### muon selection
        self.muons = op.select(t.Muon, lambda mu: op.AND(mu.pt >= 27, op.abs(mu.eta) <= 2.5, mu.tightId, mu.pfRelIso04_all < 0.15))

        plots.append(Plot.make1D(
            "nGoodMus", op.rng_len(self.muons), noSel,
            EqB(5, -0.5, 4.5), title="Number of good muons (no-selection)"))

        oneMuSel = noSel.refine("oneMuSel", cut=op.AND(op.rng_len(self.muons)==1))

        cfr.add(oneMuSel, "Exactly one muon")

        ##### electron-veto 
        self.electrons = op.select(t.Electron, lambda ele: op.AND(ele.pt > 15, op.abs(ele.eta) < 2.4, op.OR(op.abs(ele.eta + ele.deltaEtaSC) < 1.4442, op.abs(ele.eta + ele.deltaEtaSC) > 1.566), ele.cutBased > 0, ele.pfRelIso03_all < 0.25))

        plots.append(Plot.make1D(
            "nVetoEles", op.rng_len(self.electrons), oneMuSel,
            EqB(5, -0.5, 4.5), title="Number of veto electrons (1#mu-selection)"))

        eleVetoSel = oneMuSel.refine("eleVeto", cut=op.AND(op.rng_len(self.electrons)==0))

        cfr.add(eleVetoSel, "Veto electron")

        ##### jet definition
        self.jets = op.select(t.Jet, lambda j: op.AND(j.jetId == 6, op.abs(j.eta) < 2.4, j.pt > 30., op.deltaR(self.muons[0].p4, j.p4) > 0.4))

        plots.append(Plot.make1D(
            "nJets", op.rng_len(self.jets), eleVetoSel, EqB(10, -0.5, 9.5), title="Number of clean jets (ele_veto-selection)"))
        

        oneMufourJetSel = eleVetoSel.refine("oneMu4JetsSel", cut=[op.rng_len(self.jets) > 3])
        
        cfr.add(oneMufourJetSel, "At least 4 jets")

        plots.append(Plot.make1D(
            "pt_mu", self.muons[0].pt, oneMufourJetSel,
            EqB(100, 20., 120.), title="muon pt (4jet-selection)"))

        plots.append(Plot.make1D(
            "leadJetPt", self.jets[0].pt, oneMufourJetSel,
            EqB(130, 30., 420.), title="leading jet pt (4jet-selection)"))

        plots.append(Plot.make1D(
            "subleadJetPt", self.jets[1].pt, oneMufourJetSel,
            EqB(90, 30., 300.), title="subleading jet pt (4jet-selection)"))

        HT = op.rng_sum(self.jets, lambda j: j.pt)
        plots.append(Plot.make1D(
            "HT", HT, oneMufourJetSel,
            EqB(140, 100., 800.), title="HT (4j-selection)"))

        ##### b-jet definition
        self.bJetsM = op.select(self.jets, lambda bjet: bjet.btagDeepFlavB >= 0.64)
        self.bJetsM = op.sort(self.bJetsM, lambda j: -j.pt)

        plots.append(Plot.make1D(
            "nMediumBJets", op.rng_len(self.bJetsM), oneMufourJetSel, EqB(10, -0.5, 9.5), title="Number of medium b-jets (4jet-selection)"))

        oneMu4Jet1stBSel = oneMufourJetSel.refine("oneMu4Jets1BSel", cut=op.rng_len(self.bJetsM) >= 1)

        cfr.add(oneMu4Jet1stBSel, "At least one b-jets")

        plots.append(Plot.make1D(
            "leadBJetPt_1bSel", self.bJetsM[0].pt, oneMu4Jet1stBSel,
            EqB(130, 30., 420.), title="leading b-jet pt (>= 1b-lection)"))

        plots.append(Plot.make1D(
            "subleadBJetPt_1bSel", self.bJetsM[1].pt, oneMu4Jet1stBSel,
            EqB(90, 30., 300.), title="subleading b-jet pt (>= 1b-selection)"))
        
        oneMu4Jet2ndBSel = oneMu4Jet1stBSel.refine("oneMu4Jets2BSel", cut=op.rng_len(self.bJetsM) >= 2)

        cfr.add(oneMu4Jet2ndBSel, "At least two b-jets")

        plots.append(Plot.make1D(
            "leadBJetPt_2bSel", self.bJetsM[0].pt, oneMu4Jet2ndBSel,
            EqB(130, 30., 420.), title="leading b-jet pt (>= 2b-lection)"))

        plots.append(Plot.make1D(
            "subleadBJetPt_2bSel", self.bJetsM[1].pt, oneMu4Jet2ndBSel,
            EqB(90, 30., 300.), title="subleading b-jet pt (>= 2b-selection)"))

        ##### [try to do rec-soft-muons definition, 02 of march], should be moved to the next class once worked 

        genSoftMusFromBHadrons = op.select(t.GenPart, lambda mu: op.AND(op.abs(mu.pdgId) == 13, isBHadron(op.abs(mu.genPartMother.pdgId)), op.OR(op.abs(mu.genPartMother.genPartMother.pdgId) == 5 , isBHadron(op.abs(mu.genPartMother.genPartMother.pdgId))), mu.statusFlags & (0x1 << 6), op.NOT(mu.statusFlags & (0x1 << 2)), op.NOT(mu.statusFlags & (0x1 << 3)), mu.status == 1, op.rng_any(mu.ancestors, lambda a : op.abs(a.pdgId) == 6) ) )

        plots.append(Plot.make1D(
            "nGenSoftMus", op.rng_len(genSoftMusFromBHadrons), oneMu4Jet2ndBSel, EqB(10, 0., 10.), title="n gen soft mus from b-hadron decays"))

        self.soft_muons_preSel = op.select(t.Muon, lambda mu: op.AND(mu.pt >= 10, op.abs(mu.eta) <= 2.4, mu.looseId, mu.dxy < 0.3, mu.dz < 20, mu.nTrackerLayers > 5, mu.idx != self.muons[0].idx, op.rng_any(self.bJetsM, lambda bjet: bjet.idx == mu.jet.idx)))

        oneSoftMu_preSel = oneMu4Jet2ndBSel.refine("muon_4jets_2b_1softMu", cut=op.rng_len(self.soft_muons_preSel) >= 1)
        cfr.add(oneSoftMu_preSel, "with one soft mu (pre-sel)")

        plots.append(Plot.make1D(
            "DR_softMu_associatedJet", op.deltaR(t.Jet[self.soft_muons_preSel[0].jet.idx].p4,self.soft_muons_preSel[0].p4), oneSoftMu_preSel, EqB(50, 0., 5.), title="DR_softMu_associatedJet"))

        if op.AND(op.rng_len(genSoftMusFromBHadrons)==1, op.deltaR(genSoftMusFromBHadrons[0].p4,self.soft_muons_preSel[0].p4)<0.3):
         plots.append(Plot.make1D(
            "matchedsoftMu_pt_ratio", self.soft_muons_preSel[0].pt/t.Jet[self.soft_muons_preSel[0].jet.idx].pt, oneSoftMu_preSel, EqB(50, 0., 5.), title="soft_pt_ratio"))

        self.soft_muons = op.select(self.soft_muons_preSel, lambda mu: op.rng_any(self.bJetsM, lambda bjet: op.AND(bjet.idx == mu.jet.idx, mu.pt/t.Jet[bjet.idx].pt>0.3)))
        oneSoftMu = oneSoftMu_preSel.refine("finalSel", cut=op.rng_len(self.soft_muons) >= 1)
        cfr.add(oneSoftMu, "with one soft mu")

        plots.append(Plot.make1D(
            "nRecSoftMus", op.rng_len(self.soft_muons), oneSoftMu, EqB(10, 0., 10.), title="n rec soft mus from b-hadron decays"))

        return plots


class calculateBTagEff(NanoAODHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend", type=str, default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import NanoAODDescription
        return super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v5", year="2016", isMC=self.isMC(sample)), backend=self.args.backend)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo.plots import VariableBinning as VarBin
        

        era = sampleCfg.get("era") if sampleCfg else None
        ## also need to take into account TTbar reweighting and pre-firing rates
        # https://gitlab.cern.ch/gsaha/hhbbww_sl/-/blob/master/BaseHHtobbWW.py?ref_type=heads#L462
        # https://gitlab.cern.ch/gsaha/hhbbww_sl/-/blob/master/BaseHHtobbWW.py?ref_type=heads#L481
        ##noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight, autoSyst=False)
        # MET filter #
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, self.isMC(sample)) )


        noSel = noSel.refine("PV", cut=op.AND(t.PV.npvsGood>=1, t.PV.ndof>5))

        #--------------------------------------------
        #--- Muon Selection & Electron Veto ---------
        #--------------------------------------------

        self.muons = op.select(t.Muon, muonDef(era))
        self.electrons = op.select(t.Electron, eleDef(era))

        self.vetoMuons = op.select(t.Muon, vetoMuonDef)
        self.vetoElectrons = op.select(t.Electron, vetoEleDef)

        _,oneMuTriggerSel = buildMuonSelections(t, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))
        #_,oneEleTriggerSel = buildElectronSelections(t, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))

        plots.append(Plot.make1D(
            "pt_mu_1MuSel", self.muons[0].pt, oneMuTriggerSel,
            EqB(100, 20., 120.), title="muon pt (1#mu-selection)"))
        plots.append(Plot.make1D("PV_npvsGood", t.PV.npvsGood, oneMuTriggerSel, EqB(60, 0, 60), title="Number of good primary vertices (1#mu-selection)"))
        #--------------------------------------------
        #--- Jet Selection --------------------------
        #--------------------------------------------

        self.rawJets = op.select(t.Jet, jetDef)
        self.cleanedJets =cleanJets(self.rawJets, self.muons, self.electrons)

        plots.append(Plot.make1D(
            "nJets", op.rng_len(self.cleanedJets), oneMuTriggerSel, EqB(10, 0., 10.), title="Number of clean jets (1#mu-selection)"))

        oneMufourJetSel = oneMuTriggerSel.refine("oneMuonsFourJets", cut=[op.rng_len(self.cleanedJets) > 3])

        plots.append(Plot.make1D("PV_npvsGood_noPUweight", t.PV.npvsGood, oneMufourJetSel, EqB(60, 0, 60), title="Number of good primary vertices (4j-selection, w/o PUw)"))
        if self.isMC(sample):
            pileup = makePUWeight(t, era, oneMufourJetSel)
            oneMufourJetSel = oneMufourJetSel.refine("weights", weight=pileup)
        plots.append(Plot.make1D("PV_npvsGood_withPUweight", t.PV.npvsGood, oneMufourJetSel, EqB(60, 0, 60), title="Number of good primary vertices (4j-selection, w PUw)"))
        
        plots.append(Plot.make1D(
            "dxy_mu_4JSel", op.abs(self.muons[0].dxy), oneMufourJetSel,
            EqB(20, 0., 0.01), title="muon dxy (4j-selection)"))
        plots.append(Plot.make1D(
            "pt_mu_4JSel", self.muons[0].pt, oneMufourJetSel,
            EqB(65, 20., 150.), title="muon pt (4j-selection)"))
        leadjpt = Plot.make1D(
            "leadJetPT", self.cleanedJets[0].pt, oneMufourJetSel, EqB(60, 0., 300.), title="Leading jet pt")
        # leadjphi = Plot.make1D(
        #     "leadJetPHI", op.Phi_mpi_pi(jets[0].phi), oneMufourJetSel, EqB(50, -3.142, 3.142),
        #     title="Leading jet PHI")
        subleadjpt = Plot.make1D(
            "subleadJetPT", self.cleanedJets[1].pt, oneMufourJetSel, EqB(60, 0., 300.), title="subleading jet pt")
        plots += [leadjpt, subleadjpt]
#        plots.append(SummedPlot(
#            "twoLeadJetPT", [leadjpt, subleadjpt], xTitle="Leading two jet PTs"))

        HT = op.rng_sum(self.cleanedJets, lambda j: j.pt)
        plots.append(Plot.make1D(
            "HT", HT, oneMufourJetSel,
            EqB(140, 100., 800.), title="HT (4j-selection)"))

        #--------------------------------------------
        #--- B-Jet Efficiency finder ----------------
        #--------------------------------------------

        ### idea from: https://mattermost.web.cern.ch/cms-exp/pl/ic9bmru7n3gjxjqqc1nfus576o, also https://gitlab.cern.ch/swertz/ttbbRun2Bamboo/-/blob/master/python/bTagPlotter.py?ref_type=heads#L93
        bottomFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 5)
        charmFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 4)
        lightFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 0)
 
        for flav, flavJets in zip(['b', 'c', 'light'], [bottomFlavJets, charmFlavJets, lightFlavJets]):
          # b tagging efficiencies as a function of flavour/pt/|eta|
          binning = (VarBin([30, 40, 60, 80, 100, 200, 350, 1000]), EqB(5, 0, 2.5))

          pt = op.map(flavJets, lambda j: j.pt)
          eta = op.map(flavJets, lambda j: op.abs(j.eta))
          plots.append(Plot.make2D(f"Denominator_jet_pt_eta_{flav}", (pt, eta), oneMufourJetSel, binning, xTitle="Jet pt", yTitle="Jet eta"))

          for wp in ["L", "M", "T"]:
              bTagThr = bTagWorkingPoints[era]["btagDeepFlavB"][wp]
              selJets = op.select(flavJets, lambda j: j.btagDeepFlavB >= bTagThr)
              pt = op.map(selJets, lambda j: j.pt)
              eta = op.map(selJets, lambda j: op.abs(j.eta))
              plots.append(Plot.make2D(f"Numerator_jet_pt_eta_{flav}_wp{wp}", (pt, eta), oneMufourJetSel, binning, xTitle="Jet pt", yTitle="Jet eta"))

        return plots

class MinimalTopMass(NanoAODHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend", type=str, default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import NanoAODDescription
        return super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v5", year="2016", isMC=self.isMC(sample)), backend=self.args.backend)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import CutFlowReport, SummedPlot, Skim
        from bamboo.plots import EquidistantBinning as EqB
        

        era = sampleCfg.get("era") if sampleCfg else None
        ## also need to take into account TTbar reweighting and pre-firing rates
        # https://gitlab.cern.ch/gsaha/hhbbww_sl/-/blob/master/BaseHHtobbWW.py?ref_type=heads#L462
        # https://gitlab.cern.ch/gsaha/hhbbww_sl/-/blob/master/BaseHHtobbWW.py?ref_type=heads#L481
        ##noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

        cfr = CutFlowReport("yields", recursive=True)
        plots.append(cfr)
        cfr.add(noSel, "Initial")

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight, autoSyst=True)
        # MET filter #
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, self.isMC(sample)) )
        
        psISRSyst = op.systematic(op.c_float(1.), name="psISR", up=t.PSWeight[0], down=t.PSWeight[2])
        noSel = noSel.refine("psISR", weight=psISRSyst)
        mcWgts = []
        mcWgts += [op.systematic(op.c_float(1.), **{f"qcdScale{i:d}": t.LHEScaleWeight[i] for i in (0, 1, 3, 5, 7, 8)}),]
        noSel = noSel.refine("muFmuRScaleWeight", weight=mcWgts)
        
#       if self.isMC(sample):
#            if self.doJES =="merged":
#                jesUncertaintySources = ['Regrouped_Absolute', f'Regrouped_Absolute_{era_}',
#                                         'Regrouped_BBEC1', f'Regrouped_BBEC1_{era_}',
#                                         'Regrouped_EC2', f'Regrouped_EC2_{era_}',
#                                         'Regrouped_FlavorQCD',
#                                         'Regrouped_HF', f'Regrouped_HF_{era_}',
#                                         'Regrouped_RelativeBal', f'Regrouped_RelativeSample_{era_}']
#            elif self.doJES == "total":
#                jesUncertaintySources = ["Total"]
#            elif self.doJES == "all":
        jesUncertaintySources = ["Total"]

        JECs = {'2016-preVFP' : "Summer19UL16APV_V7_MC",
                    '2016-postVFP': "Summer19UL16_V7_MC",
                    '2017'        : "Summer19UL17_V5_MC",
                    '2018UL'        : "Summer19UL18_V5_MC"
                    }

        JERs = {'2016-preVFP' : "Summer20UL16APV_JRV3_MC",
                    '2016-postVFP': "Summer20UL16_JRV3_MC",
                    '2017'        : "Summer19UL17_JRV2_MC",
                    '2018UL'        : "Summer19UL18_JRV2_MC"
                    }
#        else:
#            jesUncertaintySources = None
#            JECs = {'2016-preVFP' : "Summer19UL16APV_RunBCDEF_V7_DATA",
#                    '2016-postVFP': "Summer19UL16_RunFGH_V7_DATA",
#                    '2017'        : "Summer19UL17_RunBCDEF_V5_DATA",
#                    '2018'        : "Summer19UL18_V5_DATA",
#                    }
#
#            JERs = {'2016-preVFP' : "Summer20UL16APV_JRV3_DATA",
#                    '2016-postVFP': "Summer20UL16_JRV3_DATA",
#                    '2017'        : "Summer19UL18_JRV2_DATA",
#                    '2018'        : "Summer19UL18_JRV2_DATA",
#                    }
#
        cmJMEArgs = {
                "jsonFile": localizePOGSF(era, "JME", "jet_jerc.json.gz"),
                "jsonFileSmearingTool": os.path.join("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration", "POG", "JME", "jer_smear.json.gz"),
                "jec": JECs[era],
                "smear": JERs[era],
                #"splitJER": False,
                "jesUncertaintySources": jesUncertaintySources,
                "addHEM2018Issue": (era == "2018"),
                "isMC": self.isMC(sample),
                "backend": self.args.backend,
                #"jecLevels":[], #  default : L1FastJet, L2Relative, L3Absolute, and also L2L3Residual for data
                }

        # configure corrections and variations
        from bamboo.analysisutils import configureJets
        configureJets(t.Jet, jetType="AK4PFchs", **cmJMEArgs)


        noSel = noSel.refine("PV", cut=op.AND(t.PV.npvsGood>=1, t.PV.ndof>5))

        #--------------------------------------------
        #--- Muon Selection & Electron Veto ---------
        #--------------------------------------------

        self.muons = op.select(t.Muon, muonDef(era))
        self.electrons = op.select(t.Electron, eleDef(era))

        self.vetoMuons = op.select(t.Muon, vetoMuonDef)
        self.vetoElectrons = op.select(t.Electron, vetoEleDef)

        _,oneMuTriggerSel = buildMuonSelections(t, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))
        #_,oneEleTriggerSel = buildElectronSelections(t, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))

        cfr.add(oneMuTriggerSel, "Lepton selection (1 muon, ele veto)")
        
        plots.append(Plot.make1D(
            "pt_mu_1MuSel", self.muons[0].pt, oneMuTriggerSel,
            EqB(100, 20., 120.), title="muon pt (1#mu-selection)"))
        plots.append(Plot.make1D("PV_npvsGood", t.PV.npvsGood, oneMuTriggerSel, EqB(60, 0, 60), title="Number of good primary vertices (1#mu-selection)"))
        #--------------------------------------------
        #--- Jet Selection --------------------------
        #--------------------------------------------

        self.rawJets = op.select(t.Jet, jetDef)
        self.cleanedJets =cleanJets(self.rawJets, self.muons, self.electrons)

##@@        moved this plot to be comparable with the plots after b-bweight
##@@        plots.append(Plot.make1D(
##@@            "nJets", op.rng_len(self.cleanedJets), oneMuTriggerSel, EqB(10, 0., 10.), title="Number of clean jets (1#mu-selection)"))

        plots.append(Plot.make1D(
            "nJets_noPUweight", op.rng_len(self.cleanedJets), oneMuTriggerSel, EqB(10, 0., 10.), title="Number of clean jets (1#mu-selection)"))

        oneMufourJetSel = oneMuTriggerSel.refine("oneMu4Jets", cut=[op.rng_len(self.cleanedJets) > 3])
        cfr.add(oneMufourJetSel, "At least 4 jets")

        plots.append(Plot.make1D("PV_npvsGood_noPUweight", t.PV.npvsGood, oneMufourJetSel, EqB(60, 0, 60), title="Number of good primary vertices (4j-selection, w/o PUw)"))
        if self.isMC(sample):
            pileup = makePUWeight(t, era, oneMufourJetSel)
            oneMufourJetSel = oneMufourJetSel.refine("puWeight", weight=pileup)
        plots.append(Plot.make1D("PV_npvsGood_withPUweight", t.PV.npvsGood, oneMufourJetSel, EqB(60, 0, 60), title="Number of good primary vertices (4j-selection, w PUw)"))
        plots.append(Plot.make1D("rho", t.fixedGridRhoFastjetCentralChargedPileUp, oneMufourJetSel, EqB(60, 0, 60), title="#rho from charged PF Candidates for central region"))
        
        plots.append(Plot.make1D(
            "dxy_mu_4JSel", op.abs(self.muons[0].dxy), oneMufourJetSel,
            EqB(20, 0., 0.01), title="muon dxy (4j-selection)"))
        plots.append(Plot.make1D(
            "pt_mu_4JSel", self.muons[0].pt, oneMufourJetSel,
            EqB(65, 20., 150.), title="muon pt (4j-selection)"))
        leadjpt = Plot.make1D(
            "leadJetPT", self.cleanedJets[0].pt, oneMufourJetSel, EqB(78, 30., 420.), title="leading jet pt")
        # leadjphi = Plot.make1D(
        #     "leadJetPHI", op.Phi_mpi_pi(jets[0].phi), oneMufourJetSel, EqB(50, -3.142, 3.142),
        #     title="Leading jet PHI")
        subleadjpt = Plot.make1D(
            "subleadJetPT", self.cleanedJets[1].pt, oneMufourJetSel, EqB(54, 30., 300.), title="subleading jet pt")
        plots += [leadjpt, subleadjpt]
#        plots.append(SummedPlot(
#            "twoLeadJetPT", [leadjpt, subleadjpt], xTitle="Leading two jet PTs"))

        HT = op.rng_sum(self.cleanedJets, lambda j: j.pt)
        plots.append(Plot.make1D(
            "HT", HT, oneMufourJetSel,
            EqB(70, 100., 800.), title="HT (4j-selection)"))

        #--------------------------------------------
        #--- B-Jet Selection ------------------------
        #--------------------------------------------
        self.bJetsM = bTagDef(self.cleanedJets, era, "M", "btagDeepFlavB")


##===========================
##      this piece was working properly for the talk on 16 march!
##        oneMufourJets2MBSel = oneMufourJetSel.refine("oneMu4Jets2B", cut=op.rng_len(self.bJetsM) >= 2)
##        cfr.add(oneMufourJets2MBSel, "At least 2 b-jets")
##
##        plots.append(Plot.make1D(
##            "nMBJets", op.rng_len(self.bJetsM), oneMufourJets2MBSel, EqB(10, 0., 10.), title="Number of M b-jets"))
##
##        plots.append(Plot.make1D(
##            "leadBJetPt", self.bJetsM[0].pt, oneMufourJets2MBSel,
##            EqB(78, 30., 420.), title="leading b-jet pt (>= 1b-lection)"))
##
##        plots.append(Plot.make1D(
##            "subleadBJetPt_1bSel", self.bJetsM[1].pt, oneMufourJets2MBSel,
##            EqB(54, 30., 300.), title="subleading b-jet pt (>= 1b-selection)"))
##
##        from bamboo.scalefactors import get_correction
##        params = {
##         "pt": lambda j: j.pt, "abseta": lambda j: op.abs(j.eta), "working_point": "M", "flavor": 5}
##        bjetSFs = get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet_comb", params=params,
##                              systParam="systematic", systNomName="central", sel=oneMufourJets2MBSel,
##                             defineOnFirstUse=True)
##        b_jet1_sf = bjetSFs(self.bJetsM[0])
##        b_jet2_sf = bjetSFs(self.bJetsM[1])
##        if self.isMC(sample): 
##            oneMufourJets2MBSel = oneMufourJets2MBSel.refine("btag_sf", weight= [b_jet1_sf, b_jet2_sf])

##        plots.append(Plot.make1D(
##            "nMBJets_reweighted", op.rng_len(self.bJetsM), oneMufourJets2MBSel, EqB(10, 0., 10.), title="Number of M b-jets (reweighted with SF_b1*SF_b2)"))

##        plots.append(Plot.make1D(
##            "leadBJetPT_reweighted", self.bJetsM[0].pt, oneMufourJets2MBSel, EqB(78, 30., 420.), title="Leading B jet PT (reweighted with SF_b1*SF_b2)"))

##        plots.append(Plot.make1D(
##            "subleadBJetPT_reweighted", self.bJetsM[1].pt, oneMufourJets2MBSel, EqB(54, 30., 300.), title="subleading B jet PT (reweighted with SF_b1*SF_b2)"))

##===========================

        oneMufourJets2MBSel = oneMufourJetSel.refine("sf1sf2oneMu4Jets2B", cut=op.rng_len(self.bJetsM) >= 2)

        from bamboo.scalefactors import get_bTagSF_fixWP, get_correction, makeBtagWeightMeth1a

        params = {
         "pt": lambda j: j.pt, "abseta": lambda j: op.abs(j.eta), "working_point": "M", "flavor": 5}
        bjetSFs = get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet_comb", params=params,
                              systParam="systematic", systNomName="central", sel=oneMufourJets2MBSel,
                             defineOnFirstUse=True)
        b_jet1_sf = bjetSFs(self.bJetsM[0])
        b_jet2_sf = bjetSFs(self.bJetsM[1])
        if self.isMC(sample): 
            oneMufourJets2MBSel= oneMufourJets2MBSel.refine("btag_sfMultiplication", weight= [b_jet1_sf, b_jet2_sf])        
        cfr.add(oneMufourJets2MBSel, "weighted_sf1_times_sf2")
        plots.append(Plot.make1D(
            "nBJets_sf1_times_sf2", op.rng_len(self.bJetsM), oneMufourJets2MBSel, EqB(10, 0., 10.), title="Number of b jets (weighted_sf1_times_sf2)"))
#####
#####        
#####        trueWeighted1Mu4J2B = oneMufourJetSel.refine("trueWeightoneMu4Jets2B", cut=op.rng_len(self.bJetsM) >= 2)
#####
#####        plots.append(Plot.make1D(
#####            "nBJets_doubleCheck", op.rng_len(self.bJetsM), trueWeighted1Mu4J2B, EqB(10, 0., 10.), title="Number of b jets"))


####        from bamboo.scalefactors import get_bTagSF_fixWP, get_correction, makeBtagWeightMeth1a

        trueWeighted1Mu4J2B = oneMufourJetSel.refine("trueWeightoneMu4Jets2B", cut=op.rng_len(self.bJetsM) >= 2)
        trueWeighted1Mu4J2B_sw = oneMufourJetSel.refine("trueWeightoneMu4Jets2B_sw", cut=op.rng_len(self.bJetsM) >= 2)

        bTagEff_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "../", "BTaggEffProduction_allSignalStats", "FlavIncludedBTagEff_maps_UL18.json.gz")

        effparams = {
            "pt": lambda j: j.pt, "eta": lambda j: op.abs(j.eta), "flavor": lambda j: j.hadronFlavour#, "workingPoint": "M"
        }
        btvEff = {"M": get_correction(bTagEff_file, "bTaggingEffValues", params=effparams, sel=trueWeighted1Mu4J2B, defineOnFirstUse=True)}
        btvSF = lambda wp, flav: get_bTagSF_fixWP(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet", wp, flav, trueWeighted1Mu4J2B, syst_prefix="btagSF_", decorr_eras=False, era="2018_UL", decorr_wps=False)

        ### !!!! BEGIN check the weight which is calculated based on swerts eff file !!!! ### 
        bTagEff_file_sw = os.path.join(os.path.dirname(os.path.dirname(__file__)), "../", "btagEff_deepJet_2018.json")
        effparams_sw = {
            "pt": lambda j: j.pt, "eta": lambda j: j.eta,
            "jetFlavor": lambda j: j.hadronFlavour, "workingPoint": "M"
        }
        btvEff_sw = {"M": get_correction(bTagEff_file_sw, "btagEff", params=effparams_sw, sel=trueWeighted1Mu4J2B_sw, defineOnFirstUse=True)}
        btvSF_sw = lambda wp, flav: get_bTagSF_fixWP(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet", wp, flav, trueWeighted1Mu4J2B_sw, syst_prefix="btagSF_", decorr_eras=False, era="2018_UL", decorr_wps=False)
        ### !!!! END check the weight which is calculated based on swerts eff file !!!! ### 

        btvWeight = -1.
        btvWeight_swertz = -1.
        if self.isMC(sample): 
          btvWeight = makeBtagWeightMeth1a(self.cleanedJets, "btagDeepFlavB", ["M"], {"M": 0.2783}, btvSF, btvEff)
          btvWeight_swertz = makeBtagWeightMeth1a(self.cleanedJets, "btagDeepFlavB", ["M"], {"M": 0.2783}, btvSF_sw, btvEff_sw)
          plots.append(Plot.make1D(
            "btvWeight", btvWeight, trueWeighted1Mu4J2B, EqB(60, 0.8, 1.4), title="btvWeight"))
          plots.append(Plot.make1D(
            "btvWeight_sw", btvWeight_swertz, trueWeighted1Mu4J2B_sw, EqB(60, 0.8, 1.4), title="btvWeight_sw"))
          trueWeighted1Mu4J2B = trueWeighted1Mu4J2B.refine("btag_trueWeight", weight=btvWeight)
          trueWeighted1Mu4J2B_sw = trueWeighted1Mu4J2B_sw.refine("btag_trueWeight_sw", weight=btvWeight_swertz)
        cfr.add(trueWeighted1Mu4J2B, "weighted_true_makeBtagWeightMeth1a")
        cfr.add(trueWeighted1Mu4J2B_sw, "weighted_true_makeBtagWeightMeth1a_sw")
        plots.append(Plot.make1D(
            "nJets", op.rng_len(self.cleanedJets), oneMufourJetSel, EqB(10, 0., 10.), title="Number of clean jets (no weight applied)"))
        plots.append(Plot.make1D(
            "nJets_weighted", op.rng_len(self.cleanedJets), trueWeighted1Mu4J2B, EqB(10, 0., 10.), title="Number of jets (with at least 2b)"))
        plots.append(Plot.make1D(
            "nJets_weighted_sw", op.rng_len(self.cleanedJets), trueWeighted1Mu4J2B_sw, EqB(10, 0., 10.), title="Number of jets (with at least 2b)"))
        plots.append(Plot.make1D(
            "nBJets_weighted", op.rng_len(self.bJetsM), trueWeighted1Mu4J2B, EqB(10, 0., 10.), title="Number of B jets"))
        plots.append(Plot.make1D(
            "nBJets_weighted_sw", op.rng_len(self.bJetsM), trueWeighted1Mu4J2B_sw, EqB(10, 0., 10.), title="Number of B jets"))

        if self.isMC(sample): 
         plots.append(Plot.make1D(
            "selectedBs_hadronFlavs", self.bJetsM[op.rng_len(self.bJetsM)-1].hadronFlavour, trueWeighted1Mu4J2B, EqB(10, 0., 10.), title="b-jet hadron flavor"))
         plots.append(Plot.make1D(
            "selectedBs_hadronFlavs_sw", self.bJetsM[op.rng_len(self.bJetsM)-1].hadronFlavour, trueWeighted1Mu4J2B_sw, EqB(10, 0., 10.), title="b-jet hadron flavor"))


        plots.append(Plot.make1D(
            "nMediumBJets", op.rng_len(self.bJetsM), oneMufourJetSel, EqB(10, 0., 10.), title="Number of medium b-jets (no weight applied)"))

        plots.append(Plot.make1D(
            "leadBJetPt_2bSel", self.bJetsM[0].pt, trueWeighted1Mu4J2B,
            EqB(130, 30., 420.), title="leading b-jet pt"))

        plots.append(Plot.make1D(
            "subleadBJetPt_2bSel", self.bJetsM[1].pt, trueWeighted1Mu4J2B,
            EqB(90, 30., 300.), title="subleading b-jet pt"))


        #--------------------------------------------
        #--- (gen)Soft Muon Selection ---------------
        #--------------------------------------------

        if self.isMC(sample): 
         genSoftMusFromBHadrons = op.select(t.GenPart, lambda mu: op.AND(op.abs(mu.pdgId) == 13, isBHadron(op.abs(mu.genPartMother.pdgId)), op.OR(op.abs(mu.genPartMother.genPartMother.pdgId) == 5 , isBHadron(op.abs(mu.genPartMother.genPartMother.pdgId))), mu.statusFlags & (0x1 << 6), op.NOT(mu.statusFlags & (0x1 << 2)), op.NOT(mu.statusFlags & (0x1 << 3)), mu.status == 1, op.rng_any(mu.ancestors, lambda a : op.abs(a.pdgId) == 6) ) )

         plots.append(Plot.make1D(
            "nGenSoftMus", op.rng_len(genSoftMusFromBHadrons), trueWeighted1Mu4J2B, EqB(10, 0., 10.), title="#mu_{soft}^{gen}"))


        #--------------------------------------------
        #--- (rec)Soft Muon Selection ---------------
        #--------------------------------------------


        self.soft_muons_preSel = op.select(t.Muon, lambda mu: op.AND(mu.pt >= 10, op.abs(mu.eta) <= 2.4, mu.looseId, mu.dxy < 0.3, mu.dz < 20, mu.nTrackerLayers > 5, mu.idx != self.muons[0].idx, op.rng_any(self.bJetsM, lambda bjet: bjet.idx == mu.jet.idx)))

        oneSoftMu_preSel = trueWeighted1Mu4J2B.refine("muon_4jets_2b_1softMu", cut=op.rng_len(self.soft_muons_preSel) >= 1)
        cfr.add(oneSoftMu_preSel, "with one soft mu (pre-sel)")

        plots.append(Plot.make1D(
            "DR_softMu_associatedJet", op.deltaR(t.Jet[self.soft_muons_preSel[0].jet.idx].p4,self.soft_muons_preSel[0].p4), oneSoftMu_preSel, EqB(20, 0., 0.5), title="DR_softMu_associatedJet"))
#
#        if self.isMC(sample):
#         #if op.rng_len(genSoftMusFromBHadrons)==1:
#          oneSoftMu_preSel_withoneGenMu = oneSoftMu_preSel.refine("withoneGenMu", cut=op.rng_len(genSoftMusFromBHadrons) == 1)
#          plots.append(Plot.make1D(
#            "DR_gensoftMu_recsoftMuPresel", op.deltaR(genSoftMusFromBHadrons[0].p4,self.soft_muons_preSel[0].p4), oneSoftMu_preSel_withoneGenMu, EqB(50, 0., 5.), title="DR_gensoftMu_recsoftMuPresel"))
#
#
#         #if op.AND(op.rng_len(genSoftMusFromBHadrons)==1, op.deltaR(genSoftMusFromBHadrons[0].p4,self.soft_muons_preSel[0].p4)<0.05):
#          oneSoftMu_preSel_withoneGenMu_matched = oneSoftMu_preSel_withoneGenMu.refine("withoneGenMu_matched", cut=op.deltaR(genSoftMusFromBHadrons[0].p4,self.soft_muons_preSel[0].p4)<0.05)
#    
#          plots.append(Plot.make1D(
#            "matchedsoftMu_pt_ratio", self.soft_muons_preSel[0].pt/t.Jet[self.soft_muons_preSel[0].jet.idx].pt, oneSoftMu_preSel_withoneGenMu_matched, EqB(50, 0., 1.), title="soft_pt_ratio"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_DR_softMu_associatedJet",op.deltaR(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4), oneSoftMu_preSel_withoneGenMu_matched, EqB(20, 0., 0.5), title="DR_softMu_associatedJet(events with matched gen-soft mu)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_Dphi_softMu_associatedJet",op.deltaPhi(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4), oneSoftMu_preSel_withoneGenMu_matched, EqB(128, -3.2, 3.2), title="Dphi_softMu_associatedJet(events with matched gen-soft mu)"))
#
#          dEta = op.sqrt(op.pow(op.deltaR(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4),2)-op.pow(op.deltaPhi(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4),2))
#          dTheta = 2.*op.atan(op.exp(-1.*dEta))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_Dtheta_softMu_associatedJet", dTheta, oneSoftMu_preSel_withoneGenMu_matched, EqB(64, 0., 3.2), title="Dtheta_softMu_associatedJet(events with matched gen-soft mu)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_cosDphi_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.cos(op.deltaPhi(self.soft_muons_preSel[0].p4, t.Jet[self.soft_muons_preSel[0].jet.idx].p4))), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 100.), title="point product (d-phi)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_sinDphi_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.abs(op.sin(op.deltaPhi(self.soft_muons_preSel[0].p4, t.Jet[self.soft_muons_preSel[0].jet.idx].p4)))), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 20.), title="cross product (d-phi)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_cosDtheta_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.cos(dTheta)), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 20.), title="point product (d-theta)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_sinDtheta_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.abs(op.sin(dTheta))), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 100.), title="cross product (d-theta)"))
#
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_cosDR_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.cos(op.deltaR(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4))), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 100.), title="point product (d-R)"))
#
#          plots.append(Plot.make1D(
#            "matchedsoftMu_sinDR_mu_associatedJet", op.product(self.soft_muons_preSel[0].pt, op.abs(op.sin(op.deltaR(self.soft_muons_preSel[0].p4,t.Jet[self.soft_muons_preSel[0].jet.idx].p4)))), oneSoftMu_preSel_withoneGenMu_matched, EqB(200, 0., 100.), title="cross product (d-R)"))
#

        self.soft_muons_ptRel = op.select(self.soft_muons_preSel, lambda mu: op.rng_any(self.bJetsM, lambda bjet: op.AND(bjet.idx == mu.jet.idx, op.product(mu.pt, op.cos(2.*op.atan(op.exp(-1.*op.sqrt(op.pow(op.deltaR(mu.p4,t.Jet[bjet.idx].p4),2)-op.pow(op.deltaPhi(mu.p4,t.Jet[bjet.idx].p4),2))))))>0.3)))

        oneSoftMu_ptRel = oneSoftMu_preSel.refine("finalSel_ptRel", cut=op.rng_len(self.soft_muons_ptRel) >= 1)
        cfr.add(oneSoftMu_ptRel, "one soft mu (ptRel)")

        plots.append(Plot.make1D(
            "nRecSoftMus_ptRel", op.rng_len(self.soft_muons_ptRel), oneSoftMu_ptRel, EqB(10, 0., 10.), title="#mu_{soft}^{rec} (ptRel)"))

        plots.append(Plot.make1D(
            "MlMu_ptRel", op.invariant_mass(self.muons[0].p4, self.soft_muons_ptRel[0].p4), oneSoftMu_ptRel, EqB(25, 0., 250.), title="M_{#mu^{hard},#mu^{soft}} OR M_{l#mu}"))


        self.soft_muons = op.select(self.soft_muons_preSel, lambda mu: op.rng_any(self.bJetsM, lambda bjet: op.AND(bjet.idx == mu.jet.idx, mu.pt/t.Jet[bjet.idx].pt>0.18)))

        oneSoftMu = oneSoftMu_preSel.refine("finalSel", cut=op.rng_len(self.soft_muons) >= 1)
        cfr.add(oneSoftMu, "with one soft mu")

        plots.append(Plot.make1D(
            "nRecSoftMus_ptRatio", op.rng_len(self.soft_muons), oneSoftMu, EqB(10, 0., 10.), title="#mu_{soft}^{rec} (pt ratio)"))

        plots.append(Plot.make1D(
            "MlMu_ptRatio", op.invariant_mass(self.muons[0].p4, self.soft_muons[0].p4), oneSoftMu, EqB(25, 0., 250.), title="M_{#mu^{hard},#mu^{soft}} OR M_{l#mu}"))

        plots.append(Skim(
                "eventSel_Skim", {
                    "nSoftMuons": op.static_cast("UInt_t", op.rng_len(self.soft_muons)),  # TBranch doesn't accept size_t
                    "softMuons_i": self.soft_muons.idxs
                }, oneSoftMu,
                keepOriginal=[Skim.KeepAll]))


        return plots

class test_btagWeight_swertz(NanoAODHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend", type=str, default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import NanoAODDescription
        return super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v5", year="2016", isMC=self.isMC(sample)), backend=self.args.backend)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        

        era = sampleCfg.get("era") if sampleCfg else None

        plots = []

        cfr = CutFlowReport("yields", recursive=True)
        plots.append(cfr)
        cfr.add(noSel, "Initial")

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight, autoSyst=False)
        self.rawJets = op.select(t.Jet, jetDef)

        fourJetSel = noSel.refine("fourJets", cut=[op.rng_len(self.rawJets) > 3])
        cfr.add(fourJetSel, "At least 4 jets")

        plots.append(Plot.make1D(
            "nJets", op.rng_len(self.rawJets), fourJetSel, EqB(10, 0., 10.), title="Number of raw jets (no-selection)"))

        self.bJetsM = bTagDef(self.rawJets, era, "M", "btagDeepFlavB")
        fourJet2BSel = fourJetSel.refine("fourJets2B", cut=[op.rng_len(self.bJetsM) > 1])
        cfr.add(fourJet2BSel, "At least two b-jets")

        plots.append(Plot.make1D(
            "nBJets", op.rng_len(self.bJetsM), fourJet2BSel, EqB(10, 0., 10.), title="Number of b jets"))


        from bamboo.scalefactors import get_bTagSF_fixWP, get_correction, makeBtagWeightMeth1a

        params = {
         "pt": lambda j: j.pt, "abseta": lambda j: op.abs(j.eta), "working_point": "M", "flavor": 5}
        bjetSFs = get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet_comb", params=params,
                              systParam="systematic", systNomName="central", sel=fourJet2BSel,
                             defineOnFirstUse=True)
        b_jet1_sf = bjetSFs(self.bJetsM[0])
        b_jet2_sf = bjetSFs(self.bJetsM[1])
        if self.isMC(sample): 
            fourJet2BSel_sf_weighted = fourJet2BSel.refine("btag_sf", weight= [b_jet1_sf, b_jet2_sf])        
        cfr.add(fourJet2BSel_sf_weighted, "weighted_sf1_times_sf2")
        plots.append(Plot.make1D(
            "nBJets_sf1_times_sf2", op.rng_len(self.bJetsM), fourJet2BSel_sf_weighted, EqB(10, 0., 10.), title="Number of b jets (weighted_sf1_times_sf2)"))



        bTagEff_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "../", "btagEff_deepJet_2018.json")
        effparams = {
            "pt": lambda j: j.pt, "eta": lambda j: j.eta,
            "jetFlavor": lambda j: j.hadronFlavour, "workingPoint": "M"
        }
        btvEff = {"M": get_correction(bTagEff_file, "btagEff", params=effparams, sel=fourJet2BSel, defineOnFirstUse=True)}
        btvSF = lambda wp, flav: get_bTagSF_fixWP(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet", wp, flav, fourJet2BSel, syst_prefix="btagSF_", decorr_eras=False, era="2018_UL", decorr_wps=False)

#        plots.append(Plot.make1D(
 #           "leadBSF", b_jet1_sf, fourJet2BSel, EqB(100, 0., 1.), title="leading B-jet SF"))

        btvWeight = makeBtagWeightMeth1a(self.bJetsM, "btagDeepFlavB", ["M"], {"M": 0.2783}, btvSF, btvEff)

        #bTagWeight = op.rng_product(self.bJetsM, btvWeight)
        fourJet2BSel_true_weight = fourJet2BSel.refine("btag", weight=btvWeight)
        cfr.add(fourJet2BSel_true_weight, "weighted_true_makeBtagWeightMeth1a")
        plots.append(Plot.make1D(
            "nBJets_weighted", op.rng_len(self.bJetsM), fourJet2BSel_true_weight, EqB(10, 0., 10.), title="Number of b jets"))


        return plots

class test_btagWeight_maznli(NanoAODHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend", type=str, default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import NanoAODDescription
        return super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v5", year="2016", isMC=self.isMC(sample)), backend=self.args.backend)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        

        era = sampleCfg.get("era") if sampleCfg else None

        plots = []

        cfr = CutFlowReport("yields", recursive=True)
        plots.append(cfr)
        cfr.add(noSel, "Initial")

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight, autoSyst=False)
        self.rawJets = op.select(t.Jet, jetDef)

        fourJetSel = noSel.refine("fourJets", cut=[op.rng_len(self.rawJets) > 3])
        cfr.add(fourJetSel, "At least 4 jets")

        plots.append(Plot.make1D(
            "nJets", op.rng_len(self.rawJets), fourJetSel, EqB(10, 0., 10.), title="Number of raw jets (no-selection)"))

        self.bJetsM = bTagDef(self.rawJets, era, "M", "btagDeepFlavB")
        fourJettwoBSel = fourJetSel.refine("fourJets2B", cut=[op.rng_len(self.bJetsM) > 1])
        cfr.add(fourJettwoBSel, "At least two b-jets")

        plots.append(Plot.make1D(
            "nBJets", op.rng_len(self.bJetsM), fourJettwoBSel, EqB(10, 0., 10.), title="Number of b jets"))


        from bamboo.scalefactors import get_bTagSF_fixWP, get_correction, makeBtagWeightMeth1a

        params = {
         "pt": lambda j: j.pt, "abseta": lambda j: op.abs(j.eta), "working_point": "M", "flavor": 5}
        bjetSFs = get_correction(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet_comb", params=params,
                              systParam="systematic", systNomName="central", sel=fourJettwoBSel,
                             defineOnFirstUse=True)
        b_jet1_sf = bjetSFs(self.bJetsM[0])
        b_jet2_sf = bjetSFs(self.bJetsM[1])
        if self.isMC(sample): 
            fourJet2BSel_sf_weighted = fourJettwoBSel.refine("btag_sf", weight= [b_jet1_sf, b_jet2_sf])        
        cfr.add(fourJet2BSel_sf_weighted, "weighted_sf1_times_sf2")
        plots.append(Plot.make1D(
            "nBJets_sf1_times_sf2", op.rng_len(self.bJetsM), fourJet2BSel_sf_weighted, EqB(10, 0., 10.), title="Number of b jets (weighted_sf1_times_sf2)"))

        plots.append(Plot.make1D(
            "nBJets_doubleCheck", op.rng_len(self.bJetsM), fourJettwoBSel, EqB(10, 0., 10.), title="Number of b jets"))


        bTagEff_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), "../", "BTaggEffProduction_allSignalStats", "BTagEff_maps_UL18.json.gz")

        effparams = {
            "pt": lambda j: j.pt, "eta": lambda j: op.abs(j.eta)#,
            #"jetFlavor": lambda j: j.hadronFlavour, "workingPoint": "M"
        }
        btvEff = {"M": get_correction(bTagEff_file, "bflav_ratio", params=effparams, sel=fourJettwoBSel, systParam="ValType", systNomName="sf", defineOnFirstUse=True)}
        btvSF = lambda wp, flav: get_bTagSF_fixWP(localizePOGSF(era, "BTV", "btagging.json.gz"), "deepJet", wp, flav, fourJettwoBSel, syst_prefix="btagSF_", decorr_eras=False, era="2018_UL", decorr_wps=False)

#        plots.append(Plot.make1D(
 #           "leadBSF", b_jet1_sf, fourJet2BSel, EqB(100, 0., 1.), title="leading B-jet SF"))

        btvWeight = makeBtagWeightMeth1a(self.bJetsM, "btagDeepFlavB", ["M"], {"M": 0.2783}, btvSF, btvEff)

        #bTagWeight = op.rng_product(self.bJetsM, btvWeight)
        fourJet2BSel_true_weight = fourJettwoBSel.refine("btag", weight=btvWeight)
        cfr.add(fourJet2BSel_true_weight, "weighted_true_makeBtagWeightMeth1a")
        plots.append(Plot.make1D(
            "nBJets_weighted", op.rng_len(self.bJetsM), fourJet2BSel_true_weight, EqB(10, 0., 10.), title="Number of b jets"))


        return plots
