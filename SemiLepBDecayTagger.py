import os
import copy
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class SemiLepBDecayTagger(Module):
    def __init__(self):
        self.histos={}

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ToMaAn_GenSoftMu_pt", "F")
        self.out.branch("ToMaAn_RecMu_pt", "F")

        outputFile.cd()
        self.histos={
        }
        self.histos['nEles'] = ROOT.TH1I( 'nEles', ';nr of es; Events',3,-0.5,2.5)
        self.histos['nMus'] = ROOT.TH1I( 'nMus', ';nr of #mu s; Events',3,-0.5,2.5)
        self.histos['nTaus'] = ROOT.TH1I( 'nTaus', ';nr of #tau s; Events',3,-0.5,2.5)
        self.histos['lepW_M'] = ROOT.TH1F( 'lepW_M', ';mass of leptonic W; Events',160,0,160)
        self.histos['hadW_M'] = ROOT.TH1F( 'hadW_M', ';mass of hadronic W; Events',160,0,160)
        self.histos['lepTop_M'] = ROOT.TH1F( 'lepTop_M', ';mass of leptonic Top; Events',150,100,250)
        self.histos['hadTop_M'] = ROOT.TH1F( 'hadTop_M', ';mass of hadronic Top; Events',150,100,250)
        self.histos['nSemiMuBDecays'] = ROOT.TH1I( 'nSemiMuBDecays', ';nr of semi-#mu B Decays; Events',3,-0.5,2.5)
        self.histos['nRecMus'] = ROOT.TH1I( 'nRecMus', ';nr of reconstructed #mu s; Events',6,-0.5,5.5)
        self.histos['RecMu_pt'] = ROOT.TH1F( 'RecMu_pt', ';#mu s p_T^; Events',200,0.,200.)
        self.nrSemiMuBDecays = 0


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #write histos to output file
        outputFile.cd()
        for _,h in self.histos.items():
            h.SetDirectory(outputFile)
            h.Sumw2()
            h.Write()
        print' -->> total nr of events with SemiMuBDecays are : ', self.nrSemiMuBDecays

    def isBHadron(self, pId):
        if (((pId / 100) % 10 == 5) or (pId >= 5000 and pId <= 5999)):
           return True


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        nrEles = 0
        nrMus = 0
        nrTaus = 0
        nrSemiMuBs = 0
        self.out.fillBranch( 'ToMaAn_GenSoftMu_pt' , -100. )
        self.out.fillBranch( 'ToMaAn_RecMu_pt' , -100. )

        #--------------------------------------------
        #--- GenParticles Analysis ------------------
        #--------------------------------------------


        if not hasattr( event , 'nLHEPart' ):
            return False

        lheparts = None
        try :
            lheparts=Collection(event, 'LHEPart')
        except:
            return False

	p4_lep = ROOT.TLorentzVector()
	p4_nu = ROOT.TLorentzVector()
	p4_j1 = ROOT.TLorentzVector()
	p4_j2 = ROOT.TLorentzVector()
	p4_bs = []
        pdgIds_bs = []

        for part in lheparts:
            pdg = abs(part.pdgId)

            if pdg == 11 or pdg == 13 or pdg == 15:
                       p4_lep.SetPtEtaPhiM(part.pt,part.eta,part.phi,part.mass)
                       lepPdgId = part.pdgId
            if pdg == 12 or pdg == 14 or pdg == 16:
                       p4_nu.SetPtEtaPhiM(part.pt,part.eta,part.phi,part.mass)
            if pdg == 1 or pdg == 3:
                       p4_j1.SetPtEtaPhiM(part.pt,part.eta,part.phi,part.mass)
            if pdg == 2 or pdg == 4:
                       p4_j2.SetPtEtaPhiM(part.pt,part.eta,part.phi,part.mass)
            if pdg == 5:
                       p4_bs.append(part.p4())
                       pdgIds_bs.append(part.pdgId)

            if pdg == 11:
                nrEles += 1
            if pdg == 13:
                nrMus += 1
            if pdg == 15:
                nrTaus += 1


        self.histos[ 'nEles' ].Fill( nrEles  )
        self.histos[ 'nMus' ].Fill( nrMus  )
        self.histos[ 'nTaus' ].Fill( nrTaus  )

        #--------------------------------------------
        #--- W/Top Mass using GenParticles ----------
        #--------------------------------------------

        p4_hadW = p4_j1+p4_j2
        p4_lepW = p4_lep+p4_nu

        for i in range(len(p4_bs)):
            if lepPdgId*pdgIds_bs[i]<0:
                p4_lepTop = p4_lepW + p4_bs[i]
            else:
                p4_hadTop = p4_hadW + p4_bs[i]

        self.histos[ 'lepW_M' ].Fill( p4_lepW.M()  )
        self.histos[ 'hadW_M' ].Fill( p4_hadW.M()  )
        self.histos[ 'lepTop_M' ].Fill( p4_lepTop.M()  )
        self.histos[ 'hadTop_M' ].Fill( p4_hadTop.M()  )
        
        particles = Collection(event,"GenPart")

        for p in particles:
            mothidx  = p.genPartIdxMother
            if mothidx < 0 : continue
            moth    = particles[mothidx]
            mothpid = moth.pdgId
            if abs(p.pdgId) in [13]: # electrons removed from the list because we are only interested in semi-muonic B decays
              bHadronProduct = self.isBHadron(abs(mothpid))
              if bHadronProduct:
                  gmothidx  = moth.genPartIdxMother
                  if gmothidx < 0 : continue
                  gmoth    = particles[gmothidx]
                  gmothpid = gmoth.pdgId
                  if abs(gmothpid) not in [1,2,3,4,21]:
                   ggmothidx  = gmoth.genPartIdxMother
                   if ggmothidx < 0 : continue
                   ggmoth    = particles[ggmothidx]
                   ggmothpid = ggmoth.pdgId
                   if abs(gmothpid) not in [5]:
                    if abs(ggmothpid) not in [1,2,3,4,21]:
                       nrSemiMuBs+=1
                       self.out.fillBranch("ToMaAn_GenSoftMu_pt", p.pt)
                   else:
                       nrSemiMuBs+=1
                       self.out.fillBranch("ToMaAn_GenSoftMu_pt", p.pt)

        self.histos[ 'nSemiMuBDecays' ].Fill( nrSemiMuBs  )
        if nrSemiMuBs>0:
            self.nrSemiMuBDecays+=1

        #--------------------------------------------
        #--- RecParticles Analysis ------------------
        #--------------------------------------------

        muons = Collection(event, "Muon")
        for m in muons:
            self.histos[ 'RecMu_pt' ].Fill( m.pt  )
            self.out.fillBranch("ToMaAn_RecMu_pt", m.pt)

        self.histos[ 'nRecMus' ].Fill( event.nMuon  )

        return True


