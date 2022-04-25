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
        self.out.branch("SemiMuBGenParticles_pt", "F")

        outputFile.cd()
        self.histos={
        }
        self.histos['nEles'] = ROOT.TH1I( 'nEles', ';nr of es; Events',3,-0.5,2.5)
        self.histos['nMus'] = ROOT.TH1I( 'nMus', ';nr of #mu s; Events',3,-0.5,2.5)
        self.histos['nTaus'] = ROOT.TH1I( 'nTaus', ';nr of #tau s; Events',3,-0.5,2.5)
        self.histos['nSemiMuBDecays'] = ROOT.TH1I( 'nSemiMuBDecays', ';nr of semi-#mu B Decays; Events',3,-0.5,2.5)
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

        if not hasattr( event , 'nLHEPart' ):
            return False

        lheparts = None
        try :
            lheparts=Collection(event, 'LHEPart')
        except:
            return False

        for part in lheparts:
            pdg = abs(part.pdgId)
            if pdg == 11:
                nrEles += 1
            if pdg == 13:
                nrMus += 1
            if pdg == 15:
                nrTaus += 1


        self.histos[ 'nEles' ].Fill( nrEles  )
        self.histos[ 'nMus' ].Fill( nrMus  )
        self.histos[ 'nTaus' ].Fill( nrTaus  )

        
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
                  if gmothpid not in [1,2,3,4,21]:
                   print '<<< semi-mu B-hadron decay is found >>> '
                   nrSemiMuBs+=1
                   self.nrSemiMuBDecays+=1
                   self.out.fillBranch("SemiMuBGenParticles_pt", moth.pt)

        self.histos[ 'nSemiMuBDecays' ].Fill( nrSemiMuBs  )

        return True


