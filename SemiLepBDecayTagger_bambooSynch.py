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
        self.out.branch("ToMaAn_gMothPid", "F")

        outputFile.cd()
        self.histos={
        }
        
        self.histos['nGenSoftMus'] = ROOT.TH1I( 'nGenSoftMus', ';nr of soft muons at gen-level; Events',11,-0.5,10.5)
        self.histos['grandmoth_pid_genSoftMus'] = ROOT.TH1I( 'grandmoth_pid_genSoftMus', ';grandmoth_pid_genSoftMus',25,0,25)
        self.histos['gMothPid'] = ROOT.TH1I( 'gMothPid', ';gMothPid',25,0,25)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #write histos to output file
        outputFile.cd()

        for _,h in self.histos.items():
            h.SetDirectory(outputFile)
            h.Sumw2()
            h.Write()

    def isBHadron(self, pId):
        if (((pId / 100) % 10 == 5) or (pId >= 5000 and pId <= 5999)):
           return True


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        nrSemiMuBs = 0
        self.out.fillBranch( 'ToMaAn_GenSoftMu_pt' , -100. )
        self.out.fillBranch( 'ToMaAn_gMothPid', -1.)

        #--------------------------------------------
        #--- GenParticles Analysis ------------------
        #--------------------------------------------
        particles = Collection(event,"GenPart")
        genSoftMus = []


        #print('  ')
        #print('  __________________________________________________________________________________________________________________________')
        #print('  ')
        #print('Analyzing event #{0} ... , with genTtbarId: {1}'.format(event.event, event.genTtbarId))
        #print('  ')

        for p in particles:
            mothidx  = p.genPartIdxMother
            if mothidx < 0 : continue
            moth    = particles[mothidx]
            mothpid = moth.pdgId
            if abs(p.pdgId) in [13] and ((p.statusFlags & (1 << 6))>0) and not ((p.statusFlags & (1 << 3))>0) and not ((p.statusFlags & (1 << 2))>0) and p.status==1:
              bHadronProduct = self.isBHadron(abs(mothpid))
              if bHadronProduct:
                  safe = False
                  gmothidx  = moth.genPartIdxMother
                  while gmothidx > 0:
                    gmoth    = particles[gmothidx]
                    gmothidx  = gmoth.genPartIdxMother
                    gmothpid = gmoth.pdgId
                    if abs(gmothpid) == 6: 
                        safe = True
                        break
                  if safe != True : continue 

                  if ((abs(particles[particles[mothidx].genPartIdxMother].pdgId) in [5]) or self.isBHadron(abs(particles[particles[mothidx].genPartIdxMother].pdgId))): 
                      if ((abs(particles[particles[particles[mothidx].genPartIdxMother].genPartIdxMother].pdgId) in [1,2,3,4]) or (abs(particles[particles[particles[particles[mothidx].genPartIdxMother].genPartIdxMother].genPartIdxMother].pdgId) in [1,2,3,4])):
                          continue;
                      nrSemiMuBs+=1
                      self.histos[ 'gMothPid' ].Fill( abs(gmothpid)  )
                      self.out.fillBranch("ToMaAn_gMothPid", abs(gmothpid) )
                      self.out.fillBranch("ToMaAn_GenSoftMu_pt", p.pt)
                      genSoftMus.append(p)
            
        if (nrSemiMuBs>=3): 
            print("\n\n************* >=3 nrSemiMuBs found! , nrSemiMuBs: ", nrSemiMuBs, ' and evt: ', event.event)
            for p in genSoftMus:
                print('  ')
                print('+++++++++++++++++++++++++++++ A muon from BHadron decay is found   ')
                mup4 = ROOT.TLorentzVector()
                mup4.SetPtEtaPhiM(p.pt,p.eta,p.phi,p.mass)
                print('      +++ pt: {0}, eta: {1}, phi: {2}, and energy: {3}'.format(p.pt,p.eta,p.phi,mup4.E()))
                mothidx  = p.genPartIdxMother
                print('      +++ mother pdgId: {0}, g-mother pdgId: {1}, gg-mother pdgId: {2}, ggg-mother pdgId: {3}'.format(particles[mothidx].pdgId,particles[particles[mothidx].genPartIdxMother].pdgId,particles[particles[particles[mothidx].genPartIdxMother].genPartIdxMother].pdgId,particles[particles[particles[particles[mothidx].genPartIdxMother].genPartIdxMother].genPartIdxMother].pdgId))
                print('  ')
            print('  __________________________________________________________________________________________________________________________')
        self.histos[ 'nGenSoftMus' ].Fill( nrSemiMuBs  )
        if (nrSemiMuBs>0): self.histos[ 'grandmoth_pid_genSoftMus' ].Fill( abs(particles[particles[genSoftMus[0].genPartIdxMother].genPartIdxMother].pdgId)  )

        return True
