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
        self.out.branch("ToMaAn_weight", "F")
        self.out.branch("ToMaAn_MlMu", "F")
        self.out.branch("ToMaAn_MlMu_MCTruth", "F")
        self.out.branch("ToMaAn_MlMu_softMuMatched", "F")
        self.out.branch("ToMaAn_MlMu_noQCD", "F")
        self.out.branch("ToMaAn_MlMu_sameTopSelection", "F")
        self.out.branch("ToMaAn_MlMu_softMuPtCut15", "F")
        self.out.branch("ToMaAn_goodMuIndex", "i")
        self.out.branch("ToMaAn_goodSoftMuIndex", "i")
        self.out.branch("ToMaAn_sameTop", "F")
        self.out.branch("ToMaAn_dR_hardMu_softMu_sameTop", "F")
        self.out.branch("ToMaAn_dR_hardMu_softMu_oppositeTop", "F")

        outputFile.cd()
        self.histos={
        }
        self.histos['nEles'] = ROOT.TH1I( 'nEles', ';nr of es; Events',3,-0.5,2.5)
        self.histos['nMus'] = ROOT.TH1I( 'nMus', ';nr of #mus; Events',3,-0.5,2.5)
        self.histos['nTaus'] = ROOT.TH1I( 'nTaus', ';nr of #taus; Events',3,-0.5,2.5)
        self.histos['lepW_M'] = ROOT.TH1F( 'lepW_M', ';mass of leptonic W; Events',160,0,160)
        self.histos['hadW_M'] = ROOT.TH1F( 'hadW_M', ';mass of hadronic W; Events',160,0,160)
        self.histos['lepTop_M'] = ROOT.TH1F( 'lepTop_M', ';mass of leptonic Top; Events',150,100,250)
        self.histos['hadTop_M'] = ROOT.TH1F( 'hadTop_M', ';mass of hadronic Top; Events',150,100,250)
        self.histos['nSemiMuBDecays'] = ROOT.TH1I( 'nSemiMuBDecays', ';nr of semi-#mu B Decays; Events',3,-0.5,2.5)
        self.histos['nRecMus'] = ROOT.TH1I( 'nRecMus', ';nr of reconstructed #mus; Events',6,-0.5,5.5)
        self.histos['RecMu_pt'] = ROOT.TH1F( 'RecMu_pt', ';#mus p_T; Muons',200,0.,200.)
        self.histos['RecMu_eta'] = ROOT.TH1F( 'RecMu_eta', ';#mus #eta; Muons',60,-3.,3.)
        self.histos['RecMu_phi'] = ROOT.TH1F( 'RecMu_phi', ';#mus #phi; Muons',63,-3.15,3.15)
        self.histos['RecMu_isTight'] = ROOT.TH1I( 'RecMu_isTight', ';#mus tightId; Muons',2,0,2)
        self.histos['RecMu_isTracker_beforeTightId'] = ROOT.TH1I( 'RecMu_isTracker_beforeTightId', ';#mus isTracker; Muons',2,0,2)
        self.histos['RecMu_isGlobal_beforeTightId'] = ROOT.TH1I( 'RecMu_isGlobal_beforeTightId', ';#mus isGlobal; Muons',2,0,2)
        self.histos['RecMu_nTrackerLayers_beforeTightId'] = ROOT.TH1I( 'RecMu_nTrackerLayers_beforeTightId', ';#mus nTrackerLayers; Muons',20,0,20)
        self.histos['RecMu_nStations_beforeTightId'] = ROOT.TH1I( 'RecMu_nStations_beforeTightId', ';#mus nStations; Muons',6,0,6)
        self.histos['RecMu_segmentComp_beforeTightId'] = ROOT.TH1F( 'RecMu_segmentComp_beforeTightId', ';#mus segmentComp; Muons',100,0.,1.)
        self.histos['RecMu_dxy_beforeTightId'] = ROOT.TH1F( 'RecMu_dxy_beforeTightId', ';#mus ; Muons',30,0.,3.)
        self.histos['RecMu_dz_beforeTightId'] = ROOT.TH1F( 'RecMu_dz_beforeTightId', ';#mus ; Muons',100,0.,10.)
        self.histos['RecMu_isTracker'] = ROOT.TH1I( 'RecMu_isTracker', ';#mus isTracker; Muons',2,0,2)
        self.histos['RecMu_isGlobal'] = ROOT.TH1I( 'RecMu_isGlobal', ';#mus isGlobal; Muons',2,0,2)
        self.histos['RecMu_nTrackerLayers'] = ROOT.TH1I( 'RecMu_nTrackerLayers', ';#mus nTrackerLayers; Muons',20,0,20)
        self.histos['RecMu_nStations'] = ROOT.TH1I( 'RecMu_nStations', ';#mus nStations; Muons',6,0,6)
        self.histos['RecMu_segmentComp'] = ROOT.TH1F( 'RecMu_segmentComp', ';#mus segmentComp; Muons',100,0.,1.)
        self.histos['RecMu_relIso'] = ROOT.TH1F( 'RecMu_relIso', ';#mus relIso; Muons',40,0.,2.)
        self.histos['RecMu_dxy'] = ROOT.TH1F( 'RecMu_dxy', ';#mus ; Muons',30,0.,3.)
        self.histos['RecMu_dz'] = ROOT.TH1F( 'RecMu_dz', ';#mus ; Muons',100,0.,10.)
        self.histos['RecMu_jetIdx'] = ROOT.TH1I( 'RecMu_jetIdx', ';#mus ; Muons',22,-1,21)
        self.histos['nGoodRecMus'] = ROOT.TH1I( 'nGoodRecMus', ';nr of good reconstructed #mus; Events',6,-0.5,5.5)
        self.histos['DR_recMu_vs_genMu_inSemiMuTTbar'] = ROOT.TH1F( 'DR_recMu_vs_genMu_inSemiMuTTbar', ';#Delta R(#mu_gen,#mu_rec); Muons',16,0.,0.8)
        self.histos['DR_recMu_vs_associatedGenPar_inSemiMuTTbar'] = ROOT.TH1F( 'DR_recMu_vs_associatedGenPar_inSemiMuTTbar', ';#Delta R(associatedGenPar,#mu_rec); Muons',16,0.,0.8)
        self.histos['DR_recMu_vs_associatedGenPar_inOtherTTbar'] = ROOT.TH1F( 'DR_recMu_vs_associatedGenPar_inOtherTTbar', ';#Delta R(associatedGenPar,#mu_rec); Muons',16,0.,0.8)
        self.histos['RecMu_associatedGenPartFlav_inSemiMuTTbar'] = ROOT.TH1I( 'RecMu_associatedGenPartFlav_inSemiMuTTbar', '; associatedGenPar to recMu; Muons',16,0,16)
        self.histos['RecMu_associatedGenPartFlav_inOtherTTbar'] = ROOT.TH1I( 'RecMu_associatedGenPartFlav_inOtherTTbar', '; associatedGenPar to recMu; Muons',16,0,16)

        self.histos['passedEvts'] = ROOT.TH1I( 'passedEvts', ';nr of events passed; Events',6,0,6)
        
        self.nrSemiMuBDecays = 0
        self.nrSoftRecMu = 0
        self.nrSoftRecMu_genMatched = 0
        self.nrSoftRecMu_sameTopSelection = 0
        self.nrSoftRecMu_sameTopSelection_genMatched = 0


        self.histos['DR_goodMu_jets'] = ROOT.TH1F( 'DR_goodMu_jets', ';#Delta R(good #mu,jets); Jets',30,0.,3.)

        self.histos[ 'jets_nConstituents' ] = ROOT.TH1I( 'jets_nConstituents', ';nConstituents; Jets',50,0,50) 
        self.histos[ 'jets_muEF' ] = ROOT.TH1F( 'jets_muEF', ';muEF; Jets',100,0.,1.)
        self.histos[ 'jets_chEmEF' ] = ROOT.TH1F( 'jets_chEmEF', ';chEmEF; Jets',100,0.,1.)
        self.histos[ 'jets_neEmEF' ] = ROOT.TH1F( 'jets_neEmEF', ';neEmEF; Jets',100,0.,1.)
        self.histos[ 'jets_chHEF' ] = ROOT.TH1F( 'jets_chHEF', ';chHEF; Jets',100,0.,1.)
        self.histos[ 'jets_neHEF' ] = ROOT.TH1F( 'jets_neHEF', ';neHEF; Jets',100,0.,1.)
        self.histos[ 'jets_jetId' ] = ROOT.TH1I( 'jets_jetId', ';jetId; Jets',10,0,10)
        self.histos['nGoodJets'] = ROOT.TH1I( 'nGoodJets', ';nr of good reconstructed jets; Events',11,-0.5,10.5)
        self.histos['nrBTaggedJets'] = ROOT.TH1I( 'nrBTaggedJets', ';nr of BTagged jets; Events',5,-0.5,4.5)

        self.histos['softMu_looseId'] = ROOT.TH1I( 'softMu_looseId', ';soft #mus looseId; Muons',2,0,2)
        self.histos['softMu_mediumId'] = ROOT.TH1I( 'softMu_mediumId', ';soft #mus mediumId; Muons',2,0,2)
        self.histos['softMu_softMvaId'] = ROOT.TH1I( 'softMu_softMvaId', ';soft #mus softMvaId; Muons',2,0,2)
        self.histos['softMvaId_vs_mediumId'] = ROOT.TH2I( 'softMvaId_vs_mediumId', ';soft #mus softMvaId; soft #mus mediumId',2,0,2,2,0,2)
        self.histos['looseId_vs_mediumId'] = ROOT.TH2I( 'looseId_vs_mediumId', ';soft #mus looseId; soft #mus mediumId',2,0,2,2,0,2)
        self.histos['nSoftMus'] = ROOT.TH1I( 'nSoftMus', ';nr of soft muons; Events',11,-0.5,10.5)
        self.histos['softRecMu_pt'] = ROOT.TH1F( 'softRecMu_pt', ';soft #mus p_T; Soft Muons',100,0.,100.)
        self.histos['DR_softRecMu_softGenMu'] = ROOT.TH1F( 'DR_softRecMu_softGenMu', ';#Delta R(soft rec #mu,soft gen #mu); events',60,0.,6.)
        self.histos['nGenSoftMus'] = ROOT.TH1I( 'nGenSoftMus', ';nr of soft muons at gen-level; Events',11,-0.5,10.5)
        self.histos['softMu_associatedGenPartFlav'] = ROOT.TH1I( 'softMu_associatedGenPartFlav', ';the flavour that soft muon is associated to; Muons',16,0,16)

        self.histos['dPt_softVsMatchedGenMu_OverPt'] = ROOT.TH1F( 'dPt_softVsMatchedGenMu_OverPt', ';(p_T^gen-p_T^matched)/p_T^gen; Matched Soft Muons',100,-1.,1.)
        self.histos[ 'MlMu' ] = ROOT.TH1F( 'MlMu', ';M_{#mu_{soft},#mu_{hard}}; Events',250,0.,250.)
        self.histos[ 'MlMu_MCTruth' ] = ROOT.TH1F( 'MlMu_MCTruth', ';MCTruth M_{#mu_{soft},#mu_{hard}}; Events',250,0.,250.)


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.histos[ 'passedEvts' ].Fill(0,self.nrSemiMuBDecays)
        self.histos[ 'passedEvts' ].Fill(1,self.nrSoftRecMu)
        self.histos[ 'passedEvts' ].Fill(2,self.nrSoftRecMu_genMatched)
        self.histos[ 'passedEvts' ].Fill(3,self.nrSoftRecMu_sameTopSelection)
        self.histos[ 'passedEvts' ].Fill(4,self.nrSoftRecMu_sameTopSelection_genMatched)

        #write histos to output file
        outputFile.cd()
        print' -->> total nr of events with SemiMuBDecays are : ', self.nrSemiMuBDecays
        print' -->> nr of events with SofRecMu is found to be : ', self.nrSoftRecMu
        print' -->> nr of events with gen-matched SofRecMu is : ', self.nrSoftRecMu_genMatched
        print' -->> nr of events with same-top selection cut (dr(softMu,hardMu)<2.) : ', self.nrSoftRecMu_sameTopSelection
        print' -->> nr of events with same-top selection cut plus soft muon matched : ', self.nrSoftRecMu_sameTopSelection_genMatched

        for _,h in self.histos.items():
            h.SetDirectory(outputFile)
            h.Sumw2()
            h.Write()

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
        self.out.fillBranch( 'ToMaAn_weight' , 1. )
        self.out.fillBranch( 'ToMaAn_MlMu' , -100. )
        self.out.fillBranch( 'ToMaAn_MlMu_MCTruth' , -100. )
        self.out.fillBranch( 'ToMaAn_MlMu_softMuMatched' , -100. )
        self.out.fillBranch( 'ToMaAn_MlMu_noQCD' , -100. )
        self.out.fillBranch( 'ToMaAn_MlMu_sameTopSelection' , -100. )
        self.out.fillBranch( 'ToMaAn_MlMu_softMuPtCut15' , -100. )
        self.out.fillBranch( 'ToMaAn_goodMuIndex' , 10 )
        self.out.fillBranch( 'ToMaAn_goodSoftMuIndex' , 10 )
        self.out.fillBranch( 'ToMaAn_sameTop' , 0. )
        self.out.fillBranch( 'ToMaAn_dR_hardMu_softMu_sameTop', -1.)
        self.out.fillBranch( 'ToMaAn_dR_hardMu_softMu_oppositeTop', -1.)

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

        genMus = []

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
                genMus.append(part)
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
        genSoftMus = []

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
                       genSoftMus.append(p)
                   else:
                       nrSemiMuBs+=1
                       self.out.fillBranch("ToMaAn_GenSoftMu_pt", p.pt)
                       genSoftMus.append(p)

        if len(genSoftMus) == 2 : 
            #print("======== before sort  ",genSoftMus[0].pt,genSoftMus[1].pt)
            genSoftMus.sort(key=lambda x: x.pt,reverse=True)
            #print("========  ",genSoftMus[0].pt,genSoftMus[1].pt)

        if len(genSoftMus) > 0 and len(genMus) == 1 :
            mcTruthMass = (genSoftMus[0].p4()+genMus[0].p4()).M()
            self.out.fillBranch("ToMaAn_MlMu_MCTruth", mcTruthMass)
            self.histos[ 'MlMu_MCTruth' ].Fill( mcTruthMass )
            

        self.histos[ 'nSemiMuBDecays' ].Fill( nrSemiMuBs  )
        # move the following two lines to the end of code, to make a fair comparison 
        # between nr of soft reconstructed muons with the nr of muons found in B decays
        #if nrSemiMuBs>0:
        #    self.nrSemiMuBDecays+=1

        #--------------------------------------------
        #--- RecParticles Analysis ------------------
        #--------------------------------------------

        muons = Collection(event, "Muon")
        good_muons = []
        good_muons_indices = []
        m_counter = -1
        for m in muons:
            m_counter += 1
            self.histos[ 'RecMu_pt' ].Fill( m.pt )
            self.histos[ 'RecMu_eta' ].Fill( m.eta )
            self.histos[ 'RecMu_phi' ].Fill( m.phi )


            # pt cut
            if m.pt < 27: continue

            # eta cut
            if abs(m.eta) > 2.5: continue

            self.histos[ 'RecMu_isTight' ].Fill( m.tightId )
            self.histos[ 'RecMu_isTracker_beforeTightId' ].Fill( m.isTracker )
            self.histos[ 'RecMu_isGlobal_beforeTightId' ].Fill( m.isGlobal )
            self.histos[ 'RecMu_nTrackerLayers_beforeTightId' ].Fill( m.nTrackerLayers )
            self.histos[ 'RecMu_nStations_beforeTightId' ].Fill( m.nStations )
            self.histos[ 'RecMu_segmentComp_beforeTightId' ].Fill( m.segmentComp )
            self.histos[ 'RecMu_dxy_beforeTightId' ].Fill( m.dxy )
            self.histos[ 'RecMu_dz_beforeTightId' ].Fill( m.dz )
            
            # tight id 
            if not m.tightId: continue

            self.histos[ 'RecMu_isTracker' ].Fill( m.isTracker )
            self.histos[ 'RecMu_isGlobal' ].Fill( m.isGlobal )
            self.histos[ 'RecMu_nTrackerLayers' ].Fill( m.nTrackerLayers )
            self.histos[ 'RecMu_nStations' ].Fill( m.nStations )
            self.histos[ 'RecMu_segmentComp' ].Fill( m.segmentComp )
            self.histos[ 'RecMu_dxy' ].Fill( m.dxy )
            self.histos[ 'RecMu_dz' ].Fill( m.dz )

            # the explicit cut values can be found at https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
            # with the tight id requirement, the below values are implicitly applied. 
            # isGlobal, nTrackerLayers>=6, nStations>=2, dxy<0.2, dz<0.5, numberOfValidPixelHits() > 0, numberOfValidMuonHits() > 0
            # chi2/ndof of the global-muon track fit < 10 ==> all of the quality cuts mentioned in AN2015_218 are applied

            self.histos[ 'RecMu_relIso' ].Fill( m.pfRelIso04_all )

            # tight iso
            # according to the following page:
            # https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection#Particle_Flow_isolation
            # the tight working point is equivalent to cut value < 0.15
            if not m.pfRelIso04_all < 0.15: continue


            self.histos[ 'RecMu_jetIdx' ].Fill( m.jetIdx )

            self.out.fillBranch("ToMaAn_RecMu_pt", m.pt)

            # now save good muons
            good_muons.append(m)
            good_muons_indices.append(m_counter)


        self.histos[ 'nGoodRecMus' ].Fill( len(good_muons) )
        self.histos[ 'nRecMus' ].Fill( event.nMuon  )

        hard_mu_mother = 0.
        if len(good_muons) == 1: 
            #to check the genPartFlav method: in the semi-muonic ttbar events, the selected good muon should be 
            # matched to that generated muon
            if good_muons[0].genPartIdx >= 0 :
             associatedGenPar = particles[good_muons[0].genPartIdx] 
             mothidx  = associatedGenPar.genPartIdxMother
             if mothidx >= 0 :
              moth    = particles[mothidx]
              mothpid = moth.pdgId
              #print(" ********* mothpid for associatedGenPar of good muon : ",mothpid)
              while abs(mothpid) != 6 and mothidx > 0 :
               mothidx  = moth.genPartIdxMother
               if mothidx >= 0 :
                moth    = particles[mothidx]
                mothpid = moth.pdgId
              #print(" ********* mothpid for associatedGenPar of good muon : ",mothpid)
              hard_mu_mother = mothpid
              
             if len(genMus) == 1:
                self.histos[ 'DR_recMu_vs_genMu_inSemiMuTTbar' ].Fill( genMus[0].DeltaR( good_muons[0] )  )
                self.histos[ 'RecMu_associatedGenPartFlav_inSemiMuTTbar' ].Fill( good_muons[0].genPartFlav )
                self.histos[ 'DR_recMu_vs_associatedGenPar_inSemiMuTTbar' ].Fill( associatedGenPar.DeltaR( good_muons[0] )  )
             if not len(genMus) == 1:  
                self.histos[ 'RecMu_associatedGenPartFlav_inOtherTTbar' ].Fill( good_muons[0].genPartFlav )
                self.histos[ 'DR_recMu_vs_associatedGenPar_inOtherTTbar' ].Fill( associatedGenPar.DeltaR( good_muons[0] )  )

         
        if not len(good_muons) == 1: 
            return False 
        #print(" ********* good muon index: ",good_muons_indices[0])
        self.out.fillBranch("ToMaAn_goodMuIndex", good_muons_indices[0])

        #--------------------------------------------
        #--- Jet Selection --------------------------
        #--------------------------------------------

        jets = Collection(event, "Jet")
        good_jets = []
        good_bjets_indices = []
        j_counter = -1
        for j in jets:
            j_counter += 1
            if not j.pt > 30: continue
            if not abs(j.eta) < 2.4: continue
            self.histos[ 'jets_jetId' ].Fill( j.jetId )
            # from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets, jetId==6 means 
            # pass tight and tightLepVeto ID
            # , also the values of the below identification parameters are cut over according to:
            # https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
            if j.jetId == 6:

             self.histos[ 'DR_goodMu_jets' ].Fill( j.DeltaR(good_muons[0]) )
             if not j.DeltaR(good_muons[0]) > 0.4: continue

             self.histos[ 'jets_nConstituents' ].Fill( j.nConstituents )
             self.histos[ 'jets_muEF' ].Fill( j.muEF )
             self.histos[ 'jets_chEmEF' ].Fill( j.chEmEF )
             self.histos[ 'jets_neEmEF' ].Fill( j.neEmEF )
             self.histos[ 'jets_chHEF' ].Fill( j.chHEF )
             self.histos[ 'jets_neHEF' ].Fill( j.neHEF )

             good_jets.append(j)

             # asking for b-tag identification
             # starting from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation,  the exact wps for Summer20UL16 were 
             # not available, for the moment the below cut values corresponding to a tight (misId rate ~ 0.1%) selection by deepJet 
             # b-tag algo are applied:
             # btagDeepFlavB pre-VFP 2016:  Loose, Medium, Tight: 0.0508, 0.2598, 0.6502
             # btagDeepFlavB post-VFP 2016: Loose, Medium, Tight: 0.0480, 0.2489, 0.6377

             #if not j.btagDeepFlavB > 0.25 : continue # medium criteria
             if not j.btagDeepFlavB > 0.64 : continue # tight criteria

             good_bjets_indices.append(j_counter)

        self.histos[ 'nGoodJets' ].Fill( len(good_jets) )
        self.histos[ 'nrBTaggedJets' ].Fill( len(good_bjets_indices) )

        if not len(good_jets) > 3: 
            return False 

        # put the BTagged jets requirement with cautious ! 
        if not len(good_bjets_indices) > 0: 
            return False 

        #--------------------------------------------
        #--- Soft Muon Selection --------------------
        #--------------------------------------------

        good_soft_muons = []
        good_soft_muons_indices = []
        soft_m_index = -1
        for m in muons:
            soft_m_index += 1
            # pt cut
            if m.pt < 10.: continue
            # eta cut
            if abs(m.eta) > 2.4: continue
        
            self.histos[ 'softMu_looseId' ].Fill( m.looseId )
            self.histos[ 'softMu_mediumId' ].Fill( m.mediumId )
            self.histos[ 'softMu_softMvaId' ].Fill( m.softMvaId )

            self.histos[ 'looseId_vs_mediumId' ].Fill( m.looseId , m.mediumId)
            self.histos[ 'softMvaId_vs_mediumId' ].Fill( m.softMvaId , m.mediumId)

            # id cut
            if not m.mediumId: continue


            #print("")
            #print("event: ",event.event)
            if m.jetIdx in good_bjets_indices and soft_m_index not in good_muons_indices:
             #for bj in range(len(good_bjets_indices)):
                #print("--------- b-jet with index ",good_bjets_indices[bj]," among jet collection")
                #print("          with dr : ", m.DeltaR(jets[good_bjets_indices[bj]]))
                #print("     ---- soft mu jetIdx: ",m.jetIdx)
                #print("************************  soft reconstructed muon found !!")
                #print("         whether identified as softMvaId? ",m.softMvaId)
                #print("         , and pt: ",m.pt,"  , and eta: ",m.eta)

                good_soft_muons.append(m)
                good_soft_muons_indices.append(soft_m_index)
                self.histos[ 'softMu_associatedGenPartFlav' ].Fill( m.genPartFlav )

        self.histos[ 'nGenSoftMus' ].Fill( nrSemiMuBs  ) ## there was a histogram named "nSemiMuBDecays" , filled above, before applying selection cuts. 
        self.histos[ 'nSoftMus' ].Fill( len(good_soft_muons) )

        if nrSemiMuBs > 0:
            self.nrSemiMuBDecays += 1
        soft_mu_mother = 0.
        ObservableMass = -1.
        if len(good_soft_muons) > 0:
            self.nrSoftRecMu += 1
            observable_4vec = good_soft_muons[0].p4() + good_muons[0].p4()
            self.histos[ 'MlMu' ].Fill( observable_4vec.M() )
            self.out.fillBranch("ToMaAn_MlMu", observable_4vec.M())
            ObservableMass = observable_4vec.M()
            self.out.fillBranch("ToMaAn_goodSoftMuIndex", good_soft_muons_indices[0])
            #print("   ********* good soft muon index: ",good_soft_muons_indices[0])
            if good_soft_muons[0].genPartIdx >= 0 :
             associatedGenPar = particles[good_soft_muons[0].genPartIdx] 
             mothidx  = associatedGenPar.genPartIdxMother
             if mothidx >= 0 :
              moth    = particles[mothidx]
              mothpid = moth.pdgId
              #print("   before loop ********* mothpid for associatedGenPar of good soft muon : ",mothpid)
              while abs(mothpid) != 6 and mothidx > 0:
               #print("   ********* mothpid for associatedGenPar of good soft muon : ",mothpid, "  and mother index: ", mothidx)
               mothidx  = moth.genPartIdxMother
               if mothidx >= 0 :
                moth    = particles[mothidx]
                mothpid = moth.pdgId
              #print("   ********* mothpid for associatedGenPar of good soft muon : ",mothpid)
              soft_mu_mother = mothpid
              #if hard_mu_mother*mothpid > 0 : 
               #self.out.fillBranch("ToMaAn_sameTop", +1.)


        if len(good_soft_muons) == 0:
            return False
              
        self.histos[ 'softRecMu_pt' ].Fill( good_soft_muons[0].pt )

        if len(genSoftMus) > 0:
            self.histos[ 'DR_softRecMu_softGenMu' ].Fill( good_soft_muons[0].DeltaR(genSoftMus[0]) )
            if good_soft_muons[0].DeltaR(genSoftMus[0]) < 0.3:
             self.nrSoftRecMu_genMatched += 1
             self.out.fillBranch("ToMaAn_MlMu_softMuMatched", ObservableMass)
             self.histos[ 'dPt_softVsMatchedGenMu_OverPt' ].Fill( (genSoftMus[0].pt-good_soft_muons[0].pt)/genSoftMus[0].pt )

        if hard_mu_mother*soft_mu_mother == 36 : 
         self.out.fillBranch("ToMaAn_sameTop", +1.)
         self.out.fillBranch("ToMaAn_dR_hardMu_softMu_sameTop", good_soft_muons[0].DeltaR(good_muons[0]))
        elif hard_mu_mother*soft_mu_mother == -36 : 
         self.out.fillBranch("ToMaAn_sameTop", -1.)
         self.out.fillBranch("ToMaAn_dR_hardMu_softMu_oppositeTop", good_soft_muons[0].DeltaR(good_muons[0]))

        if ObservableMass > 15. : 
            self.out.fillBranch("ToMaAn_MlMu_noQCD", ObservableMass)
        if good_soft_muons[0].DeltaR(good_muons[0]) < 2. : 
            self.out.fillBranch("ToMaAn_MlMu_sameTopSelection", ObservableMass)
            self.nrSoftRecMu_sameTopSelection += 1
            if len(genSoftMus) > 0 and good_soft_muons[0].DeltaR(genSoftMus[0]) < 0.3:
             self.nrSoftRecMu_sameTopSelection_genMatched += 1
        if good_soft_muons[0].pt > 15. : 
            self.out.fillBranch("ToMaAn_MlMu_softMuPtCut15", ObservableMass)

        return True
