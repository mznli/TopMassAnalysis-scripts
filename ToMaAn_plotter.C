////////////////////////////////////////////////////////////////////////////////////////////////////////
/// script to plot the histograms from TopMassAnalysis project /////////////////////////////////////////
/// To run the script, complie and run in a root session, e.g. /////////////////////////////////////////
/// root [0] .L ToMaAn_plotter.C++ /////////////////////////////////////////////////////////////////////
/// root [0] ToMaAn_plotter()      /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TText.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <utility> 
#include <TLorentzVector.h>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLine.h"
#include "TStyle.h"

using namespace std;

void AddOverAndUnderFlow(TH1 * Histo) {
  Histo->SetBinContent(1,
                       Histo->GetBinContent(0) + Histo->GetBinContent(1));
  Histo->SetEntries(Histo->GetEntries()-1);
  Histo->SetBinError(1,
                     sqrt(Histo->GetBinError(0)*Histo->GetBinError(0)+
                          Histo->GetBinError(1)*Histo->GetBinError(1) ));
  Histo->SetBinContent(0, 0.0);
  Histo->SetEntries(Histo->GetEntries()-1);
  Histo->SetBinContent(Histo->GetNbinsX(),
                       Histo->GetBinContent(Histo->GetNbinsX()  )+
                       Histo->GetBinContent(Histo->GetNbinsX()+1) );
  Histo->SetEntries(Histo->GetEntries()-1);
  Histo->SetBinError(Histo->GetNbinsX(),
                     sqrt(Histo->GetBinError(Histo->GetNbinsX()  )*
                          Histo->GetBinError(Histo->GetNbinsX()  )+
                          Histo->GetBinError(Histo->GetNbinsX()+1)*
                          Histo->GetBinError(Histo->GetNbinsX()+1)  ));
  Histo->SetBinContent(Histo->GetNbinsX() + 1, 0.0);
  Histo->SetEntries(Histo->GetEntries()-1);
}

	
class MakePlot {

private: 

    TString m_name;		
    TString m_title;
    int m_nBins;
    string m_path;
    bool inCutFlow;

    int sampleColors[4]  = {kRed, kGreen-2, kAzure, kOrange-3};

    TH1D * h_samples[4]; 
    TH1D * h_MCTruth = NULL; 
    THStack * myStack;

public:
    
    MakePlot(TString name, TString title, string pathPNG, bool toBeUsed = false) {
	m_name = name;
	m_title = title;
	m_path = pathPNG;
	inCutFlow = toBeUsed;

	for (int s = 0; s < 4; s++) {
		h_samples[s] = new TH1D();
		h_samples[s]->Sumw2();
		h_samples[s]->SetStats();
	    }

	myStack = new THStack("hs","");
    }

    bool useInMCTable() {
      return inCutFlow;
    }	    
    
    TString getName() {
      return m_name;
    }	    
    
    void setNBins(TH1 * h) {
      m_nBins = h->GetNbinsX();
    }	    
    
    void fillHisto(TH1 * h, TString fIn) {

	//for (int s = 0; s < 4; s++) AddOverAndUnderFlow(h_samples[s]);

	if (fIn.Contains("SemiLep")) {
		h_samples[0] = (TH1D*)h->Clone();
                h_samples[0]->SetFillColor(kRed);
                h_samples[0]->SetLineColor(kRed);
		AddOverAndUnderFlow(h_samples[0]);
	}
	else if(fIn.Contains("diLep")){
		h_samples[1] = (TH1D*)h->Clone();
                h_samples[1]->SetFillColor(kGreen);
                h_samples[1]->SetLineColor(kGreen);
		AddOverAndUnderFlow(h_samples[1]);
	}	
	else if(fIn.Contains("_top")){
		h_samples[2] = (TH1D*)h->Clone();
                h_samples[2]->SetFillColor(kAzure);
                h_samples[2]->SetLineColor(kAzure);
		AddOverAndUnderFlow(h_samples[2]);
	}
	else if(fIn.Contains("_antitop")){
		h_samples[3] = (TH1D*)h->Clone();
                h_samples[3]->SetFillColor(kYellow);
                h_samples[3]->SetLineColor(kYellow);
		AddOverAndUnderFlow(h_samples[3]);
	}
    } 

    void printHisto() {

	for (int s = 3; s >= 0; s--) {
	//for (int s = 0; s < 4; s++) 
		myStack->Add(h_samples[s]);
	}

	
	TCanvas * canvas = new TCanvas(m_name,m_name,632,112,500,502);
        canvas->SetHighLightColor(2);
        canvas->Range(0,0,1,1);
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetBorderSize(2);
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->SetLeftMargin(0.20);
        canvas->SetRightMargin(0.10);
        canvas->SetTopMargin(0.08);
        canvas->SetBottomMargin(0.13);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameBorderMode(0);

        myStack->Draw("histo");

        myStack->GetXaxis()->SetTitle(m_title);
        myStack->GetXaxis()->SetLabelFont(42);
        myStack->GetXaxis()->SetLabelOffset(0.007);
        myStack->GetXaxis()->SetLabelSize(0.045);
        myStack->GetXaxis()->SetTitleSize(0.06);
        myStack->GetXaxis()->SetTitleOffset(0.95);
        myStack->GetXaxis()->SetTitleFont(42);
	TString yTitle = h_samples[0]->GetYaxis()->GetTitle();
        myStack->GetYaxis()->SetTitle(yTitle);
        myStack->GetYaxis()->SetLabelFont(42);
        myStack->GetYaxis()->SetLabelSize(0.045);
        myStack->GetYaxis()->SetTitleSize(0.06);
        myStack->GetYaxis()->SetTitleOffset(1.45);
        myStack->GetYaxis()->SetTitleFont(42);
	
	TLatex *   tex = new TLatex(0.955,0.955,"16.8 fb^{-1} (2016, 13 TeV)");
	tex->SetNDC();
	tex->SetTextAlign(31);
	tex->SetTextFont(42);
	tex->SetTextSize(0.048);
	tex->SetLineWidth(2);
	tex->Draw();

        TLegend *leg = new TLegend(0.6129032,0.5035989,0.7983871,0.9156118,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.0446761);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);


        TLegendEntry * entry=leg->AddEntry("NULL","t#bar{t}#rightarrow l+jets","f");
        entry->SetFillColor(kRed);
        entry->SetFillStyle(1001);
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);

        entry=leg->AddEntry("NULL","t#bar{t}#rightarrow 2l+jets","f");
        entry->SetFillColor(kGreen);
        entry->SetFillStyle(1001);
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);

        entry=leg->AddEntry("NULL","#it{t}-channel","f");
        entry->SetFillColor(kAzure);
        entry->SetFillStyle(1001);
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);

        entry=leg->AddEntry("NULL","#it{#bar{t}}-channel","f");
        entry->SetFillColor(kYellow);
        entry->SetFillStyle(1001);
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);

        leg->Draw();
	/*tex = new TLatex(0.8684108,4612.434,"CMS");
	tex->SetTextAlign(31);
	tex->SetTextSize(0.05956813);
	tex->SetLineWidth(2);
	tex->Draw();
	tex = new TLatex(1.421617,4179.812,"Preliminary");
	tex->SetTextAlign(31);
	tex->SetTextFont(52);
	tex->SetTextSize(0.048);
	tex->SetLineWidth(2);
	tex->Draw();
                                                                             
	pad1->Modified();
        canvas->cd();
	*/
	canvas->cd();
	canvas->Modified();
	canvas->cd();

        canvas->SaveAs((m_path+(std::string)m_name+".png").c_str());
        canvas->SaveAs((m_path+(std::string)m_name+".C").c_str());
    }
    void makePlotMCTruth(TH1 * h) {

	gStyle->SetOptStat(0);
	h_MCTruth = (TH1D*)h->Clone();

	TCanvas * canvas = new TCanvas(m_name,m_name,632,112,500,502);
        canvas->SetHighLightColor(2);
        canvas->Range(0,0,1,1);
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetBorderSize(2);
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->SetLeftMargin(0.20);
        canvas->SetRightMargin(0.10);
        canvas->SetTopMargin(0.08);
        canvas->SetBottomMargin(0.13);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameBorderMode(0);


        h_MCTruth->GetXaxis()->SetTitle(m_title);
        h_MCTruth->GetXaxis()->SetLabelFont(42);
        h_MCTruth->GetXaxis()->SetLabelOffset(0.007);
        h_MCTruth->GetXaxis()->SetLabelSize(0.045);
        h_MCTruth->GetXaxis()->SetTitleSize(0.06);
        h_MCTruth->GetXaxis()->SetTitleOffset(0.95);
        h_MCTruth->GetXaxis()->SetTitleFont(42);
	TString yTitle = h_MCTruth->GetYaxis()->GetTitle();
        h_MCTruth->GetYaxis()->SetTitle(yTitle);
        h_MCTruth->GetYaxis()->SetLabelFont(42);
        h_MCTruth->GetYaxis()->SetLabelSize(0.045);
        h_MCTruth->GetYaxis()->SetTitleSize(0.06);
        h_MCTruth->GetYaxis()->SetTitleOffset(1.45);
        h_MCTruth->GetYaxis()->SetTitleFont(42);
	
        h_MCTruth->SetFillColor(kRed);
        h_MCTruth->SetLineColor(kRed);
        h_MCTruth->SetStats();
	canvas->cd();
	
        h_MCTruth->Draw("hist");

	TLatex *   tex = new TLatex(0.955,0.955,"16.8 fb^{-1} (2016, 13 TeV)");
	tex->SetNDC();
	tex->SetTextAlign(31);
	tex->SetTextFont(42);
	tex->SetTextSize(0.048);
	tex->SetLineWidth(2);
	tex->Draw();

        TLegend *leg = new TLegend(0.61,0.71,0.71,0.75,NULL,"brNDC");
        //TLegend *leg = new TLegend(0.6129032,0.5035989,0.7983871,0.9156118,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        //leg->SetTextSize(0.0446761);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);


        TLegendEntry * entry=leg->AddEntry("NULL","t#bar{t}#rightarrow l+jets","f");
        entry->SetFillColor(kRed);
        entry->SetFillStyle(1001);
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(0.01);
        entry->SetTextFont(42);
	
        leg->Draw();
	/*tex = new TLatex(0.8684108,4612.434,"CMS");
	tex->SetTextAlign(31);
	tex->SetTextSize(0.05956813);
	tex->SetLineWidth(2);
	tex->Draw();
	tex = new TLatex(1.421617,4179.812,"Preliminary");
	tex->SetTextAlign(31);
	tex->SetTextFont(52);
	tex->SetTextSize(0.048);
	tex->SetLineWidth(2);
	tex->Draw();
                                                                             
	pad1->Modified();
        canvas->cd();
	*/
	canvas->cd();
	canvas->Modified();
	canvas->cd();

        canvas->SaveAs((m_path+(std::string)m_name+".png").c_str());
        canvas->SaveAs((m_path+(std::string)m_name+".C").c_str());
    }

    void makeCutFlowTable() {
	
	if (inCutFlow) {

	        for (int s = 0; s < 4; s++) {
		cout<<"calling makeCutFlowTable..."<<endl;	
                cout<<fixed<<" --- entries: "<<h_samples[s]->GetName()<<" : "<<h_samples[s]->GetEntries()<<endl;
		cout<<endl;
                cout<<fixed<<" --- integral: "<<h_samples[s]->GetName()<<" : "<<h_samples[s]->Integral()<<endl;
		cout<<endl;
		double err;
		h_samples[s]->IntegralAndError(1,m_nBins,err);
		cout<<"	sample name: "<<h_samples[s]->GetName()<<",   integral: "<<h_samples[s]->Integral()<<" +- "<<err<<endl;
        }


	}		

    }	  
};

void ToMaAn_plotter() {

    std::vector<TString > inputFiles;

    // semi-lep ttbar, signal MC sample
    inputFiles.push_back("/home/zeinlai/TopMassMeasurement/inputfiles/M172p5.root");
    // di-lep ttbar
    inputFiles.push_back("/home/zeinlai/TopMassMeasurement/inputfiles/TTTo2L2Nu_172p5.root");
    // single-top t-channel
    inputFiles.push_back("/home/zeinlai/TopMassMeasurement/inputfiles/ST_t-channel_top_172p5.root");
    // single-antitop t-channel
    inputFiles.push_back("/home/zeinlai/TopMassMeasurement/inputfiles/ST_t-channel_antitop_172p5.root");
    
    std::vector<TString> sampleNames = {"SemiLepTTbar","diLepTTbar","ST_t-channel_top","ST_t-channel_antitop"};


    const double lumi = 16.81*1000; // in units of pb^-1, second part of 2016 dataset analyzed
    std::vector<double> xsec = {833.9*4/9,833.9/9,134.2,80.};

    double yields[4][4];
    for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 4; j++) {
		yields[i][j] = 0.0;
	}
    }	

    string pathPNG = "ToMaAnPlots/";
    mkdir(pathPNG.c_str(),0777);

    std::vector<string> cutFlowList = {"mcTruth/", "initial/", "hardMuSelection/", "jetSelection/", "softMuSelection/"};
    std::vector<string> combined_path;
    for (unsigned int p = 0; p < cutFlowList.size(); p++) {
	    combined_path.push_back(pathPNG+cutFlowList[p]);
	    mkdir((pathPNG+cutFlowList[p]).c_str(),0777);
    }			
/*
  KEY: TH1F	DR_recMu_vs_associatedGenPar_inOtherTTbar;1
  KEY: TH1F	DR_recMu_vs_associatedGenPar_inSemiMuTTbar;1
  KEY: TH1I	RecMu_associatedGenPartFlav_inOtherTTbar;1
  KEY: TH1I	RecMu_associatedGenPartFlav_inSemiMuTTbar;1

#mctruth
  KEY: TH1F	dPt_softVsMatchedGenMu_OverPt;1
  KEY: TH1F	DR_recMu_vs_genMu_inSemiMuTTbar;1
  */
    std::vector<MakePlot > allHistos;

    MakePlot nrEles("nEles","N_{e_{hard}^{gen}}",combined_path[0]); 
    MakePlot nrMus("nMus","N_{#mu_{hard}^{gen}}",combined_path[0]); 
    MakePlot nrTaus("nTaus","N_{#tau_{hard}^{gen}}",combined_path[0]);
    MakePlot nrSemiMuBDecays("nSemiMuBDecays","N_{#mu_{soft}^{gen}} (before selection)",combined_path[0]); 
    MakePlot nrGenSoftMus("nGenSoftMus","N_{#mu_{soft}^{gen}} (after jet selection)",combined_path[0]); 
    MakePlot MlMu_MCTruth_sameTop("MlMu_MCTruth_sameTop","M_{(#mu_{hard}^{gen}, #mu_{soft}^{gen})} from same #it{t}",combined_path[0]); 
    MakePlot MlMu_MCTruth_oppositeTop("MlMu_MCTruth_oppositeTop","M_{(#mu_{hard}^{gen}, #mu_{soft}^{gen})} from opposite #it{t}",combined_path[0]); 
    MakePlot nrRecMus("nRecMus","N_{#mu_{hard}^{rec}}",combined_path[1], true); 
    MakePlot nrGoodRecMus("nGoodRecMus","N_{good #mu_{hard}^{rec}}",combined_path[1]); 
    MakePlot muPt("RecMu_pt","p_{T}^{#mu_{hard}^{rec}}",combined_path[1]); 
    MakePlot muEta("RecMu_eta","#eta^{#mu_{hard}^{rec}}",combined_path[1]); 
    MakePlot muPhi("RecMu_phi","#phi^{#mu_{hard}^{rec}}",combined_path[1]); 
    MakePlot isTight("RecMu_isTight","isTight",combined_path[1]); 
    MakePlot isTracker_beforeTightId("RecMu_isTracker_beforeTightId","isTracker",combined_path[1]); 
    MakePlot isTracker("RecMu_isTracker","isTracker",combined_path[1]); 
    MakePlot isGlobal_beforeTightId("RecMu_isGlobal_beforeTightId","isGlobal",combined_path[1]); 
    MakePlot isGlobal("RecMu_isGlobal","isGlobal",combined_path[1]); 
    MakePlot nTrackerLayers_beforeTightId("RecMu_nTrackerLayers_beforeTightId","nTrackerLayers",combined_path[1]); 
    MakePlot nTrackerLayers("RecMu_nTrackerLayers","nTrackerLayers",combined_path[1]); 
    MakePlot nStations_beforeTightId("RecMu_nStations_beforeTightId","nStations",combined_path[1]); 
    MakePlot nStations("RecMu_nStations","nStations",combined_path[1]); 
    MakePlot segmentComp_beforeTightId("RecMu_segmentComp_beforeTightId","segmentComp",combined_path[1]); 
    MakePlot segmentComp("RecMu_segmentComp","segmentComp",combined_path[1]); 
    MakePlot dxy_beforeTightId("RecMu_dxy_beforeTightId","dxy (#it{cm})",combined_path[1]); 
    MakePlot dxy("RecMu_dxy","dxy (#it{cm})",combined_path[1]); 
    MakePlot dz_beforeTightId("RecMu_dz_beforeTightId","dz (#it{cm})",combined_path[1]); 
    MakePlot dz("RecMu_dz","dz (#it{cm})",combined_path[1]); 
    MakePlot relIso("RecMu_relIso","relIso",combined_path[1]); 
    MakePlot jetIdx("RecMu_jetIdx","jetIdx",combined_path[1]); 
    MakePlot nrGoodJets("nGoodJets","N_{good jets}",combined_path[2], true); 
    MakePlot nrBJets("nrBTaggedJets","N_{good #it{b}-jets}",combined_path[2]); 
    MakePlot jetId("jets_jetId","jetId",combined_path[2]);
    MakePlot DeltaRMuJets("DR_goodMu_jets","#Delta R(good #mu_{hard}^{rec},jets)",combined_path[2]); 
    MakePlot nConstituents("jets_nConstituents","nConstituents",combined_path[2]);
    MakePlot muEF("jets_muEF","muEF",combined_path[2]); 
    MakePlot chEmEF("jets_chEmEF","chEmEF",combined_path[2]); 
    MakePlot neEmEF("jets_neEmEF","neEmEF",combined_path[2]); 
    MakePlot chHEF("jets_chHEF","chHEF",combined_path[2]); 
    MakePlot neHEF("jets_neHEF","neHEF",combined_path[2]); 
    MakePlot nrSoftMus("nSoftMus","N_{#mu_{soft}^{rec}}",combined_path[3], true);
    MakePlot looseId("softMu_looseId","#mu_{soft} looseId",combined_path[3]); 
    MakePlot mediumId("softMu_mediumId","#mu_{soft} mediumId",combined_path[3]); 
    MakePlot softMvaId("softMu_softMvaId","#mu_{soft} softMvaId",combined_path[3]); 
    MakePlot DeltaRsoftMuAssJets("DR_softMu_associatedJet","#Delta R(#mu_{soft}^{rec}, associated-jet)",combined_path[3]); 
    MakePlot ptOverJetPt("matchedSoftMu_ptOverJetPt","p_{T}^{#mu_{soft}^{gen-matched}}/p_{T}^{associated-jet}",combined_path[3]); 
    MakePlot associatedGenPartFlav("softMu_associatedGenPartFlav","associatedGenPartFlav_{#mu_{soft}^{selected}}",combined_path[3]); 
    MakePlot softMuPt("softRecMu_pt","p_{T}^{#mu_{soft}^{rec}}",combined_path[4]); 
    MakePlot DeltaR_softRecMu_softGenMu("DR_softRecMu_softGenMu","#Delta R(#mu_{soft}^{rec}, #mu_{soft}^{gen})",combined_path[4]); 
    MakePlot M_lMu("MlMu","M_{(#mu_{hard}^{rec}, #mu_{soft}^{rec})}",combined_path[4], true); 
    MakePlot passedEvts("passedEvts","passedEvts",combined_path[4]); 
    allHistos.push_back(nrEles);
    allHistos.push_back(nrMus);
    allHistos.push_back(nrTaus);
    allHistos.push_back(nrSemiMuBDecays);
    allHistos.push_back(nrGenSoftMus);
    allHistos.push_back(MlMu_MCTruth_sameTop);
    allHistos.push_back(MlMu_MCTruth_oppositeTop);
    allHistos.push_back(nrRecMus);
    allHistos.push_back(nrGoodRecMus);
    allHistos.push_back(muPt);
    allHistos.push_back(muEta);
    allHistos.push_back(muPhi);
    allHistos.push_back(isTight);
    allHistos.push_back(isTracker_beforeTightId);
    allHistos.push_back(isTracker);
    allHistos.push_back(isGlobal_beforeTightId);
    allHistos.push_back(isGlobal);
    allHistos.push_back(nTrackerLayers_beforeTightId);
    allHistos.push_back(nTrackerLayers);
    allHistos.push_back(nStations_beforeTightId);
    allHistos.push_back(nStations);
    allHistos.push_back(segmentComp_beforeTightId);
    allHistos.push_back(segmentComp);
    allHistos.push_back(dxy_beforeTightId);
    allHistos.push_back(dxy);
    allHistos.push_back(dz_beforeTightId);
    allHistos.push_back(dz);
    allHistos.push_back(relIso);
    allHistos.push_back(jetIdx);
    allHistos.push_back(nrGoodJets);
    allHistos.push_back(nrBJets);
    allHistos.push_back(jetId);
    allHistos.push_back(DeltaRMuJets);
    allHistos.push_back(nConstituents);
    allHistos.push_back(muEF);
    allHistos.push_back(chEmEF);
    allHistos.push_back(neEmEF);
    allHistos.push_back(chHEF);
    allHistos.push_back(neHEF);
    allHistos.push_back(nrSoftMus);
    allHistos.push_back(looseId);
    allHistos.push_back(mediumId);
    allHistos.push_back(softMvaId);
    allHistos.push_back(DeltaRsoftMuAssJets);
    allHistos.push_back(ptOverJetPt);
    allHistos.push_back(associatedGenPartFlav);
    allHistos.push_back(softMuPt);
    allHistos.push_back(DeltaR_softRecMu_softGenMu);
    allHistos.push_back(M_lMu);
    allHistos.push_back(passedEvts);

    std::vector<MakePlot > MCTruth_Histos;
    MakePlot lepWmass("lepW_M","M_{leptonic W}",combined_path[0]); 
    MakePlot hadWmass("hadW_M","M_{hadronic W}",combined_path[0]); 
    MakePlot lepTmass("lepTop_M","M_{leptonic #it{t}}",combined_path[0]); 
    MakePlot hadTmass("hadTop_M","M_{hadronic #it{t}}",combined_path[0]); 
    MCTruth_Histos.push_back(lepWmass);
    MCTruth_Histos.push_back(hadWmass);
    MCTruth_Histos.push_back(lepTmass);
    MCTruth_Histos.push_back(hadTmass);
    
    int cut = -1;
    for (unsigned int s = 0; s < allHistos.size(); s++) {
     if (allHistos[s].useInMCTable()) cut++;
     //else continue;
     cout<<"	\n+++++ starting plot production for "<<allHistos[s].getName()<<endl;

     for (unsigned int file = 0; file < inputFiles.size(); file++) {

	TString fIn = inputFiles[file];
	TFile *f = TFile::Open(fIn);
	TH1I * h =(TH1I*)f->Get("nEles");
	double nEntries = h->GetEntries();

	const double weight = xsec[file]*lumi/nEntries;
	cout<<" --- sample name: "<<sampleNames[file]<<" , and weight: "<<weight<<endl;

	auto k =(TH1D*)f->Get(allHistos[s].getName());
	allHistos[s].setNBins(k);
	k->Scale(weight);
	allHistos[s].fillHisto(k,sampleNames[file]);
	
	if (allHistos[s].useInMCTable()) {
		cout<<"@@@@@@@@@@@@@ integral before overandunder "<<k->Integral()<<endl;
		AddOverAndUnderFlow(k);
		yields[cut][file] = k->Integral();
		cout<<"@@@@@@@@@@@@@ integral after overandunder "<<k->Integral()<<endl;
	}
     }

    allHistos[s].printHisto();
    allHistos[s].makeCutFlowTable();
    
    }
    TString cutList[4]  = {"initial", "==1Hard-Muon", ">=4Jets & >=1B-jets", "==1Soft-Muon"};
    double totYields[4];
		cout<< "\n\n"<< std::string(100, '=') <<endl;
		cout<<"cut flow table"<<endl;
		cout<< std::string(100, '=') <<endl;
	for (unsigned int cut = 0; cut < 4; cut++) {
		totYields[cut] = 0.;
		cout<<cutList[cut]<<" | ";
		for (unsigned int sample = 0; sample < 4; sample++) {
			totYields[cut] += yields[cut][sample];
			cout<<setprecision(5)<<sampleNames[sample]<<"  "<<fixed<<yields[cut][sample]<<" | ";
		}
		cout<<" \n \t\t(cut purity: "<<yields[cut][0]/totYields[cut]<<")\n";
	}
		cout<<" [Definition of cut purity : the contamination of signal events at each step of cut selection]"<<endl;
		cout<< "\n\n"<<std::string(100, '=') <<endl;
	        cout<<"calculating the cut efficiency (i.e. the fraction of survived events w.r.t. the previous step)"<<endl;
		cout<< std::string(100, '=') <<endl;
	        for (unsigned int cut = 0; cut < 4; cut++) {
                cout<<cutList[cut]<<" | ";

		for (unsigned int sample = 0; sample < 4; sample++) {
			if (cut!=0)
			cout<<setprecision(5)<<sampleNames[sample]<<"  "<<yields[cut][sample]/yields[cut-1][sample]<<" | ";
			else  cout<<" 1. | ";
		}
		                cout<<endl;
		}

		                cout<<"\n\n"<<endl;
					cout<<"total selection efficiency (%): "<<endl;
				for (unsigned int sample = 0; sample < 4; sample++) {
					cout<<setprecision(2)<<sampleNames[sample]<<"  "<<yields[3][sample]/yields[0][sample]*100<<" | ";
				}
                                cout<<"\n\n"<<endl;

    for (unsigned int file = 0; file < inputFiles.size(); file++) {
	if (!(sampleNames[file].Contains("SemiLep"))) continue;
	cout<<"looping over MCTruth_Histos.size "<<endl;
	TString fIn = inputFiles[file];
    	for (unsigned int s = 0; s < MCTruth_Histos.size(); s++) {

		TFile *f = TFile::Open(fIn);
		TH1I * h =(TH1I*)f->Get("nEles");
		double nEntries = h->GetEntries();

		const double weight = xsec[file]*lumi/nEntries;

		auto k =(TH1D*)f->Get(MCTruth_Histos[s].getName());
		MCTruth_Histos[s].setNBins(k);
		k->Scale(weight);
		MCTruth_Histos[s].makePlotMCTruth(k);
    
    	}
    }
}
