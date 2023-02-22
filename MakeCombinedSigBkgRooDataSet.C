#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <ios>
#include <fstream>
#include "TH2.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TColorWheel.h"
#include "TColorGradient.h"
#include "RooPlot.h"
#include "RooRandom.h"

using namespace std;
using namespace RooFit;

// ---- this script provides the instruction to add weighted signal plus background events ----
// ---- into roodatasets, as components of simultaneous fit, before calling the fit process ---


void MakeCombinedSigBkgRooDataSet() {

   	RooRealVar * mass = new RooRealVar("ToMaAn_MlMu", "ToMaAn_MlMu", 0, 250); 
	RooCategory sampleWeighted("sampleWeighted", "sampleWeighted");	
	std::vector<TTree*> trees;
	std::vector<RooFormulaVar*> weight_semilepttbar;
	std::vector<RooRealVar*> wVars;
	std::vector<RooDataSet*> datas;	
	std::vector<RooDataSet*> Wdatas;	
	stringstream fname;
        std::vector<TString> tree_names;


	std::vector<RooDataSet*> bkgdatas;	

        std::vector<TString> bkgFiles;
	bkgFiles.push_back("TTTo2L2Nu_34.root");
	bkgFiles.push_back("ST_t-channel_top_18.root");
	bkgFiles.push_back("ST_t-channel_antitop_15.root");
	bkgFiles.push_back("ST_tW_top_9.root");
	bkgFiles.push_back("ST_tW_antitop_10.root");

	std::vector<double> bkg_expected;
	bkg_expected.push_back(833.9*38.25*1000/9); // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
	bkg_expected.push_back(134.2*38.25*1000); // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef
	bkg_expected.push_back(80.*38.25*1000);
	bkg_expected.push_back(79.3*38.25*1000/2);
	bkg_expected.push_back(79.3*38.25*1000/2);


	
	// --- reading bkg root files ---	
	cout<<"\nstarting bkg inclusion ..."<<endl;

	for (unsigned int b = 0; b < bkgFiles.size(); b++) {

	
        TFile * f = TFile::Open(bkgFiles[b]);

	TH1I * h =(TH1I*)f->Get("nRecMus");
        double nEntries = h->GetEntries();
	RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);
	cout << fixed;
	cout<<"________ "<<bkgFiles[b]<<" sample: "<<nGen->getVal()<<" and nEntries: "<<nEntries<<",   weight should be: "<<bkg_expected[b]/nEntries<<endl;
        TTree* tree_bkg = (TTree*) f->Get("Events");
	RooRealVar * exp = new RooRealVar("exp","exp",bkg_expected[b]);
        RooFormulaVar * wgt = new RooFormulaVar("wgt", "(@0/@1)", RooArgSet(*exp,*nGen));

	RooDataSet * dataset = new RooDataSet("dataset", "dataset",tree_bkg, *mass, "");
	RooRealVar * weight_bkg = (RooRealVar*)dataset->addColumn(*wgt);
	RooDataSet * weightedData = new RooDataSet("weightedData", "weightedData",dataset,*dataset->get(),0, weight_bkg->GetName());

	bkgdatas.push_back(new RooDataSet(bkgFiles[b],bkgFiles[b],dataset,*dataset->get(),0, weight_bkg->GetName()));

	}


	// --- reading sig root files ---	

	cout<<"\nloop over signal samples ..."<<endl;

	for(int i = 0; i < 5; i++){
		if (i==2) {
		for(int j = -1; j < 2; j++){	
			if (j == 0) continue;
		fname.str("");
        	fname << "m" << 166 + (3*i) + j <<"p5.root";
                TFile * f = TFile::Open(fname.str().c_str());
		TH1I * h =(TH1I*)f->Get("nEles");
                double nEntries = h->GetEntries();
		RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);
		cout << fixed;
		cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nEntries: "<<(int)nEntries<<",   weight should be: "<<833.9*38.25*1000*4/(9*nEntries)<<endl;
                trees.push_back((TTree*) f->Get("Events"));
		fname.str("");
        	fname << 166 + (3*i) + j <<"p5";
		tree_names.push_back(fname.str().c_str());
                fname.str("");
                fname << "Weights" << 166 + (3*i) + j;
                weight_semilepttbar.push_back(new RooFormulaVar(fname.str().c_str(), "(833.9*38.25*1000*4/(9*@0))", RooArgSet(*nGen)));
                fname.str("");
                fname << "sample"<<     166 + (3*i) + j << "p5";
                sampleWeighted.defineType(fname.str().c_str());
		}
		}
		else {
		fname.str("");
		fname << "m" << 166 + (3*i) <<"p5.root";
		TFile * f = TFile::Open(fname.str().c_str());
		TH1I * h =(TH1I*)f->Get("nEles");
                double nEntries = h->GetEntries();
		RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);
		cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nEntries: "<<(int)nEntries<<",   weight should be: "<<833.9*38.25*1000*4/(9*nEntries)<<endl;
		trees.push_back((TTree*) f->Get("Events"));
		fname.str("");
        	fname << 166 + (3*i) <<"p5";
		tree_names.push_back(fname.str().c_str());
		fname.str("");	
		fname << "Weights" << 166 + (3*i);
                weight_semilepttbar.push_back(new RooFormulaVar(fname.str().c_str(), "(833.9*38.25*1000*4/(9*@0))", RooArgSet(*nGen)));
		fname.str("");
		fname << "sample"<<	166 + (3*i) << "p5";
		sampleWeighted.defineType(fname.str().c_str());
	}
}


	cout<<"\nloop over trees ..."<<endl;
	for(unsigned int i = 0; i < trees.size(); i++){
		fname.str("");
		fname << "data" << tree_names[i];
                cout<<"!!!!!!!! "<<fname.str().c_str()<<endl;		
		datas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),trees[i], *mass, ""));
		datas[datas.size()-1]->Print();
		wVars.push_back((RooRealVar*)datas[datas.size()-1]->addColumn(*weight_semilepttbar[i]));
		cout<<" weight_semilepttbar[i]:  "<<weight_semilepttbar[i]->evaluate()<<endl;
		Wdatas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),datas[datas.size()-1],*datas[datas.size()-1]->get(),0, wVars[wVars.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVars[wVars.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVars[wVars.size()-1]->GetName()<<endl;
		Wdatas[Wdatas.size()-1]->Print();		
		//Wdatas[i]->append(*dilep_dataset);
		//Wdatas[i]->append(*bkgdatas[0]);
		//Wdatas[Wdatas.size()-1]->Print();		
                for (unsigned int b = 0; b < bkgFiles.size(); b++) {
                        cout<<"  ... appending "<<bkgFiles[b]<<endl;
                        Wdatas[i]->append(*bkgdatas[b]);
                        Wdatas[Wdatas.size()-1]->Print();
                }

	}

	
	cout << "\nbegin making the combined weigted RooDataSet..." << endl;


	 RooDataSet combDataWeighted("combDataWeighted", "combined data", RooArgSet(*mass,*wVars[0]), Index(sampleWeighted),
            Import("sample166p5", *Wdatas[0]),WeightVar("Weights166"));
	 RooDataSet combData169p5("combDataWeighted169p51", "combined data", RooArgSet(*mass,*wVars[1]), Index(sampleWeighted),
            Import("sample169p5", *Wdatas[1]),WeightVar("Weights169"));
	 RooDataSet combData171p5("combDataWeighted171p5", "combined data", RooArgSet(*mass,*wVars[2]), Index(sampleWeighted),
            Import("sample171p5", *Wdatas[2]),WeightVar("Weights171"));
	 RooDataSet combData173p5("combDataWeighted173p51", "combined data", RooArgSet(*mass,*wVars[3]), Index(sampleWeighted),
            Import("sample173p5", *Wdatas[3]),WeightVar("Weights173"));
	 RooDataSet combData175p5("combDataWeighted175p5", "combined data", RooArgSet(*mass,*wVars[4]), Index(sampleWeighted),
            Import("sample175p5", *Wdatas[4]),WeightVar("Weights175"));
	 RooDataSet combData178p5("combDataWeighted178p51", "combined data", RooArgSet(*mass,*wVars[5]), Index(sampleWeighted),
            Import("sample178p5", *Wdatas[5]),WeightVar("Weights178"));

	combDataWeighted.append(combData169p5);
	combDataWeighted.append(combData171p5);
	combDataWeighted.append(combData173p5);
	combDataWeighted.append(combData175p5);
	combDataWeighted.append(combData178p5);

    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    combDataWeighted.Print();
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    
}

