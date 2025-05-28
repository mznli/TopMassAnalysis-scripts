void CLT(){

	std::vector<TH1D *> h ;
	std::vector<TF1 *> f ;
	double size[4] = {1,2,10,30};
	double avg = 0;
	char name[10];
	for(int l=0; l<4; l++) {
		//sprintf(name,"histo %d",l);
		h.push_back(new TH1D("CLT visual demonstration", "CLT visual demonstration", 400, 0., 20.));
		h[l]->SetLineWidth(2);
		f.push_back(new TF1("f","[0]*TMath::Gaus(x,6*sqrt(3),[1])",0.,20.));
		f[l]->SetParameter(0,1);				
		f[l]->SetParameter(1,2./sqrt(size[l]));				
		double func_norm = f[l]->Integral(-TMath::Infinity(),TMath::Infinity());
		
		for(int i=0; i<10000; i++) {
	
			if(l == 0) h[0]->Fill(gRandom->Uniform(4*sqrt(3),8*sqrt(3)));
			if(l == 1) h[1]->Fill((gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)))/2.);
                        if(l == 2) h[2]->Fill((gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)))/10.);

                        if(l == 3) h[3]->Fill((gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)) + gRandom->Uniform(4*sqrt(3),8*sqrt(3)))/30.);


		}
		f[l]->SetParameter(0,h[l]->GetEntries()*h[l]->GetBinWidth(1)/func_norm);
		f[l]->SetLineColor(l+2);
	}
h[3]->SetLineColor(41);
h[3]->Draw();	
h[2]->SetLineColor(30);
h[2]->Draw("same");	
h[1]->SetLineColor(43);
h[1]->Draw("same");	
h[0]->SetLineColor(41);
h[0]->Draw("same");	

for(int i=0; i<4; i++)
f[i]->Draw("same");


}
