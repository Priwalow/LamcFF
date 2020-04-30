{
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    
    double axisFontSize = 0.065;
    
    TString datapath = "../analysis/hmerge/";    
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    double lend=-0.25, rend=0.25;
    int Nbins=50;
    TCanvas *c1 = new TCanvas("c1","Muon neutrino invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","Squared recoil mass of #nu_{e}",Nbins,lend,rend);
    TH1D* hws = new TH1D("hws","Squared recoil mass of #nu_{e}",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","Squared recoil mass of #nu_{e} sideband",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
	TCut Mwindow = Form("rmvis*fabs(rmvis) > %lf && rmvis*fabs(rmvis) < %lf",lend,rend);
    Ntot = ch1dat -> Draw("rmvis*fabs(rmvis)>>hdat","lcch == 3 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && mlc < 2.13 && plept>0.7 && abs(lepcost)<0.7 && fox>0.2 && plc>1.4 && abs(rmx-2.2901)<0.0438*2"+Mwindow,"goff");  
    ch1dat -> Draw("rmvis*fabs(rmvis)>>hws","lcch == 300 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && mlc < 2.13 && plept>0.7 && abs(lepcost)<0.7 && fox>0.2 && plc>1.4 && abs(rmx-2.2901)<0.0438*2"+Mwindow,"goff");
    ch1dat -> Draw("rmvis*fabs(rmvis)>>hsb","lcch == 3 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && mlc < 2.13 && plept>0.7 && abs(lepcost)<0.7 && fox>0.2 && plc>1.4 && abs(rmx-2.2901)>0.0438*3 && abs(rmx-2.2901)<0.0438*5"+Mwindow,"goff"); 
    
    //"lcch == 4 && m*m-10.58*10.58+2*10.58*p-rm*rm > -2 && m*m-10.58*10.58+2*10.58*p-rm*rm < 2"
    
   
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]",lend,rend);  

    fdat -> SetParameters(100,0,0.05,80);
    fdat -> SetParLimits(1,-0.02,0.02);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,3,covBackground);
    double * par = fdat -> GetParameters();
   	
	/*
	fitResult -> PrintCovMatrix(cout);
    covFit.Print();
    covSignal.Print();
    covBackground.Print();
	*/
    
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);
    
    Nsig = fsig -> Integral(lend,rend)/binw; 
    dNsig = fsig -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;


    cout << "Total number of events: " << Ntot << endl;
    cout << "N nu_mu: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " <<  (int) dNbkg << endl;
    
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;


    
    hdat -> GetXaxis()->SetTitle("RM^{2}(X#Lambdae) [GeV^{2}]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()->SetTitle(Form("Events / ( %.2f )",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.8);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.5);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> Draw("e p");
    
    hws -> SetLineColor(2);
    hws -> SetLineWidth(2);
    //hws -> Draw("same");
    
    hsb -> SetLineColor(8);
    hsb -> SetLineWidth(4);
    hsb -> Draw("same");
    
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
	//leg->AddEntry("fsig","Signal","l");    
	//leg->AddEntry("fdat","Fit","l");
    //leg->AddEntry("fbkg","Background","l");
    //leg->AddEntry("hws","Wrong sign","l");
    leg->AddEntry("hsb","RM(X) sideband","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg->Draw("same"); 
    
     /*  
    fbkg -> SetLineStyle(4);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");	    


    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");

    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
   // fsig -> Draw("same");
 
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
	//leg->AddEntry("fsig","Signal","l");    
	leg->AddEntry("fdat","Fit","l");
    leg->AddEntry("fbkg","Background","l");	
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg->Draw("same"); 
    
  TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2))); */
    
}
 

 
