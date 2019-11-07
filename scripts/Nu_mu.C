{
    gStyle->SetOptStat(0);
    TString datapath = "../analysis/hmerge/";    
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    double lend=-1.5, rend=1.5;
    int Nbins=100;
    TCanvas *c1 = new TCanvas("c1","Muon neutrino invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","Squared recoil mass of #nu_{#mu}",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
	TCut Mwindow = Form("rm*fabs(rm) > %lf && rm*fabs(rm) < %lf",lend,rend);
    Ntot = ch1dat -> Draw("rm*fabs(rm)>>hdat","lcch == 4 && m*m-10.58*10.58+2*10.58*p-rm*rm > -2 && m*m-10.58*10.58+2*10.58*p-rm*rm < 2 && p>0.1 && mlc<2.2"+Mwindow,"goff");  
    //"lcch == 4 && m*m-10.58*10.58+2*10.58*p-rm*rm > -2 && m*m-10.58*10.58+2*10.58*p-rm*rm < 2"
    
   
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+%lf*[3]*TMath::Gaus(x,[4],[5],true)+%lf*[6]*TMath::Gaus(x,[7],[8],true)",binw,binw,binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+%lf*[3]*TMath::Gaus(x,[4],[5],true)",binw,binw),lend,rend);  

    fdat -> SetParameters(1000,0,0.05,1000,0.1,0.5,1000,0.2,1);
    //fdat -> SetParLimits(1,-0.02,0.02);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,8,covBackground);
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


    
    hdat -> GetXaxis()->SetTitle("RM^{2}(#nu_{#mu}), GeV^{2}");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.2f",binw));
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> DrawCopy("ep");

    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(3);
    fbkg -> DrawCopy("same");	    


    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(3);
    fdat-> DrawCopy("same");

    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(3);
    fsig -> Draw("same");
    
    TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
	leg->AddEntry("fsig","Signal","l");    
	leg->AddEntry("fbkg","Background","l");	
	leg->AddEntry("fdat","Signal + background","l");
    leg->Draw("same"); 
    c1 -> Draw();
    
}
 
