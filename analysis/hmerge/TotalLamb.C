{
    gStyle->SetOptStat(0);
    TChain* ch2dat = new TChain("h2"); //with pi0
    ch2dat -> Add("*.root");
    
    double lend=1.101, rend=1.131, MLambda=1.115683;
    int Nbins=100;
    TCanvas *c1 = new TCanvas("c1","Lambda invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","#Lambda invariant mass",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    ch2dat -> Draw("ml>>+hdat","","goff"); //rm*rm<0.25
    Ntot = hdat -> GetEntries();
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+%lf*[3]*TMath::Gaus(x,[4],[5],true)+%lf*[6]*TMath::Gaus(x,[7],[8],true)+[9]",binw,binw,binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+%lf*[3]*TMath::Gaus(x,[4],[5],true)+%lf*[6]*TMath::Gaus(x,[7],[8],true)",binw,binw,binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]",lend,rend); 
    
    
    
    fdat -> SetParameters(20000*3,MLambda,0.0001,20000*3,MLambda,0.001,200*3,MLambda,0.01,4000*3);
    fdat -> SetParLimits(1,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(4,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(7,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(0,0,1e7);
    fdat -> SetParLimits(3,0,1e7);
    fdat -> SetParLimits(6,0,1e7);

    
    
    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,8,covSignal);
    covFit.GetSub(9,9,covBackground);
    double * par = fdat -> GetParameters();

    //fitResult -> PrintCovMatrix(cout);
    //covFit.Print();
    //covSignal.Print();
    //covBackground.Print();
    
    
    fsig -> SetParameters(par);
    fbkg -> SetParameter(0,fdat->GetParameter(9));
 
    
   /*
    * Nsig = (fdat->GetParameter(0)+fdat->GetParameter(3)+fdat->GetParameter(6));
    dNsig = sqrt((fdat->GetParError(0)*fdat->GetParError(0)+fdat->GetParError(3)*fdat->GetParError(3)+fdat->GetParError(6)*fdat->GetParError(6)));
    *
    * Nbkg = fbkg -> GetParameter(0)*Nbins;
    dNbkg = fbkg -> GetParError(0)*Nbins;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda naive: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background naive: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
    */
    
    Nsig = fsig -> Integral(lend,rend)/binw; 
    dNsig = fsig -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+9,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
    
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
  
    
    hdat -> GetXaxis()->SetTitle("M(#Lambda), GeV");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.4f",binw));
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(3);
    hdat -> DrawCopy("ep");
    
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(3);
    fbkg -> DrawCopy("same");
    
    /*
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(3);
    fsig -> Draw("same");
    */
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(3);
    fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
    //leg->AddEntry("fsig","Signal","l");
    leg->AddEntry("fbkg","Background","l");
    leg->AddEntry("fdat","Signal + background","l");
    leg->Draw("same");

    
} 
