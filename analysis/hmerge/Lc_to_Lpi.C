{
    gStyle->SetOptStat(0);
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add("*.root");
    
    double lend=2.21, rend=2.36, MLambdac=2.28646;
    int Nbins=100;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","#Lambda_{c} invariant mass in #Lambda_{c} #rightarrow #Lambda#pi decay",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    TCut Mwindow = Form("mlc > %lf && mlc < %lf",lend,rend);
    Ntot = ch1dat -> Draw("mlc>>+hdat","lcch == 1 && ml>1.105 && ml<1.125 && rm*rm<1.44"+Mwindow,"goff"); //"lcch == 1 && ml>1.105 && ml<1.125 && rm*rm<0.25"
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.3)",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.3)",lend,rend);    
    

    fdat -> SetParameters(1000,MLambdac,0.1,3000,-700);
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,4,covBackground);
    double * par = fdat -> GetParameters();

   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);

    
    Nsig = fdat->GetParameter(0);
    dNsig = fdat->GetParError(0);
    
   
    
    Nsig = fsig -> Integral(lend,rend)/binw; 
    dNsig = fsig -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()->SetTitle("M(#Lambda_{c}), GeV");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.4f",binw));
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
    
    
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(3);
    fsig -> Draw("same");
    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(3);
    fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
    leg->AddEntry("fsig","Signal","l");
    leg->AddEntry("fbkg","Background","l");
    leg->AddEntry("fdat","Signal + background","l");
    leg->Draw("same");

    
} 
 
