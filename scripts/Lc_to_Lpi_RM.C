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
    
    double lend=1.5, rend=3., MLambdac=2.28646; //lend=2., rend=2.6
    int Nbins=30;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1024,768);
    TH1D* hdat = new TH1D("hdat","#Lambda_{c} #rightarrow #Lambda#pi",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","#Lambda_{c} #rightarrow #Lambda#pi sb",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    TCut Mwindow = Form("npi0==0 && rmx > %lf && rmx < %lf",lend,rend);
    Ntot = ch1dat -> Draw("rmx>>+hdat","lcch == 1 && tag!=2 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && abs(mlc-2.2879)<0.0153"+Mwindow,"goff"); //"lcch == 1 && ml>1.1 && ml<1.12"
    ch1dat -> Draw("rmx>>+hsb","lcch == 1 && tag!=2 && abs(ml-1.11568)<0.003 && ((tag!=11 && tag!=12) || abs(ml1-1.11568)<0.003) && abs(mlc-2.2879)>0.0255 && abs(mlc-2.2879)< 0.0408"+Mwindow,"goff");
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.287)+[5]*(x-2.287)^2",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.287)+[2]*(x-2.287)^2",lend,rend);    
    

    fdat -> SetParameters(50,MLambdac,0.1,30,300,20);
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,5,covBackground);
    double * par;
    par = fdat -> GetParameters();
  
   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);
    
    Nsig = fdat -> GetParameter(0);//fsig -> Integral(lend,rend)/binw; 
    dNsig = fdat -> GetParError(0);// -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg = fbkg -> Integral(lend,rend)/binw; 
    dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("RM(X) [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events / ( %.3f )",binw));
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
    
    hsb -> SetLineColor(8);
    hsb -> SetLineWidth(4);
    hsb -> Draw("same");
    
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");
    
    
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
   // fsig -> Draw("same");
    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.17,0.70,0.29,0.89);
    leg -> AddEntry("hdat","data","E P");
    leg -> AddEntry("hsb","sideband","l");
    //leg->AddEntry("fsig","Signal","l");
    leg -> AddEntry("fdat","fit","l");
    leg -> AddEntry("fbkg","background","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.45, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.39, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.33, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   
    
}  
