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
    
    double lend=-0.25, rend=0.25, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=25;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow #Lambda e^{+}#nu_{e}",Nbins,lend,rend);
    TH1D* hrsb = new TH1D("hrsb","right sb",Nbins,lend,rend);
    TH1D* hlsb = new TH1D("hlsb","left sb",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;

    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("rmvis*fabs(rmvis) > %lf && rmvis*fabs(rmvis) < %lf",lend,rend);
    TCut commonLcLpiCut = "lcch==3 && abs(ml-1.11568)<0.003 && ((tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))))";
    TCut rmxlsbCut = "rmx-2.29 < -0.15 && rmx-2.29 > -0.25";
    TCut rmxrsbCut = "rmx-2.29 > 0.15 && rmx-2.29 < 0.25";
    Ntot +=  ch1dat -> Draw("rmvis*fabs(rmvis)>>+hdat","abs(rmx-2.29)<0.1"+ commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("rmvis*fabs(rmvis)>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("rmvis*fabs(rmvis)>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow,"goff");
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]",lend,rend);    
    

    fdat -> SetParameters(10,0,0.05,1); //100,MLambdac,0.01,40,-300
    fdat -> SetParLimits(1,-0.1,+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,3,covBackground);
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
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("M^{2}_{recoil}(#Lambda e^{+}#bar{p}#bar{D}^{*0})||M^{2}_{recoil}(#Lambda e^{+}#bar{p}D^{*-}#pi^{+})[GeV^{2}]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetXaxis()-> SetTitleOffset(1.1);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.2f GeV^{2} ",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.7);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdat -> GetXaxis()->SetTickSize(0.04);
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1.5);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> SetMinimum(0);
    hdat -> Draw("e p");
    
    hrsb -> SetLineColor(4);
    hrsb -> SetLineWidth(4);
    hrsb -> Draw("same");
    
    hlsb -> SetLineColor(8);
    hlsb -> SetLineWidth(4);
    hlsb -> Draw("same");
    
    
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");
    
    
    //fsig -> SetLineColor(4);
    //fsig -> SetLineWidth(4);
    //fsig -> Draw("same");
    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","e p");
   // leg->AddEntry("fsig","Signal","l");
    leg -> AddEntry("fdat","fit","l");
    leg -> AddEntry("fbkg","background","l");
    leg -> AddEntry("hrsb","M(#Lambda#pi^{+}) right sb","l");
    leg -> AddEntry("hlsb","M(#Lambda#pi^{+}) left sb","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   
    
} 
 
