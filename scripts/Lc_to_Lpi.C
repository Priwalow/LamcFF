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
    
    double lend=2.236, rend=2.336, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=100;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow pK^{-}#pi^{+}",Nbins,lend,rend);
    TH1D* hrsb = new TH1D("hrsb","right sb",Nbins,lend,rend);
    TH1D* hlsb = new TH1D("hlsb","left sb",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    
    double D1_mean = 1.8651, D1_sigma = 0.0064, Dst1_mean = 2.010, Dst1_sigma = 0.0011, rsbwidth1 = 0.00445;
    TCut Dst1_D_window = Form("dch==1 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D1_mean,3*D1_sigma,Dst1_mean,3*Dst1_sigma);
    TCut Dst1_rightsb = Form("dch==1 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D1_mean,3*D1_sigma,Dst1_mean + 4*Dst1_sigma, Dst1_mean + 4*Dst1_sigma + rsbwidth1);
    
    double D2_mean = 1.8651, D2_sigma = 0.0051, Dst2_mean = 2.010, Dst2_sigma = 0.0011, rsbwidth2 = 0.0048;
    TCut Dst2_D_window = Form("dch==2 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D2_mean,3*D2_sigma,Dst2_mean,3*Dst2_sigma);
    TCut Dst2_rightsb = Form("dch==2 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D2_mean,3*D2_sigma,Dst2_mean + 4*Dst2_sigma, Dst2_mean + 4*Dst2_sigma + rsbwidth2);
    
    double D3_mean = 1.8652, D3_sigma = 0.0055, Dst3_mean = 2.010, Dst3_sigma = 0.0012, rsbwidth3 = 0.00542;
    TCut Dst3_D_window = Form("dch==3 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D3_mean,3*D3_sigma,Dst3_mean,3*Dst3_sigma);
    TCut Dst3_rightsb = Form("dch==3 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D3_mean,3*D3_sigma,Dst3_mean + 4*Dst3_sigma, Dst3_mean + 4*Dst3_sigma + rsbwidth3);
    
    double D4_mean = 1.8639, D4_sigma = 0.0098, Dst4_mean = 2.010, Dst4_sigma = 0.0015, rsbwidth4 = 0.00698;
    TCut Dst4_D_window = Form("dch==4 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D4_mean,3*D4_sigma,Dst4_mean,3*Dst4_sigma);
    TCut Dst4_rightsb = Form("dch==4 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D4_mean,3*D4_sigma,Dst4_mean + 4*Dst4_sigma, Dst4_mean + 4*Dst4_sigma + rsbwidth4);

    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("mlc > %lf && mlc < %lf",lend,rend);
    TCut commonLcLpiCut = "lcch==2 && tag==3 && dstch==1";
    TCut rmxlsbCut = "rmx-2.301 < -0.1159 && rmx-2.301 > -0.23172";
    TCut rmxrsbCut = "rmx-2.301 > 0.1159 && rmx-2.301 < 0.23172";
    Ntot +=  ch1dat -> Draw("mlc>>+hdat","abs(rmx-2.301)<0.1159"+commonLcLpiCut+Dst1_D_window+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow+Dst1_D_window,"goff");
    ch1dat -> Draw("mlc>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow+Dst1_D_window,"goff");
    
    Ntot += ch1dat -> Draw("mlc>>+hdat","abs(rmx-2.301)<0.1159"+commonLcLpiCut+Dst2_D_window+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow+Dst2_D_window,"goff");
    ch1dat -> Draw("mlc>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow+Dst2_D_window,"goff");
    
    Ntot += ch1dat -> Draw("mlc>>+hdat","abs(rmx-2.301)<0.1159"+commonLcLpiCut+Dst3_D_window+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow+Dst3_D_window,"goff");
    ch1dat -> Draw("mlc>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow+Dst3_D_window,"goff");

    Ntot += ch1dat -> Draw("mlc>>+hdat","abs(rmx-2.301)<0.1159"+commonLcLpiCut+Dst4_D_window+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow+Dst4_D_window,"goff");
     ch1dat -> Draw("mlc>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow+Dst4_D_window,"goff");
    
    /*TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.287)",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.287)",lend,rend);    
    

    fdat -> SetParameters(20,MLambdac,0.007,10,-300); //100,MLambdac,0.01,40,-300
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,4,covBackground);
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
    */
  
    
  
    hdat -> GetXaxis()-> SetTitle("M(pK^{-}#pi^{+}) [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.0f MeV ",binw*1000));
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
    
    hrsb -> SetLineColor(4);
    hrsb -> SetLineWidth(4);
    //hrsb -> Draw("same");
    
    hlsb -> SetLineColor(8);
    hlsb -> SetLineWidth(4);
    hlsb -> Draw("same");
    
    
    //fbkg -> SetLineStyle(2);
    //fbkg -> SetLineColor(12);
    //fbkg -> SetLineWidth(4);
    //fbkg -> DrawCopy("same");
    
    
    //fsig -> SetLineColor(4);
    //fsig -> SetLineWidth(4);
    //fsig -> Draw("same");
    
    
    //fdat -> SetLineColor(2);
    //fdat -> SetLineWidth(4);
    //fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","e p");
   // leg->AddEntry("fsig","Signal","l");
    //leg -> AddEntry("fdat","fit","l");
    //leg -> AddEntry("fbkg","background","l");
    //leg -> AddEntry("hrsb","M_{recoil}(pD^{*+}#pi^{-}) right sb","l");
    leg -> AddEntry("hlsb","M_{recoil}(pD^{*+}#pi^{-}) left sb","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    //tstatfit -> DrawLatex(0.67, 0.65, Form("N_tot = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
   // tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   
    
} 
