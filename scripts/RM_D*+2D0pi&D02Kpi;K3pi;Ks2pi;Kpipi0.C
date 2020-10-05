double fracsigm(double B1, double dB1, double B2, double dB2)
{
    double A, dA;
    A = B1/B2;
    dA = A*sqrt(dB1*dB1/(B1*B1)+dB2*dB2/(B2*B2));
    TRandom r; 
    r.SetSeed(0);
    double buff;
    TH1D H("H", "BrFrac", 100,A-7*dA,A+7*dA);
    
    for(int i =0; i<100000; i++)
    {
        buff = r.Gaus(B1,dB1)/r.Gaus(B2,dB2);
        H.Fill(buff);
    }
    TFitResultPtr fitResult;
    fitResult = H.Fit("gausn","L S M N Q","goff"); //L S M N Q
    return fitResult->Parameter(2);
}








void RM()
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
    
    double lend=1.8, rend=2.8, MLambdac=2.28646, 
        MDst_p=2.01026; //lend=2., rend=2.6
    int Nbins=50;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1024,768);
    TH1D* hdat = new TH1D("hdat","D^{*+} #rightarrow D^{0}#pi^{+}, D^{0}#rightarrow K^{-}#pi^{+},K^{-}#pi^{-}#pi^{+}#pi^{+},K_{S}#pi^{+}#pi^{-},K^{-}#pi^{+}#pi^{0}",Nbins,lend,rend);
    TH1D* hrsb = new TH1D("hrsb","temp",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    double D1_mean = 1.8651, D1_sigma = 0.0064, Dst1_mean = 2.0100, Dst1_sigma = 0.0011, rsbwidth1 = 0.00445;
    TCut Dst1_D_window = Form("dch==1 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D1_mean,3*D1_sigma,Dst1_mean,3*Dst1_sigma);
    TCut Dst1_rightsb = Form("dch==1 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D1_mean,3*D1_sigma,Dst1_mean + 4*Dst1_sigma, Dst1_mean + 4*Dst1_sigma + rsbwidth1);
    
    double D2_mean = 1.8651, D2_sigma = 0.0051, Dst2_mean = 2.010, Dst2_sigma = 0.0011, rsbwidth2 = 0.0048;
    TCut Dst2_D_window = Form("dch==2 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D2_mean,3*D2_sigma,Dst2_mean,3*Dst2_sigma);
    TCut Dst2_rightsb = Form("dch==2 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D2_mean,3*D2_sigma,Dst2_mean + 4*Dst2_sigma, Dst2_mean + 4*Dst2_sigma + rsbwidth2);
    
    double D3_mean = 1.8652, D3_sigma = 0.0055, Dst3_mean = 2.010, Dst3_sigma = 0.0012, rsbwidth3 = 0.00542;
    TCut Dst3_D_window = Form("dch==3 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D3_mean,3*D3_sigma,Dst3_mean,3*Dst3_sigma);
    TCut Dst3_rightsb = Form("dch==3 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D3_mean,3*D3_sigma,Dst3_mean + 4*Dst3_sigma, Dst3_mean + 4*Dst3_sigma + rsbwidth3);
    
    double D4_mean = 1.8639, D4_sigma = 0.0098, Dst4_mean = 2.0100, Dst4_sigma = 0.0015, rsbwidth4 = 0.00698;
    TCut Dst4_D_window = Form("dch==4 && abs(md-%lf)<%lf && abs(mdst-%lf)<%lf",D4_mean,3*D4_sigma,Dst4_mean,3*Dst4_sigma);
    TCut Dst4_rightsb = Form("dch==4 && abs(md-%lf)<%lf && mdst>%lf && mdst<%lf",D4_mean,3*D4_sigma,Dst4_mean + 4*Dst4_sigma, Dst4_mean + 4*Dst4_sigma + rsbwidth4);

    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("rmx > %lf && rmx < %lf",lend,rend);
    Ntot +=  ch1dat -> Draw("rmx>>+hdat","tag==3 && dstch==1"+Dst1_D_window+Mwindow,"goff");
    ch1dat -> Draw("rmx>>+hrsb","tag==3 && dstch==1"+Dst1_rightsb+Mwindow,"goff");
    
    Ntot += ch1dat -> Draw("rmx>>+hdat","tag==3 && dstch==1"+Dst2_D_window+Mwindow,"goff");
    ch1dat -> Draw("rmx>>+hrsb","tag==3 && dstch==1"+Dst2_rightsb+Mwindow,"goff");
    
    Ntot += ch1dat -> Draw("rmx>>+hdat","tag==3 && dstch==1"+Dst3_D_window+Mwindow,"goff");
    ch1dat -> Draw("rmx>>+hrsb","tag==3 && dstch==1"+Dst3_rightsb+Mwindow,"goff");

    Ntot += ch1dat -> Draw("rmx>>+hdat","tag==3 && dstch==1"+Dst4_D_window+Mwindow,"goff");
    ch1dat -> Draw("rmx>>+hrsb","tag==3 && dstch==1"+Dst4_rightsb+Mwindow,"goff");

    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.287)+[5]*(x-2.287)^2",binw),lend,2.5);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,2.5);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.287)+[2]*(x-2.287)^2",lend,2.5);    
    

    fdat -> SetParameters(4000,MLambdac,0.05,150,300,20);
    fdat -> SetParLimits(1,MLambdac-0.1,MLambdac+0.1);
    fdat -> SetParLimits(0,0,1e5);
    fdat -> SetParLimits(2,0.01,0.3);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff",lend,2.5); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,2,covSignal);
    covFit.GetSub(3,5,covBackground);
    double * par;
    par = fdat -> GetParameters();
  
   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+3);
    double mean = fdat -> GetParameter(1), sigma = fdat -> GetParameter(2);
    Nsig = fdat -> GetParameter(0);//fsig -> Integral(lend,rend)/binw; 
    dNsig = fdat -> GetParError(0);// -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg3s = fbkg -> Integral(mean-3*sigma,mean+3*sigma)/binw; 
    dNbkg3s = fbkg -> IntegralError(mean-3*sigma,mean+3*sigma,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg3s << " +- " << (int) dNbkg3s << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = (2.5-lend)/binw - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("M_{recoil}(pD^{*+}#pi^{-}) [GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events / %.2f GeV",binw));
    hdat -> GetYaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetYaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitleOffset(0.9);
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
    
   
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");   
    
    /*
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
    fsig -> Draw("same");
    */
    
    TLegend* leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry("hdat","data","ep");
	//leg->AddEntry("fsig","Signal","l");    
	leg->AddEntry("fdat","fit","l");
    leg->AddEntry("fbkg","background","l");
    leg->AddEntry("hrsb","D* right sideband","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg->Draw("same"); 
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.45, Form("S = %0.lf #pm %0.lf",Nsig, dNsig)); //
    tstatfit -> DrawLatex(0.67, 0.39, Form("B(3#sigma) = %0.lf #pm %0.lf",Nbkg3s, dNbkg3s));
    tstatfit -> DrawLatex(0.67, 0.33, Form("S/B(3#sigma) = %0.2lf #pm %0.2lf",Nsig/Nbkg3s, fracsigm(Nsig,dNsig,Nbkg3s,dNbkg3s)));
    
    // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
   // tstatfit -> DrawLatex(0.67, 0.39, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
   // tstatfit -> DrawLatex(0.67, 0.33, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   
    
}  
