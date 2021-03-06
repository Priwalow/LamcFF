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





void Lc2Lpi()
{
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    double axisFontSize = 0.065;
    
    
    
    
    TString datapath = "../analysis/hmerge/";    
    //TString datapath = "../mc_analysis/hmerge/";
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    double lend=2.24, rend=2.33, MLambdac=2.28646; //lend=2.21, rend=2.36
    int Nbins=18;
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}",Nbins,lend,rend);
    TH1D* hrsb = new TH1D("hrsb","right sb",Nbins,lend,rend);
    TH1D* hlsb = new TH1D("hlsb","left sb",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;

    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("mlc > %lf && mlc < %lf",lend,rend);
    TCut commonLcLpiCut = "lcch==1 && abs(ml-1.11568)<0.003 && ((tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))))) || (tag==1  && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==2 &&  abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))) ) && abs(pvis)<0.05 && abs(ecms-sqrt(mvis*mvis+pvis*pvis))<0.05"; //
    TCut rmxlsbCut = "rmx-2.29 < -0.15 && rmx-2.29 > -0.65";
    TCut rmxrsbCut = "rmx-2.29 > 0.15 && rmx-2.29 < 0.25";
    Ntot +=  ch1dat -> Draw("mlc>>+hdat","abs(rmx-2.29)<0.1"+commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hrsb",rmxrsbCut+commonLcLpiCut+Mwindow,"goff");
    ch1dat -> Draw("mlc>>+hlsb",rmxlsbCut+commonLcLpiCut+Mwindow,"goff");
    
    
   TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.287)",binw),lend,rend);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,rend);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.287)",lend,rend);    
    

    fdat -> SetParameters(40,MLambdac,0.007,2,-10); //100,MLambdac,0.01,40,-300
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
    double mean = fdat -> GetParameter(1), sigma = fdat -> GetParameter(2);
    Nsig = fdat -> GetParameter(0);//fsig -> Integral(lend,rend)/binw; 
    dNsig = fdat -> GetParError(0);// -> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    Nbkg3s = fbkg -> Integral(mean-3*sigma,mean+3*sigma)/binw; 
    dNbkg3s = fbkg -> IntegralError(mean-3*sigma,mean+3*sigma,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda_c: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg3s << " +- " << (int) dNbkg3s << endl;

 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    /* */
  
    
  
    hdat -> GetXaxis()-> SetTitle("M(#Lambda#pi^{+}) [GeV/c^{2}]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.0f MeV/c^{2} ",binw*1000));
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
    
    hlsb -> Scale(0.2);
    hlsb -> SetLineColor(8);
    hlsb -> SetLineWidth(4);
    hlsb -> Draw("hist same");
    
    
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
    leg -> AddEntry("fdat","gaus+pol1 fit","l");
    //leg -> AddEntry("fbkg","background","l");
    leg -> AddEntry("hrsb","M_{recoil} right sb","l");
    leg -> AddEntry("hlsb","M_{recoil} left sb","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg -> Draw("same");
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.65, Form("S = %0.lf #pm %0.lf",Nsig, dNsig));
    tstatfit -> DrawLatex(0.67, 0.53, Form("gaus #sigma = %0.1lf #pm %0.1lf MeV/c^{2}",1000*par[2], 1000*fdat->GetParError(2)));//
    tstatfit -> DrawLatex(0.67, 0.59, Form("B(3#sigma) = %.0lf #pm %.0lf",Nbkg3s, dNbkg3s)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("gaus #mu = %.1lf #pm %.1lf MeV/c^{2}",1000*par[1], 1000*fdat -> GetParError(1)));
   // tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));
     tstatfit -> DrawLatex(0.67, 0.33, Form("S/B(3#sigma) = %.1lf #pm %.1lf",Nsig/Nbkg3s, fracsigm(Nsig,dNsig,Nbkg3s,dNbkg3s)));
    
    
    double top = hdat->GetMaximum();
    tstatfit -> DrawLatex(0.67, 0.65, "#Lambda^{+}_{c} #rightarrow #Lambda#pi^{+}"); 
    TLine* line2 = new TLine(MLambdac-0.02,0,MLambdac-0.02,top*0.95);
    TLine* line3 = new TLine(MLambdac+0.02,0,MLambdac+0.02,top*0.95);
    line2 -> SetLineWidth(4);
    line3 -> SetLineWidth(4);
    line2 -> Draw("same");
    line3 -> Draw("same");
    
} 
 
