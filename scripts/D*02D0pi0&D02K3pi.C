double MDst0=2.00685, lend=1.998, rend=2.015;
int Nbins=170;
double hwidth = rend-lend, binw = hwidth/Nbins;


Double_t mybkg(Double_t* x, Double_t *par)
{
   Double_t xx = x[0], d = xx-par[0];
   if (xx>par[0])
       return sqrt(xx-par[0])*(par[1]+par[2]*d)+par[3];
   else return par[3];
}

Double_t mysig(Double_t* x, Double_t *par)
{
   Double_t xx = x[0];
   return binw*(par[0]*TMath::Gaus(xx,par[1],par[2],true)+par[3]*TMath::Gaus(xx,par[4],par[5],true));
}

Double_t myfdat(Double_t* x, Double_t *par)
{
   Double_t xx = x[0], d = xx-par[6];
   if (xx>par[6])
       return binw*(par[0]*TMath::Gaus(xx,par[1],par[2],true)+par[3]*TMath::Gaus(xx,par[4],par[5],true))+sqrt(xx-par[6])*(par[7]+par[8]*d)+par[9];
   else return par[9];
}



double fx(double x, double a, double b, double c) { return c-b*pow(x,5./2)-a*pow(x,3./2);} // вычисляемая функция
double dfx(double x, double a, double b, double c) { return -b*5/2*pow(x,3./2)-a*3/2*sqrt(x);} // производная функции

double solve(double x0,double a, double b, double c) {
  double x1  = x0 - fx(x0,a,b,c)/dfx(x0,a,b,c); // первое приближение
  while (fabs(x1-x0)>0.000001) { // пока не достигнута точность 0.000001
    x0 = x1;
    x1 = x0 - fx(x0,a,b,c)/dfx(x0,a,b,c); // последующие приближения
  }
  return x1;
}



void FitDst()
{
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetErrorX(0);
    double axisFontSize = 0.065;
    
    
    
    
    TString datapath = "../analysis/hmerge/";    
    TChain* ch1dat = new TChain("h1"); //
    ch1dat -> Add(datapath+"*.root");
    
    
    
    TCanvas *c1 = new TCanvas("c1","D^{*0} invariant mass",1600,900);
    TH1D* hdat = new TH1D("hdat","D^{*0} #rightarrow D^{0}#pi^{0}, D^{0} #rightarrow K^{-}#pi^{-}#pi^{+}#pi^{+}",Nbins,lend,rend);
    TH1D* hsb = new TH1D("hsb","D^{*0} #rightarrow D^{0}#pi^{0}, D^{0} #rightarrow K^{-}#pi^{-}#pi^{+}#pi^{+}",Nbins,lend,rend);
    
    
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    TCut Mwindow = Form("tag==4 && dstch==1 && dch==2 && mdst > %lf && mdst < %lf && abs(md-1.86483)<0.015",lend,rend);
    Ntot = ch1dat -> Draw("mdst>>hdat",Mwindow,"goff"); //-md+1.86483
    //ch1dat -> Draw("mdst-md+2.01026>>hsb",Mwindow,"goff"); //abs(rmx-2.2969)>0.0468*3 && abs(rmx-2.2969)<0.0468*5
    
    
    
    
    TF1* fsig = new TF1("fsig",mysig,lend,rend,6);
    TF1* fbkg = new TF1("fbkg",mybkg,lend,rend,4);
    TF1* fdat = new TF1("fdat",myfdat,lend,rend,10);
    double initpar[] = {10000,MDst0,0.0005,10000,MDst0,0.001,1.98,1675,25000,2};//-1500000};
    fdat -> SetParameters(initpar); 
    fdat -> SetParLimits(1,MDst0-0.002,MDst0+0.002);
    fdat -> SetParLimits(2,0.0001,0.001);
    fdat -> SetParLimits(4,MDst0-0.002,MDst0+0.002);
    fdat -> SetParLimits(5,0.001,0.002);
    fdat -> SetParLimits(0,0,1e6);
    fdat -> SetParLimits(6,1.998,2.003);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff"); //L S M N Q
    TMatrixDSym covFit = fitResult->GetCovarianceMatrix();
    TMatrixDSym covSignal, covBackground;
    covFit.GetSub(0,5,covSignal);
    covFit.GetSub(6,9,covBackground);
    double * par;
    par = fdat -> GetParameters();
  
   
    fsig -> SetParameters(par);
    fbkg -> SetParameters(par+6);
    
    Nsig = fsig -> Integral(lend,rend)/binw; //fdat -> GetParameter(0); fsig -> Integral(lend,rend)/binw;  
    dNsig = fsig-> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw; //fdat -> GetParError(0); fsig-> IntegralError(lend,rend,par,covSignal.GetMatrixArray())/binw;
    //Nbkg = fbkg -> Integral(lend,rend)/binw; 
    //dNbkg = fbkg -> IntegralError(lend,rend,par+3,covBackground.GetMatrixArray())/binw;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N D*ch: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    //cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;
 
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
  
    
  
    hdat -> GetXaxis()-> SetTitle("M(D^{0}#pi^{0})[GeV]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events /  %.1f MeV ",binw*1000));
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
    
    //hsb -> SetLineColor(8);
   // hsb -> SetLineWidth(4);
    //hsb -> Draw("same");
    
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(4);
    fbkg -> DrawCopy("same");
    
    
    fsig -> SetLineColor(4);
    fsig -> SetLineWidth(4);
    //fsig -> Draw("same");
    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(4);
    fdat-> DrawCopy("same");
    
    
    TLegend* leg = new TLegend(0.77,0.75,0.89,0.95);
    leg -> AddEntry("hdat","data","E P");
   // leg->AddEntry("fsig","Signal","l");
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
    tstatfit -> DrawLatex(0.67, 0.65, Form("N_{sig} = %0.lf #pm %0.lf",Nsig, dNsig)); //
   // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
    tstatfit -> DrawLatex(0.67, 0.59, Form("Mean_{sig1} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
    tstatfit -> DrawLatex(0.67, 0.53, Form("#sigma_{sig2} = %0.4lf #pm %0.4lf",par[5], fdat -> GetParError(5)));

    double widthNbkg3s = fbkg -> Integral(par[1]-3*par[5],par[1]+3*par[5]);
    double C = widthNbkg3s + par[7]*2/3*pow(par[1]+4*par[5]-par[6],3./2)+par[8]*2/5*pow(par[1]+4*par[5]-par[6],5./2), A = par[7]*2/3, B = par[8]*2/5;
    cout << "Useful width = " << solve(par[1]+10*par[5]-par[6],A,B,C)-par[1]-4*par[5]+par[6] << endl;
} 
