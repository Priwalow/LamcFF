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
    //TString datapath = "../mc_analysis/hmerge/";
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add(datapath+"*.root");
    
    double lend=1.6, rend=2.55, MLambdac=2.28646, 
        MDst_p=2.01026, RightFitLimit = 2.5, LeftFitLimit=lend; //lend=2., rend=2.6
    int Nbins=95;
    TH1D* hdat = new TH1D("hdat","RM(D^{*0})|| RM(D^{*+})",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    


    
    double Ntot=0, Nsig, dNsig, Nbkg3s, dNbkg3s;
    TCut Mwindow = Form("rmx > %lf && rmx < %lf",lend,rend);
    Ntot +=  ch1dat -> Draw("rmx>>+hdat","((tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))))) || (tag==1  && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==2 &&  abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))) )","goff");

    // full (tag==1  && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==2 &&  abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))) || (tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))))
    
    //Only D0 tag==1  && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))
    
    //Only D+ tag==2 &&  abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3))
    
    //Only D*0 tag==4 && dstch==1 && abs(mdst-2.00685)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6))
    
    //Only D*+ tag==3 && ((dstch==1 && abs(mdst-2.01026)<0.002 && abs(md-1.86483)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=3 && dch!=6)))) || (dstch==2 && abs(mdst-2.01026)<0.002 && abs(md-1.86965)<0.015 && (abs(mks-0.497611)<0.0075 || (dch!=2 && dch!=3)))
    
    
    TF1* fdat = new TF1("fdat",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.287)+[5]*(x-2.287)^2",binw),lend,RightFitLimit);
    TF1* fsig = new TF1("fsig",Form("%lf*[0]*TMath::Gaus(x,[1],[2],true)",binw),lend,RightFitLimit);
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.287)+[2]*(x-2.287)^2",lend,RightFitLimit);    
    

    fdat -> SetParameters(4000,MLambdac,0.05,5250,5350,2200);
    fdat -> SetParLimits(1,MLambdac-0.05,MLambdac+0.05);
    fdat -> SetParLimits(0,0,1e5);
    fdat -> SetParLimits(2,0.01,0.3);

    TFitResultPtr fitResult;
    fitResult = hdat -> Fit("fdat","L S M N","goff",lend,RightFitLimit); //L S M N Q
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
    int ndf = (RightFitLimit-lend)/binw - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;
    
    TCanvas *c1 = new TCanvas("c1","Lambda_c invariant mass",1600,900);
    TPad *pad1 = new TPad("pad1","This is pad1",0,0.3,1,1);
    TPad *pad2 = new TPad("pad2","This is pad2",0,0,1,0.3);
    pad1->Draw();
    pad2->Draw();
    
    pad1->cd();
    hdat -> GetXaxis()-> SetTitle("M_{recoil}(X) [GeV/c^{2}]");
    hdat -> GetXaxis()-> SetTitleSize(axisFontSize);
    hdat -> GetXaxis()-> SetLabelSize(axisFontSize);
    hdat -> GetYaxis()-> SetTitle(Form("Events / %.2f GeV/c^{2}",binw));
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
	//leg->AddEntry("fsig","signal","l");    
	leg->AddEntry("fdat","gaus+pol2 fit","l");
    leg->AddEntry("fbkg","pol2 background","l");
    leg -> SetBorderSize(0);
    leg -> SetTextSize(axisFontSize);
    leg->Draw("same"); 
    
    TLatex *tstatfit = new TLatex();
    tstatfit -> SetNDC();
    tstatfit -> SetTextColor(1);
    tstatfit -> SetTextSize(axisFontSize);
    tstatfit -> SetTextAngle(0);
    tstatfit -> DrawLatex(0.67, 0.45, Form("S = %0.lf #pm %0.lf",Nsig, dNsig)); //
    tstatfit -> DrawLatex(0.67, 0.33, Form("#chi^{2}/ndf = %.2lf",chisq/ndf));
    tstatfit -> DrawLatex(0.67, 0.57, Form("gaus #sigma = %.1lf #pm %.1lf MeV/c^{2}",1000*par[2],1000*fdat->GetParError(2)));
    tstatfit -> DrawLatex(0.67, 0.33, Form("S/B(3#sigma) = %0.3lf #pm %0.3lf",Nsig/Nbkg3s, fracsigm(Nsig,dNsig,Nbkg3s,dNbkg3s)));
    
    tstatfit -> SetTextColor(8);
    tstatfit -> DrawLatex(0.67, 0.21, "left sideband");
    
    tstatfit -> SetTextColor(4);
    tstatfit -> DrawLatex(0.67, 0.09, "right");
    tstatfit -> DrawLatex(0.67, 0.0, "sideband");
    // tstatfit -> DrawLatex(0.67, 0.39, Form("B(3#sigma) = %0.lf #pm %0.lf",Nbkg3s, dNbkg3s));
    
    
    // tstatfit -> DrawLatex(0.67, 0.59, Form("N_{bkg} = %0.lf #pm %0.lf",Nbkg, dNbkg)); //
   // tstatfit -> DrawLatex(0.67, 0.39, Form("Mean_{sig} = %0.4lf #pm %0.4lf",par[1], fdat -> GetParError(1)));
   // tstatfit -> DrawLatex(0.67, 0.33, Form("#sigma_{sig} = %0.4lf #pm %0.4lf",par[2], fdat -> GetParError(2)));

   
    double top = hdat->GetMaximum();
    TLine* line1 = new TLine(2.29-0.1,0,2.29-0.1,top*1.);
    TLine* line2 = new TLine(2.29+0.1,0,2.29+0.1,top*1.);
    line1 -> SetLineWidth(4);
    line2 -> SetLineWidth(4);
    
    TLine* line3 = new TLine(2.29-0.65,0,2.29-0.65,top*1.);
    TLine* line4 = new TLine(2.29-0.15,0,2.29-0.15,top*1.);
    line3 -> SetLineWidth(4);
    line3 -> SetLineColor(8);
    line4 -> SetLineWidth(4);
    line4 -> SetLineColor(8);
    
    TLine* line5 = new TLine(2.29+0.15,0,2.29+0.15,top*1.);
    TLine* line6 = new TLine(2.29+0.25,0,2.29+0.25,top*1.);
    line5 -> SetLineWidth(4);
    line5 -> SetLineColor(4);
    line6 -> SetLineWidth(4);
    line6 -> SetLineColor(4);
    
    line1 -> Draw("same");
    line2 -> Draw("same");
    line3 -> Draw("same");
    line4 -> Draw("same");
    line5 -> Draw("same");
    line6 -> Draw("same");
    
    pad2->cd();
    TH1D* hdiff = new TH1D("hdiff","difference",Nbins,lend,rend);
    hdiff->SetName("hdiff");         
    hdat -> Sumw2();
    double firstbin = 1+(LeftFitLimit-lend)/binw, lastbin = (RightFitLimit-LeftFitLimit)/binw;
    for ( int i = firstbin; i <= 1+lastbin; i++ ) { hdiff->SetBinContent(i,(hdat->GetBinContent(i)-fdat->Eval(hdat->GetBinCenter(i)))/(hdat->GetBinError(i)));
    } 
    
    hdiff -> GetXaxis()-> SetLabelSize(0.1);
    hdiff -> GetYaxis()-> SetTitle("Pull");
    hdiff -> GetYaxis()-> SetTitleSize(0.1);
    hdiff -> GetYaxis()-> SetLabelSize(0.1);
    hdiff -> GetYaxis()-> SetTitleOffset(0.2);
    //hdat -> GetYaxis()->CenterTitle(true);
    hdiff -> GetXaxis()->SetTickSize(0.04);
    hdiff -> SetMarkerStyle(20);
    hdiff -> SetMarkerSize(1.5);
    hdiff -> SetMarkerColor(1);
    hdiff -> SetLineColor(1);
    hdiff -> SetLineWidth(2);
    hdiff -> SetMaximum(+3);
    hdiff -> SetMinimum(-3);
    hdiff -> Draw("p");
    
    TF1* fzero = new TF1("fzero","0",lend,RightFitLimit);
    fzero -> SetLineColor(2);
    fzero -> SetLineWidth(4);
    fzero-> DrawCopy("same");  
}  
