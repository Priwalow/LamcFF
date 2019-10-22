{
    gStyle->SetOptStat(0);
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add("*.root");
    
    double lend=0, rend=1;
    int Nbins=100;
    TCanvas *c1 = new TCanvas("c1","Muon neutrino invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","#nu_{mu} invariant mass",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    ch1dat -> Draw("rm*rm>>hdat","lcch == 4","goff"); 
    
    hdat -> Scale(1,"width");
    
    /*TF1* fdat = new TF1("fdat","[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*(x-2.3)",lend,rend);
    

    fdat -> SetParameters(600,MLambdac,0.1,3000,-700);
    fdat -> SetParLimits(1,MLambdac-0.2,MLambdac+0.2);

    hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q */

 

    hdat -> GetXaxis()->SetTitle("M(#mu_{nu}), GeV");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.4f GeV",binw));
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> DrawCopy("ep");
   
    
    /*Nsig = fdat->GetParameter(0);
    dNsig = fdat->GetParError(0);
    
    double a = fdat -> GetParameter(3), b = fdat -> GetParameter(4), da = fdat -> GetParError(3), db = fdat -> GetParError(4); 
    Nbkg = (a-b*2.3)*hwidth+b*(rend*rend-lend*lend)/2;
    dNbkg = hwidth*sqrt(da*da+ db*db*((rend+lend)*(rend+lend)/4+2.3*2.3));
    
    //cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;

    
  
    
    TF1* fbkg = new TF1("fbkg","[0]+[1]*(x-2.3)",lend,rend);    
    fbkg -> SetParameters(fdat->GetParameter(3),fdat->GetParameter(4));
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(3);
    fbkg -> DrawCopy("same");
    
    /*
     * TF1* fsig = new TF1("fsig","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)+[6]*TMath::Gaus(x,[7],[8],true)",lend,rend);
    fsig -> SetParameters(fdat -> GetParameter(0), fdat -> GetParameter(1), fdat -> GetParameter(2), fdat -> GetParameter(4), fdat -> GetParameter(5),
                           fdat -> GetParameter(6), fdat -> GetParameter(7), fdat -> GetParameter(8), fdat -> GetParameter(9));
    fsig -> SetLineColor(3);
    fsig -> SetLineWidth(3);
    fsig -> Draw("same");
    
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(3);
    fdat-> DrawCopy("same");
    
    TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
    leg->AddEntry("fdat","Signal + background","l");
    leg->AddEntry("fbkg","Background","l");
    leg->Draw("same"); */
    
} 
 
