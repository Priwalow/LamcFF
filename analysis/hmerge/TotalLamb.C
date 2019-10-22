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
    
    hdat -> Scale(1,"width");
    
    TF1* fdat = new TF1("fdat","[0]*TMath::Gaus(x,[1],[2],true)+[3]+[4]*TMath::Gaus(x,[5],[6],true)+[7]*TMath::Gaus(x,[8],[9],true)",lend,rend);
    

    fdat -> SetParameters(20000*3,MLambda,0.0001,4000*3,20000*3,MLambda,0.001,200*3,MLambda,0.01);
    fdat -> SetParLimits(1,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(5,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(8,MLambda-0.1,MLambda+0.1);
    fdat -> SetParLimits(0,0,1e7);
    fdat -> SetParLimits(4,0,1e7);
    fdat -> SetParLimits(7,0,1e7);

    hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q

 

    hdat -> GetXaxis()->SetTitle("M(#Lambda), GeV");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.4f GeV",binw));
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(3);
    hdat -> DrawCopy("ep");
   
    
    Nsig = (fdat->GetParameter(0)+fdat->GetParameter(4)+fdat->GetParameter(7));///binw;//*Ntot;
    dNsig = sqrt((fdat->GetParError(0)*fdat->GetParError(0)+fdat->GetParError(4)*fdat->GetParError(4)+fdat->GetParError(7)*fdat->GetParError(7)));///binw;//*Ntot;
    Nbkg = fdat -> GetParameter(3)*hwidth;
    dNbkg = fdat -> GetParError(3)*hwidth;
    
    cout << "Total number of events: " << Ntot << endl;
    cout << "N Lambda: " << (int) Nsig <<" +- " << (int) dNsig << endl;
    cout << "N background: " << (int) Nbkg << " +- " << (int) dNbkg << endl;

    
  
    
    TF1* fbkg = new TF1("fbkg","[0]",lend,rend);    
    fbkg -> SetParameter(0,fdat->GetParameter(3));
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
    */
    
    fdat -> SetLineColor(2);
    fdat -> SetLineWidth(3);
    fdat-> DrawCopy("same");
    
    TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->AddEntry("hdat","Data","lep");
    leg->AddEntry("fdat","Signal + background","l");
    leg->AddEntry("fbkg","Background","l");
    leg->Draw("same");
    
} 
