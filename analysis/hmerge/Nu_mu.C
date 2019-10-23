{
    gStyle->SetOptStat(0);
    TChain* ch1dat = new TChain("h1"); //without pi0
    ch1dat -> Add("*.root");
    
    double lend=-0.5, rend=0.5;
    int Nbins=25;
    TCanvas *c1 = new TCanvas("c1","Muon neutrino invariant mass",740,600);
    TH1D* hdat = new TH1D("hdat","Squared recoil mass of #nu_{#mu}",Nbins,lend,rend);
    double hwidth = rend-lend, binw = hwidth/Nbins;
    
    double Ntot, Nsig, dNsig, Nbkg, dNbkg;
    ch1dat -> Draw("rm*fabs(rm)>>hdat","lcch == 4 && m*m-10.58*10.58+2*10.58*p-rm*rm > -2 && m*m-10.58*10.58+2*10.58*p-rm*rm < 2  && p>0.3","goff"); 
    
    hdat -> Scale(1,"width");
    
    TF1* fdat = new TF1("fdat","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)+[6]*TMath::Gaus(x,[7],[8],true)",lend,rend);
    

    fdat -> SetParameters(58,0.02,0.1,13,-0.005,1.59,6.1,0.2,0.5);
    //fdat -> SetParLimits(1,-0.02,0.02);

    hdat -> Fit("fdat","L S M N Q","goff"); //L S M N Q

 

    hdat -> GetXaxis()->SetTitle("RM^{2}(#nu_{#mu}), GeV");
    hdat -> GetYaxis()->SetTitle(Form("Events / %.2f",binw));
    hdat -> SetMarkerStyle(20);
    hdat -> SetMarkerSize(1);
    hdat -> SetMarkerColor(1);
    hdat -> SetLineColor(1);
    hdat -> SetLineWidth(2);
    hdat -> DrawCopy("ep");
   
    
    Nsig = fdat->GetParameter(0);
    dNsig = fdat->GetParError(0);
    
    Nbkg = fdat -> GetParameter(3)+fdat -> GetParameter(6);
    dNbkg = sqrt(fdat -> GetParError(3)*fdat -> GetParError(3)+fdat -> GetParError(6)*fdat -> GetParError(6));
    
    //cout << "Total number of events: " << Ntot << endl;
    cout << "N nu_mu: " << Nsig <<" +- " <<  dNsig << endl;
    cout << "N background: " << Nbkg << " +- " <<  dNbkg << endl;
    
    double chisq = fdat->GetChisquare(); 
    int ndf = Nbins - fdat ->GetNumberFreeParameters(); // nbins - npar
    cout << "ChiSquare / ndf: " << chisq << " / " << ndf << " = " << chisq/ndf << endl;

    TF1* fbkg = new TF1("fbkg","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)",lend,rend);    
    fbkg -> SetParameters(fdat->GetParameter(3),fdat->GetParameter(4),fdat->GetParameter(5),
			fdat->GetParameter(6), fdat->GetParameter(7), fdat->GetParameter(8));
    fbkg -> SetLineStyle(2);
    fbkg -> SetLineColor(12);
    fbkg -> SetLineWidth(3);
    fbkg -> DrawCopy("same");
    
    /* TF1* fsig = new TF1("fsig","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)+[6]*TMath::Gaus(x,[7],[8],true)",lend,rend);
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
 
