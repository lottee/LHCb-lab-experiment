#include <"TMath.h">
#include <"TF1.h">


void drawOutput() {
   const int DalitzBins(14); // Number of bins used currently in the dalitz plot
  const Int_t xNBINS = 21;
  const Int_t yNBINS = 24;
  
  //////////////////////
  // Example Root Macro for third year B->hhh Lab
  ////////////////////// 

  // Open the root file that was produced by running the example program
  TFile *f = new TFile("outputDataAll.root");
  //TFile *f = new TFile("outputDataMagnetDown.root");

  //Section 5.5 Histograms
  TH1F *MassBPos = (TH1F*)f->Get("h_MassBPos");
  TH1F *MassBNeg = (TH1F*)f->Get("h_MassBNeg");
  TH2F *Dalitz21_Pos = (TH2F*)f->Get("h_Dalitz21_Pos");
  TH2F *Dalitz23_Pos = (TH2F*)f->Get("h_Dalitz23_Pos");
  TH2F *Dalitz31_Pos = (TH2F*)f->Get("h_Dalitz31_Pos");
  TH2F *Dalitz21_Neg = (TH2F*)f->Get("h_Dalitz21_Neg");
  TH2F *Dalitz23_Neg = (TH2F*)f->Get("h_Dalitz23_Neg");  
  TH2F *Dalitz31_Neg = (TH2F*)f->Get("h_Dalitz31_Neg");
  
  TH1F *TwoBody1_Pos = (TH1F*)f->Get("h_MassM1_Pos");
  TH1F *TwoBody2_Pos = (TH1F*)f->Get("h_MassM2_Pos");
  TH1F *TwoBody3_Pos = (TH1F*)f->Get("h_MassM3_Pos");
  TH1F *TwoBody1_Neg = (TH1F*)f->Get("h_MassM1_Neg");
  TH1F *TwoBody2_Neg = (TH1F*)f->Get("h_MassM2_Neg");
  TH1F *TwoBody3_Neg = (TH1F*)f->Get("h_MassM3_Neg");
  
  TH2F *Dalitz21_Pos_BG = (TH2F*)f->Get("h_Dalitz21_Pos_BG");
  TH2F *Dalitz23_Pos_BG = (TH2F*)f->Get("h_Dalitz23_Pos_BG");
  TH2F *Dalitz31_Pos_BG = (TH2F*)f->Get("h_Dalitz31_Pos_BG");
  TH2F *Dalitz21_Neg_BG = (TH2F*)f->Get("h_Dalitz21_Neg_BG");
  TH2F *Dalitz23_Neg_BG = (TH2F*)f->Get("h_Dalitz23_Neg_BG");
  TH2F *Dalitz31_Neg_BG = (TH2F*)f->Get("h_Dalitz31_Neg_BG");
  
  
  
  // Fitting functions initially to get the parameters, need specific ranges
  
  /// B+ Getting parameters
  TF1 *backExp_Pos = new TF1("backExp_Pos", "expo", 5100, 5600);
  TF1 *backGauss_Pos = new TF1("backGauss_Pos", "gausn", 5010, 5200);
  TF1 *peakGauss_Pos = new TF1("peakGauss_Pos", "gausn", 5220, 5340);
  TF1 *totalFit_Pos = new TF1("totalFit_Pos", "expo(0)+gausn(2)+gausn(5)", 5010, 5600);

  Double_t par_pos[8];
  
  MassBPos->Fit(backExp_Pos,"R"); // no R+ as first
  MassBPos->Fit(backGauss_Pos,"R+");
  MassBPos->Fit(peakGauss_Pos,"R+");
  //MassBPos->Fit(sumFit,"R");
  
    //extract parameters from:
  backExp_Pos->GetParameters(&par_pos[0]);
  backGauss_Pos->GetParameters(&par_pos[2]);
  peakGauss_Pos->GetParameters(&par_pos[5]);
  //totalFit_Pos->SetParameters(par_pos); // set the sum of the fits
  
  totalFit_Pos->SetParameter(3,5.05378e+03);//mean background gauss
  totalFit_Pos->SetParameter(4,1.12873e+02);//width of backfground gasu
  totalFit_Pos->SetParameter(6,5.28216e+03);//mean peak
  totalFit_Pos->SetParameter(7,2.94521e+01);//width peak`
  
  MassBPos->Fit(totalFit_Pos,"R");
  
  Double_t FinPar_Pos[8];
  
  totalFit_Pos->GetParameters(FinPar_Pos);
  
  std::cout<<"Reduced Chi Square B+ ="<<(totalFit_Pos->GetChisquare()) / (totalFit_Pos -> GetNDF())<<std::endl << std::endl;
  
  //sumFit->GetParameters(&par[0]);
    
  // Make new functions to with the previous parameters found
  // Replaced user fitted functions with predefined functipons with parameters found from the last // 14-03-18
  TF1 *ExpBG_Pos = new TF1("ExpBG_Pos", "expo", 5010, 5600);
  TF1 *FourBodyBG_Pos = new TF1("FourBodyBG_Pos", "gausn", 5010, 5600);
  TF1 *BPosPeak = new TF1("BPosPeak_Pos", "gausn", 5010, 5600);
  
  /// B- Getting parameters
  TF1 *backExp_Neg = new TF1("backExp_Neg", "expo", 5200, 5600);
  TF1 *backGauss_Neg = new TF1("backGauss_Neg", "gausn", 5010, 5200);
  TF1 *peakGauss_Neg = new TF1("peakGauss_Neg", "gausn", 5220, 5360);
  TF1 *totalFit_Neg = new TF1("totalFit_Neg", "expo(0)+gausn(2)+gausn(5)", 5010, 5600);

  Double_t par_neg[8];
  
  MassBNeg->Fit(backExp_Neg,"R"); // no R+ as first
  MassBNeg->Fit(backGauss_Neg,"R+");
  MassBNeg->Fit(peakGauss_Neg,"R+");
  //MassBPos->Fit(sumFit,"R");
  
    //extract parameters from:
  backExp_Neg->GetParameters(&par_neg[0]);
  backGauss_Neg->GetParameters(&par_neg[2]);
  peakGauss_Neg->GetParameters(&par_neg[5]);
  //totalFit_Neg->SetParameters(&par_neg[0],&par_neg[1],&par_neg[2],&par_neg[3],&par_neg[4],&par_neg[5],&par_neg[6],&par_neg[7]); // set the sum of the fits
  
  totalFit_Neg->SetParameter(3,5.03811e+03);//mean background gauss
  totalFit_Neg->SetParameter(4,1.16857e+02);//width of backfground gasu
  totalFit_Neg->SetParameter(6,5.28371e+03);//mean peak
  totalFit_Neg->SetParameter(7,3.10289e+01);//width peak`
  
  
  MassBNeg->Fit(totalFit_Neg,"R");
  
  Double_t FinPar_neg[8];
  
  totalFit_Neg->GetParameters(FinPar_neg);
  
  std::cout<<"Reduced Chi Square B- ="<<(totalFit_Neg->GetChisquare()) / (totalFit_Neg -> GetNDF())<<std::endl << std::endl;
  
  // Make new functions to with the previous parameters found
  // Replaced user fitted functions with predefined functipons with parameters found from the last // 14-03-18
  TF1 *ExpBG_Neg = new TF1("ExpBG_Neg", "expo", 5010, 5600);
  TF1 *FourBodyBG_Neg = new TF1("FourBodyBG_Neg", "gausn", 5010, 5600);
  TF1 *BNegPeak = new TF1("BNegPeak_Neg", "gausn", 5010, 5600);
  
//  ExpBG_Neg->SetParameters(&FinPar_neg[0]);
//  FourBodyBG_Neg->SetParameters(&FinPar_neg[2]);
//  BNegPeak->SetParameters(&FinPar_neg[5]);
  
 
  /// Calculating Global Asymmetry
  const int nbins = 300;
  double bin_range = 560; 
  double bin_width = bin_range / nbins;
  // Calculating N+
  Double_t N_Pos = par_pos[5] / bin_width;
  
  // Calculating N-
  Double_t N_Neg = par_neg[5] / bin_width;
  
  // Calculating A
  Double_t Asym = (N_Neg - N_Pos)/(N_Neg + N_Pos);
  Double_t A_error = sqrt((1-Asym*Asym)/(N_Neg+N_Pos));
  std::cout<<"N_Neg=" << N_Neg <<" N_Pos="<< N_Pos<<" Asym="<<Asym<<"+/-"<<A_error<< std::endl;
  

  
  
  /// B+ canvas
  TCanvas *c9 = new TCanvas("c9","B+",600,400);
  c9->SetGrid();
  
  for(int bin(1); bin <= nbins; bin++){
  // Adding bin errors
  double binContent(0);
  binContent = MassBPos->GetBinContent(bin);
  MassBPos->SetBinError(bin,TMath::Sqrt(binContent));
  MassBPos->SetMarkerStyle(20);
  MassBPos->SetMarkerSize(0.4);
  }
  
  MassBPos->Draw(); // Draw the data onto the graph
  MassBPos->SetStats(0);
  MassBPos->GetXaxis()->SetTitle("Invariant mass [MeV/c^2]");
  MassBPos->GetYaxis()->SetTitle("Events");
  MassBPos->SetTitle("Global Asymmetry B^{+} #rightarrow K^{+} #pi^{-} #pi^{+}");
  ExpBG_Pos->SetParameters(&FinPar_Pos[0]);
  FourBodyBG_Pos->SetParameters(&FinPar_Pos[2]);
  BPosPeak->SetParameters(&FinPar_Pos[5]);
  
  // Formatting the fits
  ExpBG_Pos->SetLineColor(8);
  ExpBG_Pos->SetLineWidth(2);
  ExpBG_Pos->SetLineStyle(7);
  

  FourBodyBG_Pos->SetLineColor(9);
  FourBodyBG_Pos->SetLineWidth(2);
  FourBodyBG_Pos->SetLineStyle(2);
  FourBodyBG_Pos->SetFillStyle(3004);
  FourBodyBG_Pos->SetFillColor(9);
  
  BPosPeak->SetLineColor(6);
  BPosPeak->SetFillColorAlpha(6,0.9);
  BPosPeak->SetFillStyle(3003);
  
  totalFit_Pos->SetLineColor(2);
  totalFit_Pos->SetLineWidth(2);
  
  //Draw
  ExpBG_Pos->Draw("same");
  FourBodyBG_Pos->Draw("same");
  BPosPeak->Draw("same");
  totalFit_Pos->Draw("same");

  
  
 //Draw Legend:
  TLegend *legend1 = new TLegend(0.6,0.65,0.88,0.85);
  legend1->SetTextFont(72);
  legend1->SetTextSize(0.04);
  legend1->SetFillStyle(0);
  legend1->AddEntry(MassBPos,"Data","lpe");
  legend1->AddEntry(totalFit_Pos,"Total", "L");
  legend1->AddEntry(FourBodyBG_Pos,"B^{+} #rightarrow 4 Body", "L");
  legend1->AddEntry(BPosPeak,"B^{+} #rightarrow K^{+} #pi^{-} #pi^{+}", "L");
  legend1->AddEntry(ExpBG_Pos,"Background", "L");
  legend1->Draw();

  c9->SaveAs("InvMassBmeson_Pos.pdf");
  
  /// B- CANVAS
  
  TCanvas *c10 = new TCanvas("c10","B-",600,400);
  c10->SetGrid();
  
  
  
  for(int bin(1); bin <= nbins; bin++){
  // Adding bin errors
  double binContent(0);
  binContent = MassBNeg->GetBinContent(bin);
  MassBNeg->SetBinError(bin,TMath::Sqrt(binContent));
  MassBNeg->SetMarkerStyle(20);
  MassBNeg->SetMarkerSize(0.4);
  
  }
  
  
  MassBNeg->Draw(); // Draw the data onto the graph
  MassBNeg->SetStats(0);
  MassBNeg->GetXaxis()->SetTitle("Invariant mass [MeV/c^2]");
  MassBNeg->GetYaxis()->SetTitle("Events");

  MassBNeg->SetTitle("Global Asymmetry B^{-} #rightarrow K^{-} #pi^{-} #pi^{+}");
  ExpBG_Neg->SetParameters(&FinPar_neg[0]);
  FourBodyBG_Neg->SetParameters(&FinPar_neg[2]);
  BNegPeak->SetParameters(&FinPar_neg[5]);
  //Total_Neg->SetParameters(par_neg);
  
  // Formatting the fits
  ExpBG_Neg->SetLineColor(8);
  ExpBG_Neg->SetLineWidth(2);
  ExpBG_Neg->SetLineStyle(7);
  
  FourBodyBG_Neg->SetLineColor(9);
  FourBodyBG_Neg->SetLineWidth(2);
  FourBodyBG_Neg->SetLineStyle(2);
  FourBodyBG_Neg->SetFillStyle(3004);
  FourBodyBG_Neg->SetFillColor(9);
  
  BNegPeak->SetLineColor(6);
  BNegPeak->SetFillColorAlpha(6,0.9);
  BNegPeak->SetFillStyle(3003);
  
  totalFit_Neg->SetLineColor(2);
  totalFit_Neg->SetLineWidth(2);
  
  //Draw
  ExpBG_Neg->Draw("same");
  FourBodyBG_Neg->Draw("same");
  BNegPeak->Draw("same");
  totalFit_Neg->Draw("same");
  
  legend1->Draw();

 //Draw Legend:
 
  TLegend *legend2 = new TLegend(0.6,0.65,0.88,0.85);
  legend2->SetTextFont(72);
  legend2->SetTextSize(0.04);
  legend2->SetFillStyle(0);
  legend2->AddEntry(MassBNeg,"Data","lpe");
  legend2->AddEntry(totalFit_Neg,"Total", "L");
  legend2->AddEntry(FourBodyBG_Neg,"B^{-} #rightarrow 4 Body", "L");
  legend2->AddEntry(BNegPeak,"B^{-} #rightarrow K^{-} #pi^{-} #pi^{+}", "L");
  legend2->AddEntry(ExpBG_Neg,"Background", "L");
  legend2->Draw();
  
  c10->SaveAs("InvMassBmeson_Neg.pdf");
  
  /// Merging the two canvas' for B- and B+
  TCanvas *c2 = new TCanvas("c2","multipads",900,700);
  
  gStyle->SetOptStat(0);
  c2->Divide(2,0);

  c2->cd(1);
  c2->cd(1)->SetGrid();
  gPad->SetTicky(2);
  
  // Formatting the fits
  ExpBG_Pos->SetLineColor(8);
  ExpBG_Pos->SetLineWidth(2);
  ExpBG_Pos->SetLineStyle(7);
  
  FourBodyBG_Pos->SetLineColor(9);
  FourBodyBG_Pos->SetLineWidth(2);
  FourBodyBG_Pos->SetLineStyle(2);
  FourBodyBG_Pos->SetFillStyle(3004);
  FourBodyBG_Pos->SetFillColor(9);
  
  BPosPeak->SetLineColor(6);
  BPosPeak->SetFillColorAlpha(6,0.9);
  BPosPeak->SetFillStyle(3003);
  
  totalFit_Pos->SetLineColor(2);
  totalFit_Pos->SetLineWidth(2);
  
  MassBPos->Draw();
  ExpBG_Pos->Draw("same");
  FourBodyBG_Pos->Draw("same");
  BPosPeak->Draw("same");
  totalFit_Pos->Draw("same");


  c2->cd(2);
  c2->cd(2)->SetGrid();
  gPad->SetTicky(2);
  
  // Formatting the fits
  ExpBG_Neg->SetLineColor(8);
  ExpBG_Neg->SetLineWidth(2);
  ExpBG_Neg->SetLineStyle(7);
  
  FourBodyBG_Neg->SetLineColor(9);
  FourBodyBG_Neg->SetLineWidth(2);
  FourBodyBG_Neg->SetLineStyle(2);
  FourBodyBG_Neg->SetFillStyle(3004);
  FourBodyBG_Neg->SetFillColor(9);
  
  BNegPeak->SetLineColor(6);
  BNegPeak->SetFillColorAlpha(6,0.9);
  BNegPeak->SetFillStyle(3003);
  
  totalFit_Neg->SetLineColor(2);
  totalFit_Neg->SetLineWidth(2);
  
  MassBNeg->GetYaxis()->SetLabelOffset(0.01);
  MassBNeg->Draw();
  ExpBG_Neg->Draw("same");
  FourBodyBG_Neg->Draw("same");
  BNegPeak->Draw("same");
  totalFit_Neg->Draw("same");

  c2->SaveAs("InvMassBmeson_Multipad.pdf");

  /// Chi-sqr for b+ and B-
  //std::cout<<"chisqr_B+ ="<<totalFit_Pos->GetChiSquare()<<std::endl;
  //std::cout<<"chisqr_B- ="<<totalFit_Neg->GetChiSquare()<<std::endl;
/// New canvas for the Dalitz plot
// Positve Dalitza:

  TCanvas *c1 = new TCanvas("c11","B+ K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Pos->SetTitle("B+ K#Pi VS #Pi#Pi");
  Dalitz23_Pos->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Pos->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Pos->SetStats(0);
  Dalitz23_Pos->Draw("colz"); // makes the plot colour scaled
  c11->SaveAs("B+_Signal_Dalitz.pdf");
/*
  TCanvas *c2 = new TCanvas("c12","B+ #Pi#Pi VS K#Pi",600,400);
  Dalitz31_Pos->SetTitle("B+ #Pi#Pi VS K#Pi");
  Dalitz31_Pos->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos->SetStats(0);
  Dalitz31_Pos->Draw("colz"); // makes the plot colour scaled
  c12->SaveAs("B+ KPi Vs PiPi.pdf");

  TCanvas *c3 = new TCanvas("c13","B+ K#Pi#Pi",600,400);
  Dalitz21_Pos->SetTitle("B+ K#Pi VS K#Pi");
  Dalitz21_Pos->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Pos->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Pos->SetStats(0);
  Dalitz21_Pos->Draw("colz"); // makes the plot colour scaled
  c13->SaveAs("B+ KPi Vs PiK.pdf");
*/

  // Negative Dalitz:
  TCanvas *c4 = new TCanvas("c17","B- K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Neg->SetTitle("B- K#Pi VS #Pi#Pi");
  Dalitz23_Neg->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Neg->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Neg->SetStats(0);
  Dalitz23_Neg->Draw("colz"); // makes the plot colour scaled
  c17->SaveAs("B-_Signal_Dalitz.pdf");
/*
  TCanvas *c5 = new TCanvas("c18","B- #Pi#Pi VS K#Pi",600,400);
  Dalitz31_Neg->SetTitle("B- #Pi#Pi VS K#Pi");
  Dalitz31_Neg->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg->SetStats(0);
  Dalitz31_Neg->Draw("colz"); // makes the plot colour scaled
  c18->SaveAs("B-KPi Vs PiPi.pdf");

  TCanvas *c6 = new TCanvas("c19","B- K#Pi#Pi",600,400);
  Dalitz21_Neg->SetTitle("B- K#Pi VS K#Pi");
  Dalitz21_Neg->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Neg->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Neg->SetStats(0);
  Dalitz21_Neg->Draw("colz"); // makes the plot colour scaled
  c19->SaveAs("B- KPi Vs PiK.pdf");
  
  TCanvas *c14 = new TCanvas("c14","Pi+ K+",600,400);
  TwoBody2_Pos->SetTitle("M_{#Pi+K+}");
  TwoBody2_Pos->GetXaxis()->SetTitle("M_{#Pi+K+}");
  TwoBody2_Pos->SetStats(0);
  TwoBody2_Pos->SetLineColor(kRed);
  TwoBody2_Pos->Draw(); 
  c14->SaveAs("TwoBody PI_K_Pos.pdf");
  
  TCanvas *c15 = new TCanvas("c15","Pi-K+",600,400);
  TwoBody1_Pos->SetTitle("M_{#Pi-K+}");
  TwoBody1_Pos->GetXaxis()->SetTitle("M_{#Pi-K+}");
  TwoBody1_Pos->SetStats(0);
  TwoBody1_Pos->SetLineColor(kRed);
  TwoBody1_Pos->Draw(); 
  c15->SaveAs("TwoBody K_PI_Pos.pdf");
*/
  TCanvas *c16 = new TCanvas("c16","Pi- K+",600,400);
  TwoBody3_Pos->SetTitle("M_{#Pi-K+}");
  TwoBody3_Pos->GetXaxis()->SetTitle("M_{#Pi-K+}");
  TwoBody3_Pos->SetStats(0);
  TwoBody3_Pos->SetLineColor(kRed);
  TwoBody3_Pos->Draw(); 
  c16->SaveAs("TwoBody PI_PI_Pos.pdf");
  /*
  TCanvas *c7 = new TCanvas("c7","Pi+ K-",600,400);
  TwoBody2_Neg->SetTitle("M_{#Pi-K-}");
  TwoBody2_Neg->GetXaxis()->SetTitle("M_{#Pi-K-}");
  TwoBody2_Neg->SetStats(0);
  TwoBody2_Neg->SetLineColor(kRed);
  TwoBody2_Neg->Draw(); 
  c7->SaveAs("TwoBody PI_K_Neg.pdf");
  
  TCanvas *c8 = new TCanvas("c8","Pi+ K-",600,400);
  TwoBody1_Neg->SetTitle("M_{#Pi+K-}");
  TwoBody1_Neg->GetXaxis()->SetTitle("M_{#Pi+K-}");
  TwoBody1_Neg->SetStats(0);
  TwoBody1_Neg->SetLineColor(kRed);
  TwoBody1_Neg->Draw(); 
  c8->SaveAs("TwoBody K_PI_Negs.pdf");
*/
  TCanvas *c9 = new TCanvas("c9","Pi^{-} Pi^{+} 2-Body Resonance Invariant mass",600,400);
  TwoBody3_Neg->SetTitle("M_{#Pi+#Pi-}");
  TwoBody3_Neg->GetXaxis()->SetTitle("M_{#Pi+#Pi-}");
  TwoBody3_Neg->SetStats(0);
  TwoBody3_Neg->SetLineColor(kRed);
  TwoBody3_Neg->Draw(); 
  c9->SaveAs("TwoBody PI_PI_Neg.pdf");
 
 /// Merging the two canvas' for B- and B+
/* TCanvas *c2 = new TCanvas("c2","multipads",900,700);
  
  gStyle->SetOptStat(0);
  c2->Divide(2,0);

  c2->cd(1);
  c2->cd(1)->SetGrid();
  gPad->SetTicky(2);
  TwoBody2_Neg->SetTitle("M_{#Pi^{+} K^{-}} 2-Body Resonance Invariant mass Negative (no cuts)");
  TwoBody2_Neg->GetXaxis()->SetTitle("Invariant M_{#Pi^{+} K^{-}}");
  TwoBody2_Neg->GetYaxis()->SetTitle("Events");
  TwoBody2_Neg->SetStats(0);
  TwoBody2_Neg->SetLineColor(kRed);
  TwoBody2_Neg->Draw();


  c2->cd(2);
  c2->cd(2)->SetGrid();
  gPad->SetTicky(2);
  TwoBody2_Pos->SetTitle("M_{K^{+} #Pi^{-}} 2-Body Resonance Invariant mass Positive (no cuts)");
  TwoBody2_Pos->GetXaxis()->SetTitle("Invariant M_{K^{+} #Pi^{-}}");
  TwoBody2_Pos->GetYaxis()->SetTitle("Events");
  TwoBody2_Pos->SetStats(0);
  TwoBody2_Pos->SetLineColor(kBlue);
  TwoBody2_Pos->GetYaxis()->SetLabelOffset(0.01);
  TwoBody2_Pos->Draw();
  

  c2->SaveAs("2bodyResonanceJPsiErr.pdf");
*/ 
 
  TCanvas *c2 = new TCanvas("c2","multipads",900,700);
  
  gStyle->SetOptStat(0);
  c2->Divide(2,0);

  c2->cd(1);
  c2->cd(1)->SetGrid();
  gPad->SetTicky(2);
  TwoBody3_Neg->SetTitle("M_{#Pi^{+} #Pi^{-}} 2-Body Resonance Invariant mass Negative");
  TwoBody3_Neg->GetXaxis()->SetTitle("Invariant M_{#Pi^{+} #Pi^{-}}");
  TwoBody3_Neg->GetYaxis()->SetTitle("Events");
  TwoBody3_Neg->SetStats();
  TwoBody3_Neg->SetLineColor(kRed);
  TwoBody3_Neg->Draw();


  c2->cd(2);
  c2->cd(2)->SetGrid();
  gPad->SetTicky(2);
  TwoBody3_Pos->SetTitle("M_{#Pi^{+} #Pi^{-}} 2-Body Resonance Invariant mass Positive");
  TwoBody3_Pos->GetXaxis()->SetTitle("Invariant M_{#Pi^{+} #Pi^{-}}");
  TwoBody3_Pos->GetYaxis()->SetTitle("Events");
  TwoBody3_Pos->SetStats();
  TwoBody3_Pos->SetLineColor(kBlue);
  TwoBody3_Pos->GetYaxis()->SetLabelOffset(0.01);
  TwoBody3_Pos->Draw();
  

  c2->SaveAs("2bodyResonanceJPsiErr.pdf");



/// Dalitz canvas' for the plotting the background

  TCanvas *c11 = new TCanvas("c11_BG","B+_BG K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Pos_BG->SetTitle("B+_BG K#Pi VS #Pi#Pi");
  Dalitz23_Pos_BG->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Pos_BG->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Pos_BG->SetStats(0);
  Dalitz23_Pos_BG->Draw("colz"); // makes the plot colour scaled
  c11->SaveAs("B+_Background.pdf");
/*
  TCanvas *c12 = new TCanvas("c12_BG","B+_BG #Pi#Pi VS K#Pi",600,400);
  Dalitz31_Pos_BG->SetTitle("B+_BG #Pi#Pi VS K#Pi");
  Dalitz31_Pos_BG->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos_BG->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos_BG->SetStats(0);
  Dalitz31_Pos_BG->Draw("colz"); // makes the plot colour scaled
  c12->SaveAs("B+DalitzPlot_ KPi Vs PiPi_BG.pdf");

  TCanvas *c13 = new TCanvas("c13_BG","B+_BG K#Pi VS K#Pi",600,400);
  Dalitz21_Pos_BG->SetTitle("B+_BG K#Pi VS K#Pi");
  Dalitz21_Pos_BG->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Pos_BG->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Pos_BG->SetStats(0);
  Dalitz21_Pos_BG->Draw("colz"); // makes the plot colour scaled
  c13->SaveAs("B+DalitzPlot_ KPi Vs PiK_BG.pdf");
*/
  // Negative Dalitz:
  TCanvas *c17 = new TCanvas("c17_BG","B-_BG K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Neg_BG->SetTitle("B-_BG_K#Pi VS #Pi#Pi");
  Dalitz23_Neg_BG->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Neg_BG->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Neg_BG->SetStats(0);
  Dalitz23_Neg_BG->Draw("colz"); // makes the plot colour scaled
  c17->SaveAs("B-_Background.pdf");
/*
  TCanvas *c18 = new TCanvas("c18_BG","B-_BG #Pi#Pi VS K#Pi",600,400);
  Dalitz31_Neg_BG->SetTitle("B-_BG K#Pi VS #Pi#Pi");
  Dalitz31_Neg_BG->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg_BG->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg_BG->SetStats(0); Dalitz23_Pos->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Pos->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz31_Neg_BG->Draw("colz"); // makes the plot colour scaled
  c18->SaveAs("B-DalitzPlot_ KPi Vs PiPi_BG.pdf"); 

  TCanvas *c19 = new TCanvas("c19_BG","B-_BG K#Pi VS K#Pi",600,400);
  Dalitz21_Neg_BG->SetTitle("B-_BG K#Pi VS K#Pi");
  Dalitz21_Neg_BG->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Neg_BG->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Neg_BG->SetStats(0);
  Dalitz21_Neg_BG->Draw("colz"); // makes the plot colour scaled
  c19->SaveAs("B-DalitzPlot_ KPi Vs PiK.pdf_BG");
   */

  /// Subtracting the background from the signal peak
 
 
  /// Calculate the ratio for the scaling factor so the background can be successfully subtracted from the signal peak
  double Area_R1_Pos = backExp_Pos->Integral(5228, 5343);
  double Area_R2_Pos = backExp_Pos->Integral(5400, 5515);
  double Ratio_A_Pos = Area_R1_Pos / Area_R2_Pos;

  double Area_R1_Neg = backExp_Neg->Integral(5228, 5343);
  double Area_R2_Neg = backExp_Neg->Integral(5400, 5515);
  double Ratio_A_Neg = Area_R1_Neg / Area_R2_Pos;



  TCanvas *c20 = new TCanvas("c17_NoBG","B-_NoBG K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Neg->Sumw2();
  Dalitz23_Neg->SetTitle("B-_NoBG K#Pi VS #Pi#Pi");
  Dalitz23_Neg->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Neg->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Neg->SetStats(0);
  Dalitz23_Neg->Add(Dalitz23_Neg_BG,-Ratio_A_Neg);
    
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= xNBINS; binX++){
    for(int binY(1); binY <= yNBINS; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz23_Neg->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz23_Neg->SetBinContent(binX,binY,0);
      }
    }
  }		
  
  Dalitz23_Neg->Draw("colz"); // makes the plot colour scaled
  c20->SaveAs("B-_NoBackground.pdf");

  /*
  
  TCanvas *c21 = new TCanvas("c18_NoBG","B-_NoBG K#Pi VS K#Pi",600,400);
  Dalitz31_Neg->SetTitle("B-_NoBG K#Pi VS K#Pi");
  Dalitz31_Neg->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Neg->SetStats(0);
  Dalitz31_Neg->Add(Dalitz31_Neg_BG,-Ratio_A_Neg);
  
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= DalitzBins; binX++){
    for(int binY(1); binY <= DalitzBins; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz31_Neg->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz31_Neg->SetBinContent(binX,binY,0);
      }
    }
  }		
  
  Dalitz31_Neg->Draw("colz"); // makes the plot colour scaled
  c21->SaveAs("B-DalitzPlot_ KPi Vs PiPi.pdf");


  TCanvas *c22 = new TCanvas("c19_NoBG","B-_NoBG K#Pi Vs #PiK",600,400);
  Dalitz21_Neg->SetTitle("B-_NoBG K#Pi VS K#Pi");
  Dalitz21_Neg->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Neg->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Neg->SetStats(0);
  Dalitz21_Neg->Add(Dalitz21_Neg_BG,-Ratio_A_Neg);
  
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= DalitzBins; binX++){
    for(int binY(1); binY <= DalitzBins; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz21_Neg->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz21_Neg->SetBinContent(binX,binY,0);
      }
    }
  }		
  Dalitz21_Neg->Draw("colz"); // makes the plot colour scaled
  c22->SaveAs("B-DalitzNoBGPlot_ KPi Vs PiK.pdf_BG");
  */
  
  TCanvas *c23 = new TCanvas("c11_NoBG","B+_NoBG K#Pi VS #Pi#Pi",600,400);
  Dalitz23_Pos->Sumw2();
  Dalitz23_Pos->SetTitle("B+_NoBG K#Pi VS #Pi#Pi");
  Dalitz23_Pos->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz23_Pos->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Dalitz23_Pos->SetStats(0);
  Dalitz23_Pos->Add(Dalitz23_Pos_BG,-Ratio_A_Pos);
    
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= xNBINS; binX++){
    for(int binY(1); binY <= yNBINS; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz23_Pos->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz23_Pos->SetBinContent(binX,binY,0);
      }
    }
  }		
  
  //Dalitz23_Pos->removeNeg();
  Dalitz23_Pos->Draw("colz");
 // Dalitz23_Pos->Draw("colz"); // makes the plot colour scaled
  c23->SaveAs("B+_NoBackground.pdf");

  /*
  TCanvas *c24 = new TCanvas("c12_NoBG","B+_NoBG #Pi#Pi VS K#Pi",600,400);
  Dalitz31_Pos->SetTitle("B+_NoBG #Pi#Pi VS K#Pi");
  Dalitz31_Pos->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz31_Pos->SetStats(0);
  Dalitz31_Pos->Add(Dalitz31_Pos_BG,-Ratio_A_Pos);
      
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= DalitzBins; binX++){
    for(int binY(1); binY <= DalitzBins; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz31_Pos->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz31_Pos->SetBinContent(binX,binY,0);
      }
    }
  }		
  
  Dalitz31_Pos->Draw("colz"); // makes the plot colour scaled
  c24->SaveAs("B+DalitzPlot_ KPi Vs PiPi.pdf");
  
  TCanvas *c25 = new TCanvas("c13_NoBG","B+_NoBG K#Pi Vs K#Pi",600,400);
  Dalitz21_Pos->SetTitle("B+_NoBG K#Pi VS K#Pi");
  Dalitz21_Pos->GetXaxis()->SetTitle("M_{#Pi#K}^{2} (GeV)^{2}");
  Dalitz21_Pos->GetYaxis()->SetTitle("M_{K #Pi}^{2} (GeV)^{2}");
  Dalitz21_Pos->SetStats(0);
  Dalitz21_Pos->Add(Dalitz21_Pos_BG,-Ratio_A_Pos);
  
  // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= DalitzBins; binX++){
    for(int binY(1); binY <= DalitzBins; binY++){
      // Adding bin errors
      double binContent(0);
      binContent = Dalitz21_Pos->GetBinContent(binX,binY);
      if (binContent < 0){
	Dalitz21_Pos->SetBinContent(binX,binY,0);
      }
    }
  }		
  
  Dalitz21_Pos->Draw("colz"); // makes the plot colour scaled
  c25->SaveAs("B+DalitzNoBGPlot_ KPi Vs PiK.pdf_BG");
  
  */
  TH2F *Asym_Local_Dalitz = (TH2F*) Dalitz23_Neg ->GetAsymmetry(Dalitz23_Pos);
  
  TCanvas *c27 = new TCanvas("27","B-Meson Asymmetry Plot",600,400);
  Asym_Local_Dalitz->Sumw2();
  Asym_Local_Dalitz->SetTitle("B-Meson Asymmetry Plot");
  Asym_Local_Dalitz->GetXaxis()->SetTitle("M_{#Pi#Pi}^{2} (GeV)^{2}");
  Asym_Local_Dalitz->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Asym_Local_Dalitz->SetStats(0);
  Asym_Local_Dalitz->SetMaximum(1);
  Asym_Local_Dalitz->SetMinimum(-1);
  
  // Deleting bins which are due to statistical fluctuations
  
  Asym_Local_Dalitz->SetBinContent(8,15,-11);
  Asym_Local_Dalitz->SetBinContent(9,15,-11);
  Asym_Local_Dalitz->SetBinContent(7,17,-11);
  Asym_Local_Dalitz->SetBinContent(6,18,-11);
  Asym_Local_Dalitz->SetBinContent(12,12,-11);
  
   for(int binX(1); binX <= xNBINS; binX++){
    for(int binY(1); binY <= yNBINS; binY++){
      // Adding bin errors
      double binContent(0);
      double binError(0);
      binContent = Asym_Local_Dalitz->GetBinContent(binX,binY);
      binError = Asym_Local_Dalitz->GetBinError(binX,binY);
      
      if (binContent == 1 || binContent == -1){
	// Gets rid of background colours
	Asym_Local_Dalitz->SetBinContent(binX,binY,-11);// Setting to -11 as out of range of plot -> goes to white
      } 
 
      if (binContent == 0){
	// Gets rid of background colours
	Asym_Local_Dalitz->SetBinContent(binX,binY,-11);// Setting to -11 as out of range of plot -> goes to white
      }   
     /// Returing the local asymmetries
     // Asym_Local_Dalitz->SetBinContent(binX,binY,-11);
     double Asym_local(0);
     double A_error_local(0);
      double Asym_local = Asym_Local_Dalitz->GetBinContent(binX,binY); 
      double A_error_local = Asym_Local_Dalitz->GetBinError(binX,binY);
      
      if (binX == 6 && binY == 16){
	std::cout<<"Local Asymetries of reigons:"<<"("<<binX<<","<<binY<<") :" << Asym_local <<"+/-"<<A_error_local<< std::endl;
	Asym_Local_Dalitz->SetBinContent(binX,binY,-11);
      }
      
    }
  }
  

  
  Asym_Local_Dalitz->Draw("colz"); // makes the plot colour scaled
  c27->SaveAs("B-Meson Asymmetry Plot.pdf");
  
  TH2F *Asym_Sig_Local_Dalitz = (TH2F*) Asym_Local_Dalitz;
  

  
   // Loop over the bins to get rid of the negative ones
  for(int binX(1); binX <= xNBINS; binX++){
    for(int binY(1); binY <= yNBINS; binY++){
      // Adding bin errors
      double binContent(0);
      double binError(0);
      double SignificanceOfBinAsym(0);
      binContent = Asym_Sig_Local_Dalitz->GetBinContent(binX,binY); // 1 / percentage error = bin / error
      binError = Asym_Sig_Local_Dalitz->GetBinError(binX,binY);
      
      SignificanceOfBinAsym = binContent / binError;

      // Deleting bins 
	Asym_Sig_Local_Dalitz->SetBinContent(8,15,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(9,15,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(7,17,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(6,18,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(12,12,-11);
	// extra ones which significant plot only contains
	Asym_Sig_Local_Dalitz->SetBinContent(9,14,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(8,16,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(11,11,-11);
	Asym_Sig_Local_Dalitz->SetBinContent(11,4,-11);
   
      
      if (SignificanceOfBinAsym == 0){
	Asym_Sig_Local_Dalitz->SetBinContent(binX,binY,-11);
      }
      
      
      Asym_Sig_Local_Dalitz->SetBinContent(binX,binY,SignificanceOfBinAsym);
    }
  }	
  
  /// Local Asymmetries in Bplus
  TCanvas *c26 = new TCanvas("26","B-Meson Asymmetry Significance Plot",600,400);
  Asym_Sig_Local_Dalitz->Sumw2();
  Asym_Sig_Local_Dalitz->SetTitle("B-Meson Asymmetry Significance Plot");
  Asym_Sig_Local_Dalitz->SetMaximum(5);
  Asym_Sig_Local_Dalitz->SetMinimum(-5);
  Asym_Sig_Local_Dalitz->GetXaxis()->SetTitle("M_{#Pi #Pi}^{2} (GeV)^{2}");
  Asym_Sig_Local_Dalitz->GetYaxis()->SetTitle("M_{#Pi K}^{2} (GeV)^{2}");
  Asym_Sig_Local_Dalitz->SetStats(0);

  Asym_Sig_Local_Dalitz->Draw("colz"); // makes the plot colour scaled
  c26->SaveAs("B-Meson Asymmetry Significance Plot.pdf");
  
  /*
  /// Merging the two canvas' for B- and B+
  TCanvas *c3 = new TCanvas("c3","AsymmetryAndSignificant_multipads",900,700);
  
  gStyle->SetOptStat(0);
  c3->Divide(2,0);

  c3->cd(1);
  gPad->SetTicky(2);
  
  Asym_Local_Dalitz->Draw("colz");

  c3->cd(2);
  //gPad->SetTicky(2);
  
  Asym_Sig_Local_Dalitz->Draw("colz");
  
  c3->SaveAs("AsymmetryAndSignificant_multipads.pdf");
    */
}