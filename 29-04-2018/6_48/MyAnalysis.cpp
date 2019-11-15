#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include<cmath>

///Define global variables
const double ProbKThreshold = 0.9;
const double ProbNotKThreshold = 0.2;
const double ProbPiThreshold = 0.6; // 0.6
const double ProbNotPiThreshold = 0.4;
const double KAON = 493.677;
const double PION = 139.57018;
const double D_0 = 1869.62;
double Mass1, Mass2, Mass3;
double M3_Pos, M2_Pos, M1_Pos;
double M3_Neg, M2_Neg, M1_Neg;
double M_Bmeson(0); // initialise the mass of the BMeson
double M_BmesonPos(0);
double M_BmesonNeg(0);

///Functions definitions

// Momentum modulus
double mvMod(double P_x, double P_y, double P_z){ 
  double modulus (0);
  modulus = pow(pow(P_x,2) + pow(P_y,2) + pow(P_z,2),0.5);  
  return modulus;
}

// Particle 4-Energy
double Energy_calc(double inv_mass, double P_x, double P_y, double P_z){
  return pow(pow(inv_mass,2) + pow(P_x,2) + pow(P_y,2) + pow(P_z,2),0.5); 
}

//3-Momentum dot product
double mv_3DotProd(double P1_x, double P1_y, double P1_z, double P2_x, double P2_y, double P2_z){
  double product(0);
  product = P1_x * P2_x + P1_y * P2_y + P1_z * P2_z;
  return product;
}

//Invariant Mass of B-Meson (3 body)
double inv_mass_Bmeson(double E1, double E2, double E3, double P12,double P13,double P23){
  double BmassSqr (0);
  BmassSqr = pow(KAON,2) + 2*pow(PION,2) + 2*(E1*E2 + E2*E3 + E1*E3 - P12 - P13 - P23);
  return sqrt(BmassSqr);
}

//Invariant Mass of D-Meson (2-body)
double invMass_D(double Mass1, double Mass2, double E_1, double E_2, double p_12){
  return sqrt(pow(Mass1,2) + pow(Mass2,2) + 2*(E_1*E_2 - p_12));
}

// A function to calculate the Chisqr(IP), the chi-squared of the impact parameter

/// Histograms

void MyAnalysis::BookHistos() {
// This function is only called once at the start of the program.
// Book your histograms here. The format is object_name, histogram_name, number_of_bins, minimum, maximum
// For a 2D histogram, use TH2F with first the number of bins and limits for the x axis and then for the y axis

/// Section 5.2 (Momentum X,Y,Z; Probablity Kaon; Probablity Pion; 2D Momenutum gradient.
/*
//Define Histograms
h_PX = new TH1F("h_PX", "", 100, -1e4, 1e4);
h_PY = new TH1F("h_PY", "", 100, -1e4, 1e4);
h_PZ = new TH1F("h_PZ", "", 100, -1e5, 1e5);
h_ProbK = new TH1F("h_ProbK", "", 100, 0, 1);
h_ProbPi = new TH1F("h_ProbPi", "", 100, 0, 1);
h_TXTY = new TH2F("h_TXTY", "", 100, -1, 1, 100, -1, 1);

//Add Histograms to a vector
v_Histos.push_back(h_PX);
v_Histos.push_back(h_PY);
v_Histos.push_back(h_PZ);
v_Histos.push_back(h_ProbK);
v_Histos.push_back(h_ProbPi);
v_Histos.push_back(h_TXTY);

*/

/// Section 5.3 B-Meson (3-body) Invariant mass
//B Meson (3-body) Invariant mass peak 5370.18,5479.28

// Define Histogram
//h_MassB = new TH1F("h_MassB", "", 100, 4028, 6364); 

//Add Histogram to a Vector
//v_Histos.push_back(h_MassB);


/// Section 5.4 D-Meson (2-body) Invariant mass 
//Define Histograms
h_MassM1_Pos = new TH1F("h_MassM1_Pos", "", 100, 100, 6000);
h_MassM2_Pos = new TH1F("h_MassM2_Pos", "", 100, 100, 6000);
h_MassM3_Pos = new TH1F("h_MassM3_Pos", "", 100, 100, 6000);

// Add all histograms to a vector
v_Histos.push_back(h_MassM1_Pos);
v_Histos.push_back(h_MassM2_Pos);
v_Histos.push_back(h_MassM3_Pos);

//Define Histograms
h_MassM1_Neg = new TH1F("h_MassM1_Neg", "", 100, 100, 6000);
h_MassM2_Neg = new TH1F("h_MassM2_Neg", "", 100, 100, 6000);
h_MassM3_Neg = new TH1F("h_MassM3_Neg", "", 100, 100, 6000);

// Add all histograms to a vector
v_Histos.push_back(h_MassM1_Neg);
v_Histos.push_back(h_MassM2_Neg);
v_Histos.push_back(h_MassM3_Neg);

// Define Histogram
h_MassB = new TH1F("h_MassB", "", 100, 4028, 6364); 

//Add Histogram to a Vector
v_Histos.push_back(h_MassB);

/// Section 5.5 Global Asymmetry (3-body) B+/- Invariant mass

// Define Histogram
h_MassBPos = new TH1F("h_MassBPos", "", 40, 5040, 5600);
h_MassBNeg = new TH1F("h_MassBNeg", "", 40, 5040, 5600); 

//Add Histogram to a Vector
v_Histos.push_back(h_MassBPos);
v_Histos.push_back(h_MassBNeg);

/// Dalitz bins
  const Int_t xNBINS = 90;// 23
  //Double_t xedges[xNBINS + 1] = {-1., 0., 0.7, 1, 1.2,2.,6.,10,15.,20.,24.,26.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.};
  const Int_t yNBINS = 2;
 // Double_t yedges[yNBINS + 1] = {-1.,0.,0.2,0.4,0.5,1.,1.2,1.7,1.9,3.0,3.8,5.,12.5,14.,15.5,22.,25.,27.,31.,32.,33.,34.,35.,36.,37.};
/*
/// Section 5.6 Dalitz plots (for the signal peak)
h_Dalitz31_Pos = new TH2F("h_Dalitz31_Pos", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz31_Pos);
h_Dalitz23_Pos = new TH2F("h_Dalitz23_Pos", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz23_Pos);
h_Dalitz21_Pos = new TH2F("h_Dalitz21_Pos", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz21_Pos);

h_Dalitz31_Neg = new TH2F("h_Dalitz31_Neg", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz31_Neg);
h_Dalitz23_Neg = new TH2F("h_Dalitz23_Neg", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz23_Neg);
h_Dalitz21_Neg = new TH2F("h_Dalitz21_Neg", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz21_Neg);

h_MassBPos_BG = new TH1F("h_MassBPos_BG", "", 300, 5010, 6000);
h_MassBNeg_BG = new TH1F("h_MassBNeg_BG", "", 300, 5010, 6000); 

//Add Histogram to a Vector
v_Histos.push_back(h_MassBPos);
v_Histos.push_back(h_MassBNeg);

/// Section 5.6 Dalitz Plots for the background subtraction
h_Dalitz31_Pos_BG = new TH2F("h_Dalitz31_Pos_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz31_Pos_BG);
h_Dalitz23_Pos_BG = new TH2F("h_Dalitz23_Pos_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz23_Pos_BG);
h_Dalitz21_Pos_BG = new TH2F("h_Dalitz21_Pos_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz21_Pos_BG);

h_Dalitz31_Neg_BG = new TH2F("h_Dalitz31_Neg_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz31_Neg_BG);
h_Dalitz23_Neg_BG = new TH2F("h_Dalitz23_Neg_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz23_Neg_BG);
h_Dalitz21_Neg_BG = new TH2F("h_Dalitz21_Neg_BG", "", xNBINS, xedges, yNBINS, yedges);
v_Histos.push_back(h_Dalitz21_Neg_BG);
*/

/// Section 5.6 Dalitz plots (for the signal peak)
h_Dalitz31_Pos = new TH2F("h_Dalitz31_Pos", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz31_Pos);
h_Dalitz23_Pos = new TH2F("h_Dalitz23_Pos", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz23_Pos);
h_Dalitz21_Pos = new TH2F("h_Dalitz21_Pos", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz21_Pos);

h_Dalitz31_Neg = new TH2F("h_Dalitz31_Neg", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz31_Neg);
h_Dalitz23_Neg = new TH2F("h_Dalitz23_Neg", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz23_Neg);
h_Dalitz21_Neg = new TH2F("h_Dalitz21_Neg", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz21_Neg);

h_MassBPos_BG = new TH1F("h_MassBPos_BG", "", 300, 5010, 6000);
h_MassBNeg_BG = new TH1F("h_MassBNeg_BG", "", 300, 5010, 6000); 

//Add Histogram to a Vector
v_Histos.push_back(h_MassBPos);
v_Histos.push_back(h_MassBNeg);

/// Section 5.6 Dalitz Plots for the background subtraction
h_Dalitz31_Pos_BG = new TH2F("h_Dalitz31_Pos_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz31_Pos_BG);
h_Dalitz23_Pos_BG = new TH2F("h_Dalitz23_Pos_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz23_Pos_BG);
h_Dalitz21_Pos_BG = new TH2F("h_Dalitz21_Pos_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz21_Pos_BG);

h_Dalitz31_Neg_BG = new TH2F("h_Dalitz31_Neg_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz31_Neg_BG);
h_Dalitz23_Neg_BG = new TH2F("h_Dalitz23_Neg_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz23_Neg_BG);
h_Dalitz21_Neg_BG = new TH2F("h_Dalitz21_Neg_BG", "", xNBINS, -1, 30, yNBINS, -1, 30);
v_Histos.push_back(h_Dalitz21_Neg_BG);
}
/// Data selection

Bool_t MyAnalysis::Cut() {
///resonant Identifier 

  if ( H3_PZ < 0 ) return false;
 
   
  if (H2_Charge == +1 && H2_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if ((H1_Charge == -1 && H3_Charge == +1) || (H3_Charge == -1 && H1_Charge == +1)){
      // par	ticle 2 is a kaon, other 2 are pions
      Mass2 = KAON;
      Mass3 = Mass1 = PION;
      
     // std::cout<<"Pos1:"<<H1_Charge<<" "<<H2_Charge<<" "<<H3_Charge<<std::endl;
      
      
      double P1_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E1_2 = Energy_calc(KAON,H2_PX, H2_PY,H2_PZ);
      M1_Pos = invMass_D(KAON, PION,E1_1,E1_2, P1_12); 

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E2_1 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(KAON,H2_PX, H2_PY, H2_PZ);
      M2_Pos = invMass_D(PION,KAON,E2_1,E2_2, P2_12);

      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double E3_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      M3_Pos = invMass_D(PION,PION,E3_1, E3_2, P3_12);
    }
  }

  // Negative
  else if (H2_Charge == -1 && H2_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if ((H1_Charge == -1 && H3_Charge == +1) || (H3_Charge == -1 && H1_Charge == +1)){
      // particle 2 is a kaon, other 2 are pions
      Mass2 = KAON;
      Mass3 = Mass1 = PION;
      //std::cout<<"Neg1:"<<H1_Charge<<" "<<H2_Charge<<" "<<H3_Charge<<std::endl;
      double P1_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E1_2 = Energy_calc(KAON,H2_PX, H2_PY,H2_PZ);
      M1_Neg = invMass_D(KAON, PION,E1_1,E1_2, P1_12); 

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E2_1 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(KAON,H2_PX, H2_PY, H2_PZ);
      M2_Neg = invMass_D(PION,KAON,E2_1,E2_2, P2_12);

      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double E3_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      M3_Neg = invMass_D(PION,PION,E3_1, E3_2, P3_12);

    }
  }

  // Positive 
  else if (H1_Charge == +1 && H1_ProbK > ProbKThreshold && H2_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if ((H2_Charge == -1 && H3_Charge == +1) || (H3_Charge == -1 && H2_Charge == +1)){
      // check that particle 1 is a kaon and particles 2 and 3 are pions
      
      Mass1 = KAON; // set the mass of particle 1 to that of a kaon
      Mass2 = Mass3 = PION; // " " pion
     // std::cout<<"Pos2:"<<H1_Charge<<" "<<H2_Charge<<" "<<H3_Charge<<std::endl;
      double P1_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
      double E1_2 = Energy_calc(PION,H2_PX, H2_PY,H2_PZ);
      M1_Pos = invMass_D(KAON, PION,E1_1,E1_2, P1_12);

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E2_1 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(PION,H2_PX, H2_PY, H2_PZ);
      M3_Pos = invMass_D(PION,PION,E2_1,E2_2, P2_12);

      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double E3_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      M2_Pos = invMass_D(KAON,PION,E3_1, E3_2, P3_12);
    }
  }

  // Negative
  else if (H1_Charge == -1 && H1_ProbK > ProbKThreshold && H2_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if ((H2_Charge == -1 && H3_Charge == +1) || (H3_Charge == -1 && H2_Charge == +1)){
      // check that particle 1 is a kaon and particles 2 and 3 are pions
      
      Mass1 = KAON; // set the mass of particle 1 to that of a kaon
      Mass2 = Mass3 = PION; // " " pion
      // std::cout<<"Neg2:"<<H1_Charge<<" "<<H2_Charge<<" "<<H3_Charge<<std::endl;
      double P1_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
      double E1_2 = Energy_calc(PION,H2_PX, H2_PY,H2_PZ);
      M1_Neg = invMass_D(KAON, PION,E1_1,E1_2, P1_12);

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E2_1 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(PION,H2_PX, H2_PY, H2_PZ);
      M3_Neg = invMass_D(PION,PION,E2_1,E2_2, P2_12);

      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double E3_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
      M2_Neg = invMass_D(KAON,PION,E3_1, E3_2, P3_12);
    }
  }

  //(3) Cut for the J/Psi Meson 
  // Positive
  else if (H3_Charge == +1 && H3_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H2_ProbPi > ProbPiThreshold){
    if ((H1_Charge == -1 && H2_Charge == +1) || (H2_Charge == -1 && H1_Charge == +1)){
      Mass3 = KAON;
      Mass1 = Mass2 = PION;
      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E3_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H2_PX, H2_PY, H2_PZ);
      M3_Pos = invMass_D(PION,PION,E3_1, E3_2, P3_12);

      double P1_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(KAON,H3_PX, H3_PY, H3_PZ);
      double E1_2 = Energy_calc(PION,H2_PX, H2_PY,H2_PZ);
      M1_Pos = invMass_D(KAON, PION,E1_1,E1_2, P1_12);

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H1_PX, H1_PY, H1_PZ);
      double E2_1 = Energy_calc(KAON,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      M2_Pos = invMass_D(PION,KAON,E2_1,E2_2, P2_12);
    }  
  }

  // Negative
  else if (H3_Charge == -1 && H3_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H2_ProbPi > ProbPiThreshold){
    if ((H1_Charge == -1 && H2_Charge == +1) || (H2_Charge == -1 && H1_Charge == +1)){

      Mass3 = KAON;
      Mass1 = Mass2 = PION;
      double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double E3_1 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      double E3_2 = Energy_calc(PION,H2_PX, H2_PY, H2_PZ);
      M3_Neg = invMass_D(PION,PION,E3_1, E3_2, P3_12);

      double P1_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
      double E1_1 = Energy_calc(KAON,H3_PX, H3_PY, H3_PZ);
      double E1_2 = Energy_calc(PION,H2_PX, H2_PY,H2_PZ);
      M1_Neg = invMass_D(KAON, PION,E1_1,E1_2, P1_12);

      double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H1_PX, H1_PY, H1_PZ);
      double E2_1 = Energy_calc(KAON,H3_PX, H3_PY, H3_PZ);
      double E2_2 = Energy_calc(PION,H1_PX, H1_PY, H1_PZ);
      M2_Neg = invMass_D(PION,KAON,E2_1,E2_2, P2_12);
    }	  
  }
  else{
    return false;
  }

  /// Section 5.4: Implementing resonant cuts
  // Cutting 3 standard deviations around the mean of the D meson 
  /*
  if((M2_Pos < 3000) || (M2_Pos > 15000){
    return false
  }
  */
  /*
  if ((M1_Pos > 1806.22 && M1_Pos < 1905.80) || (M2_Pos > 1732.22 && M2_Pos < 1905.80)){
    return false;
  }
  */
  if ((M1_Pos > 1806.22 && M1_Pos < 1905.80) || (M2_Pos > 1732.22 && M2_Pos < 1905.80)){
    return false;
  }

  if ((M1_Neg > 1806.22 && M1_Neg < 1905.80) || (M2_Neg > 1732.22 && M2_Neg < 1905.80)){ 
    return false;
  }
  

  // J/Psi inv mass cut
  //else if (M3_D > 3015.76 && M3_D < 3156.90){
  //	return false;
  // 	}

  //not implement 4 body decay as it was insignificant 
  //      else if(M3_D > 3615 && M3_D < 3777){
  //	return false;
  //      }

  /// Muon cut      
  if (H1_isMuon != 0 || H2_isMuon != 0 || H3_isMuon != 0){
      return false;
    }

/// Section 5.5 Global CP Asymmetry

// Positive B Meson invariant mass 
// Check one particle is positive, 2 are negative

  if (H1_Charge == +1 && H1_ProbK > ProbKThreshold && H2_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if (H2_Charge == -1 && H3_Charge == +1 || H3_Charge == -1 && H2_Charge == +1){
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      
      double E1 = Energy_calc(Mass1, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

      M_BmesonPos = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);

    }
  }

  else if (H2_Charge == +1 && H2_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if (H1_Charge == -1 && H3_Charge == +1 || H3_Charge == -1 && H1_Charge == +1){
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      double E1 = Energy_calc(Mass2, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass1, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

      M_BmesonPos = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);
    }
  }

  else if (H3_Charge == +1 && H3_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H2_ProbPi > ProbPiThreshold){
    if (H1_Charge == -1 && H2_Charge == +1 || H2_Charge == -1 && H1_Charge == +1){
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      double E1 = Energy_calc(Mass3, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass1, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

      M_BmesonPos = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);
    }
  }

  // Negative B meson invariant mass

  else if (H1_Charge == -1 && H1_ProbK > ProbKThreshold && H2_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if (H2_Charge == -1 && H3_Charge == +1 || H3_Charge == -1 && H2_Charge == +1){
      
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      double E1 = Energy_calc(Mass1, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

      M_BmesonNeg = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);
    }
  }

  else if (H2_Charge == -1 && H2_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H3_ProbPi > ProbPiThreshold){
    if (H1_Charge == -1 && H3_Charge == +1 || H3_Charge == -1 && H1_Charge == +1){
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      double E1 = Energy_calc(Mass2, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass1, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

      M_BmesonNeg = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);
    }
  }

  else if (H3_Charge == -1 && H3_ProbK > ProbKThreshold && H1_ProbPi > ProbPiThreshold && H2_ProbPi > ProbPiThreshold){
    if (H1_Charge == -1 && H2_Charge == +1 || H2_Charge == -1 && H1_Charge == +1){
      Mass1 = KAON;
      Mass3 = Mass2 = PION;
      double E1 = Energy_calc(Mass3, H1_PX, H1_PY, H1_PZ);
      double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
      double E3 = Energy_calc(Mass1, H3_PX, H3_PY, H3_PZ);
      double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
      double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
      double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);
      M_BmesonNeg = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23);
    }
  }

  else{
    return false;
  }   

// Making the data accurate to the pre-selection cuts

if (M_BmesonPos < 5005 || M_BmesonPos > 6300){
return false;
}

if (M_BmesonNeg < 5005 || M_BmesonNeg > 6300){
return false;
}

/*
  /// Dalitz plots: Making cuts around the BMeson peak, gets rid of noise
  if (M_BmesonPos < 5228 || M_BmesonPos > 5343){
    return false;
  }

  if (M_BmesonNeg < 5224 || M_BmesonNeg > 5345){
    return false;
  }


/// Dalitz plots: Making cuts around a section of the background, to be subtracted from under the signal peak
  if (M_BmesonPos_BG < 5400 || M_BmesonPos_BG > 5515){ // 5360 - 5475
    return false;
  }

  if (M_BmesonNeg_BG < 5400 || M_BmesonNeg_BG > 5521){ // 5360 - 5481
    return false;
  }
*/
//Note: Do not include if we want to look at all channels
// else{
//return false;}

/// Section 5.7: producing dalitz plots to our best display

// Selecting only the reigon of the dalitz plot which we believe is ok

  if (((M3_Pos*M3_Pos)/1e+6) >0.8){
    return false;
  }
  
  if (((M2_Pos*M2_Pos)/1e+6)>22) {
    return false;
  }

  /*
  if (((M3_Neg*M3_Neg)/1e+6)>0.7){
    return false;
  }

  if (((M2_Neg*M2_Neg)/1e+6)>12.5 || (((M2_Neg*M2_Neg)/1e+6)<3.0)){
    return false;
  }
*/
 /* if (M3_Neg >800 || M3_Neg <750){
    return false;
  }
  if (M3_Pos >800 || M3_Pos <750){
    return false;
  }
 */ 
 /*
 /// Cutting M3 (the phi)
   if ((M3_Pos < 1030 )&& (M3_Pos > 1010)){ // 950-1100
    return false;
  }
   if ((M3_Neg >1010) && (M3_Neg <1030)){
    return false;
  }
  */
  return true;
}

/// Fill the histograms with data

void MyAnalysis::Execute() {
  //Description:
  // This method gets called on every event.
  // Call the Cut function to decide whether to plot this event or not
  // it returns if the cut function returns false
  if ( !Cut() ) return;

  // Fill your histograms below.
  /*  //Section 5.2
  h_PX->Fill( H1_PX );
  h_PX->Fill( H2_PX );
  h_PX->Fill( H3_PX );
  // the PY of all three particles
  h_PY->Fill( H1_PY );
  h_PY->Fill( H2_PY );
  h_PY->Fill( H3_PY );
  // the PZ of all three particles
  h_PZ->Fill( H1_PZ );
  h_PZ->Fill( H2_PZ );
  h_PZ->Fill( H3_PZ );
  // 2D histogram of PX/PZ vs PY/PZ
  h_TXTY->Fill( H1_PX / H1_PZ, H1_PY / H1_PZ );
  // the probability of Kaon of all 3 particles
  h_ProbK->Fill( H1_ProbK );
  h_ProbK->Fill( H2_ProbK );
  h_ProbK->Fill( H3_ProbK );
  // the probability of Kaon of all 3 particles
  h_ProbPi->Fill( H1_ProbPi );
  h_ProbPi->Fill( H2_ProbPi );
  h_ProbPi->Fill( H3_ProbPi );

  */

  /// Section 5.3 : Three Body B-Meson spectrum
  /*
  double E1 = Energy_calc(Mass1, H1_PX, H1_PY, H1_PZ);
  double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
  double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
  double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
  double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
  double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);
  */

  ///Section 5.3  B-Meson histo fit

  // M_Bmeson = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23); // change these to the values to execute
  // the invariant mass of the B meson
  //  h_MassB->Fill( M_Bmeson );

  ///Section 5.4 D-meson histo fit
  // the invariant mass of the resonant meson
  h_MassM1_Pos->Fill(M1_Pos);// M1_Pos is the Pi-, K+
  h_MassM2_Pos->Fill(M2_Pos);// M2_Pos is the Pi+, K+
  h_MassM3_Pos->Fill(M3_Pos);// M3_Pos is the Pi-, Pi+
  h_MassM1_Neg->Fill(M1_Neg);// M1_Neg is the Pi+, K-
  h_MassM2_Neg->Fill(M2_Neg);// M2_Neg is the Pi-, K-
  h_MassM3_Neg->Fill(M3_Neg);// M3_Neg is the Pi-, Pi+
  
  ///Section 5.4 B-Meson without D-Meson and J/PsiMeson
  /*  
  double E1 = Energy_calc(Mass1, H1_PX, H1_PY, H1_PZ);
  double E2 = Energy_calc(Mass2, H2_PX, H2_PY, H2_PZ);
  double E3 = Energy_calc(Mass3, H3_PX, H3_PY, H3_PZ);
  double P12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
  double P13 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
  double P23 = mv_3DotProd(H2_PX, H2_PY, H2_PZ, H3_PX, H3_PY, H3_PZ);

  M_Bmeson = inv_mass_Bmeson(E1, E2, E3, P12, P13, P23); // change these to the values to execute
  // the invariant mass of the B meson
  h_MassB->Fill( M_Bmeson );
  */  
  /// Section 5.5 Global Asymmetry (3-body) B+/- Invariant mass


  h_MassBPos->Fill( M_BmesonPos );

  h_MassBNeg->Fill( M_BmesonNeg );


  /// Section 5.6 Dalitz plots
  //h_Dalitz->Fill((M3_D*M3_D)/1e+6, (M1_D*M1_D)/1e+6);


  /*WIP: double P1_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H2_PX, H2_PY, H2_PZ);
  double E1_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
  double E1_2 = Energy_calc(PION,H2_PX, H2_PY,H2_PZ);
  M1_D = invMass_D(KAON, PION,E1_1,E1_2, P1_12);

  double P2_12 = mv_3DotProd(H3_PX, H3_PY, H3_PZ, H2_PX, H2_PY, H2_PZ);
  double E2_1 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
  double E2_2 = Energy_calc(PION,H2_PX, H2_PY, H2_PZ);
  M2_D = invMass_D(PION,PION,E2_1,E2_2, P2_12);

  double P3_12 = mv_3DotProd(H1_PX, H1_PY, H1_PZ, H3_PX, H3_PY, H3_PZ);
  double E3_1 = Energy_calc(KAON,H1_PX, H1_PY, H1_PZ);
  double E3_2 = Energy_calc(PION,H3_PX, H3_PY, H3_PZ);
  M3_D = invMass_D(KAON,PION,E3_1, E3_2, P3_12);
  */
  
  /// Section 5.6 Dalitz plots
  //signal peaks 
  if (M_BmesonPos > 5224.8 && M_BmesonPos < 5342.4){ // deltaX = 117.6
    h_Dalitz31_Pos->Fill((M3_Pos*M3_Pos)/1e+6, (M1_Pos*M1_Pos)/1e+6);
    h_Dalitz23_Pos->Fill((M3_Pos*M3_Pos)/1e+6, (M2_Pos*M2_Pos)/1e+6);
    h_Dalitz21_Pos->Fill((M2_Pos*M2_Pos)/1e+6, (M3_Pos*M3_Pos)/1e+6);
  }

  if (M_BmesonNeg > 5225.1 && M_BmesonNeg < 5343.2){ // deltaX = 118.1
    h_Dalitz31_Neg->Fill((M3_Neg*M3_Neg)/1e+6, (M1_Neg*M1_Neg)/1e+6);
    h_Dalitz23_Neg->Fill((M3_Neg*M3_Neg)/1e+6, (M2_Neg*M2_Neg)/1e+6);
    h_Dalitz21_Neg->Fill((M2_Neg*M2_Neg)/1e+6, (M3_Neg*M3_Neg)/1e+6);
  }
  
  // Approximating the background
   if (M_BmesonPos > 5350 && M_BmesonPos < 5467.6){ // 5360 - 5475
    h_Dalitz31_Pos_BG->Fill((M3_Pos*M3_Pos)/1e+6, (M1_Pos*M1_Pos)/1e+6);
    h_Dalitz23_Pos_BG->Fill((M3_Pos*M3_Pos)/1e+6, (M2_Pos*M2_Pos)/1e+6);
    h_Dalitz21_Pos_BG->Fill((M2_Pos*M2_Pos)/1e+6, (M3_Pos*M3_Pos)/1e+6);
  }

  if (M_BmesonNeg > 5350 && M_BmesonNeg < 5468.1){ // 5360 - 5481
    h_Dalitz31_Neg_BG->Fill((M3_Neg*M3_Neg)/1e+6, (M1_Neg*M1_Neg)/1e+6);
    h_Dalitz23_Neg_BG->Fill((M3_Neg*M3_Neg)/1e+6, (M2_Neg*M2_Neg)/1e+6);
    h_Dalitz21_Neg_BG->Fill((M2_Neg*M2_Neg)/1e+6, (M3_Neg*M3_Neg)/1e+6);
  }
  
}