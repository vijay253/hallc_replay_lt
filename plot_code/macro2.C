// Code for plotting the data
// Author... VIJAY KUMAR..... 
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"

int calib()

{
  gROOT->ForceStyle();

  TCanvas * plot = new TCanvas("plot", "Consistency of parameters"); plot->SetGrid();
  TMultiGraph* mg = new TMultiGraph(); mg->SetTitle("PMT 1");  // Add title of the plot

  TGraphErrors *gr1 = new TGraphErrors("input_pmt4_guess1.txt","%lg %lg %lg");
  gr1->SetMarkerColor(kRed);
  gr1->SetMarkerStyle(kFullCircle);
  gr1->SetLineColor(kBlue);

  TGraphErrors *gr2 = new  TGraphErrors("input_pmt4_guess2.txt","%lg %lg %lg");
  gr2->SetMarkerColor(kBlue);    
  gr2->SetMarkerStyle(kFullSquare);
  gr2->SetLineColor(kRed);

  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("apl");
  mg->GetXaxis()->SetTitle("Run Numbers");                // Add title of the X axis
  mg->GetYaxis()->SetTitle("Calibration Constants");     // Add title of the Y axis
  mg->GetHistogram()->SetMaximum(8);                    // Set the range of Y axis           
  mg->GetHistogram()->SetMinimum(5.2);

  TLegend* leg = new TLegend(.6,.7,.9,.9);            // Draw the Legend
  leg->SetFillColor(0);
  leg->AddEntry(&*gr1,"First Guess");
  leg->AddEntry(&*gr2,"Second Guess");
  leg->Draw("Same");                 
  plot->SaveAs("pmt4.png");                        // Save the plot
  return 0;                                       // Finish the function



}
