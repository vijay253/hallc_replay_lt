
#define calibration_cxx
// The class definition in calibration.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("calibration.C+")
// Root > T->Process("calibration.C+","some options")
//


#include "calibration.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TSpectrum.h>
#include <TList.h>
#include <TPolyMarker.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TPaveText.h>
#include <TLatex.h>

using namespace TMath;

void calibration::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  printf("\n\n");

  TString option = GetOption();
  TString report_option = option(0,option.Length()-79);
  Info("Begin", "Script will fail unless 'calibration.C+' is used");
  Info("Begin", "Starting calibration process with option: %s", report_option.Data());
  Info("Begin", "To see details of calibration, use option showall");
  Info("Begin", "To calibrate using TrackFired leaf, use option trackfired");
  Info("Begin", "Default is no particle cut, use option cut if desired");
  Info("Begin", "Default particle ID is electrons, use option pions if desired");
  printf("\n\n");

  //Check option
  if (option.Contains("showall")) fFullShow = kTRUE;
  if (option.Contains("trackfired")) fTrack = kTRUE;
  if (option.Contains("pions") || option.Contains("pion")) fPions = kTRUE;
  if (option.Contains("cut") || fPions || option.Contains("cuts")) fCut = kTRUE; 
}

void calibration::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  printf("\n\n");
  TString option = GetOption();
   
  //Check option
  if (option.Contains("showall")) fFullShow = kTRUE;
  if (option.Contains("trackfired")) fTrack = kTRUE;
  if (option.Contains("pions") || option.Contains("pion")) fPions = kTRUE;
  if (option.Contains("cut") || fPions || option.Contains("cuts")) fCut = kTRUE;
  
  Info("SlaveBegin", "'%s' showing", (fFullShow ? "full" : "minimal"));                           
  Info("SlaveBegin", "'%s' strategy", (fTrack ? "tracking" : "quadrant"));
  Info("SlaveBegin", "cuts %s performed", (fCut ? "are" : "are not"));
  if (fCut) Info("SlaveBegin", "cutting for '%s'", (fPions ? "pions" : "electrons"));

  // Inintialize the histograms. Note they are binned per ADC channel which will be changed in the calibration analysis.
  Int_t ADC_min;
  Int_t ADC_max;
  Int_t bins;

  ADC_min = 0;
  ADC_max = 200;
  bins = 2*(abs(ADC_min) + abs(ADC_max));

  fPulseInt = new TH1F*[4];
  fPulseInt_quad = new TH1F**[4];
  fPulseInt_quad[0] = new TH1F*[4];
  fPulseInt_quad[1] = new TH1F*[4];
  fPulseInt_quad[2] = new TH1F*[4];
  fPulseInt_quad[3] = new TH1F*[4];

  for (Int_t ipmt = 0; ipmt < 4; ipmt++)
    {

      fPulseInt[ipmt] = new TH1F(Form("PulseInt_PMT%d",ipmt+1),Form("Pulse Integral PMT%d; ADC Channel (pC); Counts",ipmt+1), bins, ADC_min, ADC_max);
      GetOutputList()->Add(fPulseInt[ipmt]);

      for (Int_t iquad = 0; iquad < 4; iquad++)
	{
	  fPulseInt_quad[iquad][ipmt] = new TH1F(Form("PulseInt_quad%d_PMT%d",iquad+1,ipmt+1),Form("Pulse Integral PMT%d quad%d; ADC Channel (pC); Counts",ipmt+1,iquad+1),bins,ADC_min,ADC_max);
	  GetOutputList()->Add(fPulseInt_quad[iquad][ipmt]);
     
	}
    }


  fTim1 = new TH1F("Timing_PMT1", "ADC TDC Diff PMT1 ; Time (ns) ;Counts", 400, -40.0, 40.0);
  GetOutputList()->Add(fTim1);
  fTim2 = new TH1F("Timing_PMT2", "ADC TDC Diff PMT2 ; Time (ns) ;Counts", 400, -40.0, 40.0);
  GetOutputList()->Add(fTim2);
  fTim3 = new TH1F("Timing_PMT3", "ADC TDC Diff PMT3 ; Time (ns) ;Counts", 400, -40.0, 40.0);
  GetOutputList()->Add(fTim3);
  fTim4 = new TH1F("Timing_PMT4", "ADC TDC Diff PMT4 ; Time (ns) ;Counts", 400, -40.0, 40.0);
  GetOutputList()->Add(fTim4);

  //Timing and Beta cut visualizations
  fBeta_Cut = new TH1F("Beta_Cut", "Beta cut used for 'good' hits;Beta;Counts", 100, -0.1, 1.5);
  GetOutputList()->Add(fBeta_Cut); 
  
  fBeta_Full = new TH1F("Beta_Full", "Full beta for events;Beta;Counts", 100, -0.1, 1.5);
  GetOutputList()->Add(fBeta_Full);

  // fTim = new TH1F("Tim", "Time ; Time ;Counts", 100, -50.0, 0.0);
  // GetOutputList()->Add(fTim);

  fTiming_Cut = new TH1F("Timing_Cut", "Timing cut used for 'good' hits;Time (ns);Counts", 500, -50, 50);
  GetOutputList()->Add(fTiming_Cut);

  //fTiming_Cut = new TH2F("Timing_Cut", "Timing cut used for good hits ; Time (ns); Counts", 500, -50, 50, 500, 0.0, 200);
  //GetOutputList()->Add(fTiming_Cut);

  fTiming_Full = new TH1F("Timing_Full", "Full timing information for events;Time (ns);Counts", 500, -50, 50);
  GetOutputList()->Add(fTiming_Full);

  // fTiming_Full = new TH2F("Timing_Full", "Full timing information for events ;Time (ns); Counts", 500, -30, 40, 500, 0.0, 200);
  // GetOutputList()->Add(fTiming_Full);

  //Particle ID cut visualization

  fCut_everything = new TH2F("Cut_everything", "Visualization of no cuts; Normalized  Energy; Pre-Shower Energy (GeV)", 250, 0, 1.5, 250, 0, 1.0);
  GetOutputList()->Add(fCut_everything);
  fCut_enorm = new TH1F("Cut_enorm", "Visualization of normalized energy cuts; Normalized Energy; Counts", 200, 0, 2.0);
  GetOutputList()->Add(fCut_enorm);
  fCut_electron = new TH2F("Cut_electron", "Visualization of electron cut; Calorimeter Energy (GeV); Pre-Shower Energy (GeV)", 250, 0, 1.5, 250, 0, 1.0);
  GetOutputList()->Add(fCut_electron);
  fCut_pion = new TH2F("Cut_pion", "Visualization of pion cut; Normalized Energy; Pre-Shower Energy (GeV)", 250, 0, 1.0, 250, 0, 1.0);
  GetOutputList()->Add(fCut_pion);

  printf("\n\n");
}

Bool_t calibration::Process(Long64_t entry) 
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either calibration::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);
  
  //Output to verify script is working, and store the total number of events
  if (entry % 100000 == 0) printf("Processing Entry number %lld\n",entry);                          

  //Define quantities to loop over
  Int_t fpmts;
  fpmts = fhgc_pmts;   //Note HGC & NGC have the same # of PMTS

  //Require only one good track reconstruction for the event                         
  if (*Ndata_P_tr_beta != 1) return kTRUE;
  
  //Redundant, but useful if multiple tracks are eventually allowed
  for (Int_t itrack = 0; itrack < *Ndata_P_tr_beta; itrack++) 
    {
      //Require loose cut on particle velocity                                     
      fBeta_Full->Fill(P_tr_beta[itrack]);
      if (TMath::Abs(P_tr_beta[itrack] - 1.0) > 0.4) return kTRUE;
      fBeta_Cut->Fill(P_tr_beta[itrack]);

      //Filling the histograms
      for (Int_t ipmt = 0; ipmt < fpmts; ipmt++) 
	{	  
	  //Perform a loose timing cut    
	  fTiming_Full->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);       //doubt here

	  if(ipmt ==0){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >13 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 9) continue;

	    fTim1->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);

	  }

	  if(ipmt ==1){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >12 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 7) continue;

	    fTim2->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	  }
	  if(ipmt ==2){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >12 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 7) continue;

	    fTim3->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	  }
	  if(ipmt ==3){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >12 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 8) continue;

	    fTim4->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	  }
	  // cut modified by VK, 24/05/19

	  // if(P_hgcer_goodAdcTdcDiffTime[ipmt] >20 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 4) continue;
	   
	  // fTiming_Cut->Fill(P_hgcer_xAtCer[ipmt],P_hgcer_yAtCer[ipmt]);

	  // fTiming_Cut->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
	   
	  //Cuts to remove entries corresponding to a PMT not registering a hit    
	  if (P_hgcer_goodAdcPulseInt[ipmt] == 0.0) continue;
	 	  
	  //For quadrant cut strategy with particle ID cuts. In this case electrons are selected
	  if (!fTrack && fCut && !fPions)
	    {
	      //Retrieve particle ID information

	      //  Float_t central_p = 8.035;             // old value 6.0530


	      //  Float_t p = ((P_gtr_dp[0]/100.0)*central_p) + central_p;

	      Double_t p = *P_gtr_p;

	      //Fill histogram visualizaing the electron selection
	      fCut_everything->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);                  
	      fCut_enorm->Fill(*P_cal_etotnorm);

	      //Cut on Shower vs preshower is a tilted ellipse, this requires an angle of rotation (in radians), x/y center, semimajor and semiminor axis
	      //Float_t eangle = 3.0*3.14159/4.0;
	      //Float_t ex_center = 0.66;
	      //Float_t ey_center = 0.35;
	      //Float_t esemimajor_axis = 0.28;
	      //Float_t esemiminor_axis = 0.04;

	      Float_t eangle = 3.0*3.14159/4.0;                             
	      Float_t ex_center = 0.65;                                  // Old value 0.375
	      Float_t ey_center = 0.35;                                 // old value 0.360
	      Float_t esemimajor_axis = 0.30;
	      Float_t esemiminor_axis = 0.08;
	      if (pow((*P_cal_fly_earray/p - ex_center)*cos(eangle) + (*P_cal_pr_eplane/p - ey_center)*sin(eangle),2)/pow(esemimajor_axis,2) + 
		  pow((*P_cal_fly_earray/p - ex_center)*sin(eangle) - (*P_cal_pr_eplane/p - ey_center)*cos(eangle),2)/pow(esemiminor_axis,2) < 1
		  /* P_cal_etotnorm > 0.4*/)
		{
		  //Fill histogram visualizing the electron selection
		  fCut_electron->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

		  //Fill histogram of the full PulseInt spectra for each PMT
		  fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Fill histograms of what each PMT registers from each quadrant, this requires tracking the particle from the focal plane. Each quadrant is defined from the parameter files
		  Float_t y_pos = P_tr_y[0] + P_tr_ph[0]*fhgc_zpos;                              
		  Float_t x_pos = P_tr_x[0] + P_tr_th[0]*fhgc_zpos;
		  
		  //Condition for quadrant 1 mirror
		  if (y_pos >= 4.6 && x_pos >= 9.4) fPulseInt_quad[0][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 2 mirror
		  if (y_pos < 4.6 && x_pos >= 9.4) fPulseInt_quad[1][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 3 mirror
		  if (y_pos >= 4.6 && x_pos < 9.4) fPulseInt_quad[2][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 4 mirror
		  if (y_pos < 4.6 && x_pos < 9.4) fPulseInt_quad[3][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
		}
	    }//Marks end of electron selection condition


	  //For quadrant cut strategy with particle ID cuts. In this case pions are selected
	  if (!fTrack && fCut && fPions)
	    {
	      //Retrieve particle ID information
	      // Float_t central_p = 6.0530;
	      // Float_t p = ((P_gtr_dp[0]/100.0)*central_p) + central_p;  //

	      //Fill histogram visualizaing the pion selection

              Double_t p = *P_gtr_p;

	      fCut_everything->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);
	      fCut_enorm->Fill(*P_cal_etotnorm);

	      //Cut on Shower vs preshower is a tilted ellipse, this requires an angle of rotation (in radians), x/y center, semimajor and semiminor axis
	      Float_t piangle = 0.0;
	      Float_t pix_center = 0.3;
	      Float_t piy_center = 0.03;
	      Float_t pisemimajor_axis = 0.3;
	      Float_t pisemiminor_axis = 0.02;
	      if (pow((*P_cal_fly_earray/p - pix_center)*cos(piangle) + (*P_cal_pr_eplane/p - piy_center)*sin(piangle),2)/pow(pisemimajor_axis,2) + 
		  pow((*P_cal_fly_earray/p - pix_center)*sin(piangle) - (*P_cal_pr_eplane/p - piy_center)*cos(piangle),2)/pow(pisemiminor_axis,2) < 1)
		{
		  //Fill histogram visualizaing the pion selection
		  fCut_pion->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

		  //Fill histogram of the full PulseInt spectra for each PMT
		  fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Fill histograms of what each PMT registers from each quadrant, this requires tracking the particle from the focal plane. Each quadrant is defined from the parameter files
		  Float_t y_pos = P_tr_y[0] + P_tr_ph[0]*fhgc_zpos;                        
		  Float_t x_pos = P_tr_x[0] + P_tr_th[0]*fhgc_zpos;
		  
		
		  //Condition for quadrant 1 mirror                                                                        
		  if (y_pos >= 4.6 && x_pos >= 9.4) fPulseInt_quad[0][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 2 mirror
		  if (y_pos < 4.6 && x_pos >= 9.4) fPulseInt_quad[1][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 3 mirror
		  if (y_pos >= 4.6 && x_pos < 9.4) fPulseInt_quad[2][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Condition for quadrant 4 mirror
		  if (y_pos < 4.6 && x_pos < 9.4) fPulseInt_quad[3][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
		}
	    }   
	  //Marks end of pion selection condition
		      
	   //For quadrant cut strategy with no particle ID cut
	  if (!fTrack && !fCut)
	    {
	      //Fill histogram of the full PulseInt spectra for each PMT
	      fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

	      //Retrieve information for particle tracking from focal plane

	      //Fill histograms of what each PMT registers from each quadrant, this requires tracking the particle from the focal plane. Each quadrant is defined from the parameter files
	      Float_t y_pos = P_tr_y[0] + P_tr_ph[0]*fhgc_zpos;
	      Float_t x_pos = P_tr_x[0] + P_tr_th[0]*fhgc_zpos;
		  
	      //Condition for quadrant 1 mirror
	      if (y_pos >= 4.6 && x_pos >= 9.4) fPulseInt_quad[0][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

	      //Condition for quadrant 2 mirror
	      if (y_pos < 4.6 && x_pos >= 9.4) fPulseInt_quad[1][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

	      //Condition for quadrant 3 mirror
	      if (y_pos >= 4.6 && x_pos < 9.4) fPulseInt_quad[2][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

	      //Condition for quadrant 4 mirror
	      if (y_pos < 4.6 && x_pos < 9.4) fPulseInt_quad[3][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
	    }//Marks end of no particle ID strategy 
	  	  
	    //For TracksFired cut strategy with no particle ID cut
	  if (fTrack && !fCut)
	    {
	      //Fill histogram of the full PulseInt spectra for each PMT
	      fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

	      //Fill histograms with TracksFired cut, note that quadrant cuts are included
	      for (Int_t iregion = 0; iregion < 4; iregion++)
		{
		  if (P_hgcer_numTracksFired[iregion] == (iregion + 1)) fPulseInt_quad[iregion][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
		}
	    }//Marks end of tracksfired strategy with no particle ID

	  //For TracksFired cut strategy selecting electrons
	  if (fTrack && fCut && !fPions)
	    {
	      //Retrieve particle ID information
	      // Float_t central_p = 6.0530;
	      // Float_t p = ((P_gtr_dp[0]/100.0)*central_p) + central_p;

	      //Fill histogram visualizaing the electron selection

	      Double_t p = *P_gtr_p;

	      fCut_everything->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

	      //Cut on Shower vs preshower is a tilted ellipse, this requires an angle of rotation (in radians), x/y center, semimajor and semiminor axis
	      Float_t eangle = 3.0*3.14159/4;
	      Float_t ex_center = 0.375;
	      Float_t ey_center = 0.360;
	      Float_t esemimajor_axis = 0.38;
	      Float_t esemiminor_axis = 0.05;
	      if (pow((*P_cal_fly_earray/p - ex_center)*cos(eangle) + (*P_cal_pr_eplane/p - ey_center)*sin(eangle),2)/pow(esemimajor_axis,2) + 
		  pow((*P_cal_fly_earray/p - ex_center)*sin(eangle) - (*P_cal_pr_eplane/p - ey_center)*cos(eangle),2)/pow(esemiminor_axis,2) < 1)
		{
		  //Fill histogram visualizing the electron selection
		  fCut_electron->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

		  //Fill histogram of the full PulseInt spectra for each PMT
		  fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Fill histograms with TracksFired cut, note that quadrant cuts are included so any off quadrant histograms will be empty
		  for (Int_t iregion = 0; iregion < 4; iregion++)
		    {
		      if (P_hgcer_numTracksFired[iregion] != 0.0) fPulseInt_quad[iregion][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
		    }
		}
		}//Marks end of tracksfired with electrons

	  //For TracksFired cut strategy selecting pions
	  if (fTrack && fCut && fPions)
	    {
	      //Retrieve particle ID information
	      // Float_t central_p = 6.0530;
	      // Float_t p = ((P_gtr_dp[0]/100.0)*central_p) + central_p;

	      //Fill histogram visualizaing the electron selection

              Double_t p = *P_gtr_p;

	      fCut_everything->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

	      //Cut on Shower vs preshower is a tilted ellipse, this requires an angle of rotation (in radians), x/y center, semimajor and semiminor axis
	      Float_t piangle = 0.0;
	      Float_t pix_center = 0.26;
	      Float_t piy_center = 0.03;
	      Float_t pisemimajor_axis = 0.1;
	      Float_t pisemiminor_axis = 0.02;
	      if (pow((*P_cal_fly_earray/p - pix_center)*cos(piangle) + (*P_cal_pr_eplane/p - piy_center)*sin(piangle),2)/pow(pisemimajor_axis,2) + 
		  pow((*P_cal_fly_earray/p - pix_center)*sin(piangle) - (*P_cal_pr_eplane/p - piy_center)*cos(piangle),2)/pow(pisemiminor_axis,2) < 1)
		{
		  //Fill histogram visualizing the electron selection
		  fCut_pion->Fill(*P_cal_fly_earray/p, *P_cal_pr_eplane/p);

		  //Fill histogram of the full PulseInt spectra for each PMT
		  fPulseInt[ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);

		  //Fill histograms with TracksFired cut, note that quadrant cuts are included
		  for (Int_t iregion = 0; iregion < 4; iregion++)
		    {
		      if (P_hgcer_numTracksFired[iregion] != 0.0) fPulseInt_quad[iregion][ipmt]->Fill(P_hgcer_goodAdcPulseInt[ipmt]);
		    }
		}
		}//Marks end of tracksfired with electrons
	  
	}//Marks end of loop over PMTs

  	  
    }//Marks end of loop over tracks
  
  return kTRUE;
}

void calibration::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void calibration::Terminate()
{
 // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  printf("\n");
  Info("Terminate", "'%s' showing", (fFullShow ? "full" : "minimal"));
  Info("Terminate", "'%s' strategy", (fTrack ? "tracking" : "quadrant"));
  Info("Terminate", "cuts %s performed", (fCut ? "are" : "are not"));
  if (fCut) Info("Terminate", "cutting for '%s'", (fPions ? "pions" : "electrons"));
  printf("\n");
  Info("Terminate", "Histograms formed, now starting calibration.\n'Peak Buffer full' is a good warning!\n");
  printf("\n");

  gStyle->SetOptStat(1000000001);

//Have to extract the histograms from the OutputList
  TH1F* PulseInt[4];
  TH1F* PulseInt_quad[4][4];
  for (Int_t ipmt = 0; ipmt < 4; ipmt++)
    {
      PulseInt[ipmt] = dynamic_cast<TH1F*> (GetOutputList()->FindObject(Form("PulseInt_PMT%d",ipmt+1)));
      for (Int_t iquad = 0; iquad < 4; iquad++)
	{
	  PulseInt_quad[iquad][ipmt] = dynamic_cast<TH1F*> (GetOutputList()->FindObject(Form("PulseInt_quad%d_PMT%d",iquad+1,ipmt+1)));
	} 
    }

  //Rebin the histograms, add functionality to bin HGC & NGC independently
  if (fTrack) {
    for (Int_t ipmt=0; ipmt < (fhgc_pmts); ipmt++)
      {
	for (Int_t iquad=0; iquad<4; iquad++)
	  {
	    PulseInt_quad[iquad][ipmt]->Rebin(4);
	  }
	PulseInt[ipmt]->Rebin(4);
      }
  }

 //Canvases to display cut information
  if (fFullShow)
    {
      //Canvas to show beta cut information
      TCanvas *Beta;
      Beta = new TCanvas("Beta", "Beta information for events");
      Beta->Divide(2,1);
      Beta->cd(1);
      fBeta_Full->Draw();
      Beta->cd(2);
      fBeta_Cut->Draw();
      Beta->SaveAs("vijay1.pdf");

      //Canvas to show timing cut information
      TCanvas *Timing;
      Timing = new TCanvas("Timing", "Timing information for events");
      Timing->Divide(2,1);
      Timing->cd(1);
      fTiming_Full->Draw("Colz");
      Timing->cd(2);
      fTiming_Cut->Draw("Colz");
      Timing->SaveAs("vijay.pdf");
      TCanvas *Timing1;
      Timing1 = new TCanvas("Timing1","time info.");
      Timing1->Divide(2,2);
      Timing1->cd(1);
      fTim1->Draw();
      Timing1->cd(2);
      fTim2->Draw();
      Timing1->cd(3);
      fTim3->Draw();
      Timing1->cd(4);
      fTim4->Draw();
      
      TCanvas *pmt1_2;
      pmt1_2 = new TCanvas("pmt1_2","Bet. PMT1 &PMT2");
      pmt1_2->Divide(2,1);
      pmt1_2->cd(1);
      //  fPMT1_2->Draw("Colz");
   
    } 

if (fCut)
    {
      TCanvas *cut_enorm = new TCanvas("cut_enorm", "Visualization of etotnorm");
      fCut_enorm->Draw();
      TLatex text(1.5, 120000, "VIJAY");

      // TLatex *pt = new TLatex(0,2.4,100,120);
  
     text.DrawClone("SAME");

      TCanvas *cut_visualization = new TCanvas("cut_visualization", "Visualization of the particle ID cuts performed");
      cut_visualization->Divide(2,1);
      cut_visualization->cd(1);
      fCut_everything->Draw("Colz");
      cut_visualization->cd(2);
      fPions ? fCut_pion->Draw("Colz") : fCut_electron->Draw("Colz");
    }
   gStyle->SetOptFit(111);
//Single Gaussian to find mean of SPE
  TF1 *Gauss1 = new TF1("Gauss1",gauss,100,3,3);
  Gauss1->SetParNames("Amplitude","Mean","Std. Dev.");

  //Sum of two Gaussians to determine SPE with minimal systematics
  TF1 *Gauss2 = new TF1("Gauss2",gauss,100,6,6);
  Gauss2->SetParNames("Amplitude 1","Mean 1","Std. Dev. 1","Amplitude 2","Mean 2","Std. Dev. 2");

  //Sum of three Gaussians to determine NPE spacing
  TF1 *Gauss3 = new TF1("Gauss3",gauss,20,3.5,9);
  Gauss3->SetParNames("Amplitude 1","Mean 1","Std. Dev. 1","Amplitude 2","Mean 2","Std. Dev. 2","Amplitude 3","Mean 3","Std. Dev. 3");

  //Poisson distribution to remove high NPE background
  TF1 *Poisson = new TF1("Poisson",poisson,0.0,30,2);
  Poisson->SetParNames("Mean", "Amplitude");

  //Note about Poisson background, the mean varies between detectors/operating conditions so this quantity may require user input
  Double_t Poisson_mean;
  Poisson_mean = 5.5;  

  //Linear function used to determine goodness-of-fit for NPE spacing
  TF1 *Linear = new TF1("Linear",linear,0,4,2);
  Linear->SetParNames("Slope", "Intercept");
      
  //An array is used to store the means for the SPE, and to determine NPE spacing
  Double_t mean[3];
  Double_t SD[3];
  Double_t mean_err[3];
  Double_t x_npe[3], y_npe[3], x_err[3], y_err[3];
  Double_t RChi2[3];
  Bool_t GoodFit[6];

  //Two more arrays are used to store the estimates for the calibration constants and another two to store goodness of calibration
  Double_t calibration_mk1[4], calibration_mk1Err[4], calibration_mk2[4], calibration_mk2Err[4], pmt_calib[4], pmt_calib_mk2[4];

  TPaveText *GoodFitText = new TPaveText (0.65, 0.15, 0.85, 0.2, "NDC");
  GoodFitText->SetTextColor(kGreen);
  GoodFitText->AddText("Good fit");
  TPaveText *BadFitText = new TPaveText (0.65, 0.15, 0.85, 0.2, "NDC");  
  BadFitText->SetTextColor(kRed);
  BadFitText->AddText("Bad fit");

  TString outputpdf = "PMT_Fits.pdf";    //Name of the pdf file 
                      
  //Array to hold the Poisson character of the calibrations
  Double_t Pois_Chi[2];
  Pois_Chi[0] = 0.0, Pois_Chi[1] = 0.0;
  gStyle->SetOptFit(111);

 //Main loop for calibration
  for (Int_t ipmt=0; ipmt < (fhgc_pmts); ipmt++)
    {  
      //Initialize the various arrays (calibration arrays are explicitly filled)
      for (Int_t i=0; i<3; i++)
	{
	  mean[i] = 0.0;
	  x_npe[i] = 0, y_npe[i] = 0, x_err[i] = 0, y_err[i] = 0;
	  RChi2[i] = 0;
	}  //Begin strategy for quadrant cut calibration
      if (!fTrack)
	{
	  //TSpectrum class is used to find the SPE peak using the search method
	  TSpectrum *s = new TSpectrum(2);  

	  //Create Canvas to see the search result for the SPE  
	  if (fFullShow) quad_cuts[ipmt] = new TCanvas(Form("quad_cuts_%d",ipmt), Form("First Photoelectron peaks PMT%d",ipmt+1));
	  if (fFullShow) quad_cuts[ipmt]->Divide(3,1);  
	  
	  Int_t ipad = 1; //Variable to draw over pads correctly
	 
	  for (Int_t iquad=0; iquad<4; iquad++)
	    { 
	      if (iquad == ipmt) continue; //ignore a PMT looking at its own quadrant
	      if (fFullShow) quad_cuts[ipmt]->cd(ipad);

	      if (PulseInt_quad[iquad][ipmt]->GetEntries() > 0) 
		{
		  //Perform search for the SPE and save the peak into the array xpeaks
		  fFullShow ? s->Search(PulseInt_quad[iquad][ipmt], 2.5, "nobackground",0.001) : s->Search(PulseInt_quad[iquad][ipmt], 2.5, "nobackground&&nodraw",0.001);
                    
		  TList *functions = PulseInt_quad[iquad][ipmt]->GetListOfFunctions(); 
		  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
		  
		  if ( pm == nullptr)               
		    {
		      cout << "pm is null!!!\n\n ";                                   
		      cout << "ipmt = " << ipmt << " and iquad = " << iquad <<endl;
		      continue;
		    }
               		 		  
		  Double_t *xpeaks = pm->GetX();
		  // if (xpeaks[1] < xpeaks[0]) xpeaks[1] = xpeaks[0];
		  
		  //Use the peak to fit the SPE with a Gaussian to determine the mean
		  if(fFullShow) PulseInt_quad[iquad][ipmt]->Draw("E");
		  Gauss2->SetRange(0,17);
		  Gauss2->SetParameter(1, xpeaks[0]);
		  Gauss2->SetParameter(2, 5.0);
		  Gauss2->SetParameter(4, xpeaks[1]);	
		  Gauss2->SetParameter(5, 2.0);
		  Gauss2->SetParLimits(0, 0., PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])));
		  Gauss2->SetParLimits(1, xpeaks[0]-1, xpeaks[0]+1);
		  Gauss2->SetParLimits(2, 0.5, 10.);
		  Gauss2->SetParLimits(3, 0., PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[1])));
		  Gauss2->SetParLimits(4, xpeaks[1]-1, xpeaks[1]+1);
		  Gauss2->SetParLimits(5, 0.5, 10.);
		  fFullShow ? PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQ") : PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQN");
		  if (fFullShow) PulseInt_quad[iquad][ipmt]->GetXaxis()->SetRangeUser(0,30);
		  
		  // cout<< "  " << PulseInt_quad[iquad][ipmt]->GetBinContent(20) << "   " << PulseInt_quad[iquad][ipmt]->GetBinError(20) << endl;// mean[ipad-1] = Gauss2->GetParameter(1);

		  //if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) mean[ipad-1] = Gauss2->GetParameter(1); 
		  //if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) SD[ipad-1] = Gauss2->GetParameter(2); 
		  //if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) RChi2[ipad-1] = Gauss2->GetChisquare()/Gauss2->GetNDF();
		  mean[ipad-1] = Gauss2->GetParameter(1);
		  SD[ipad-1] = Gauss2->GetParameter(2);
		  RChi2[ipad-1] = Gauss2->GetChisquare()/Gauss2->GetNDF();
		  mean_err[ipad-1] = Gauss2->GetParError(1);
		 
		  
		  if (RChi2[ipad-1] < 0.5 || RChi2[ipad-1] > 10) {
		    GoodFit[ipad-1] = kFALSE; // Set Boolean of whether fit is good or not here
		    BadFitText->Draw("same");
		    cout << "BAD FIT RCHI2 "  << RChi2[ipad-1] << endl;
		  } 
		   else  if (RChi2[ipad-1] > 0.5 && RChi2[ipad-1] < 10){
		   GoodFit[ipad-1] = kTRUE;
		   GoodFitText->Draw("same");
		   cout << "GOOD FIT RCHI2 "  << RChi2[ipad-1] << endl;
		   }
		  
		  // cout << xpeaks[0] <<endl;
		  // cout<< " Amplitude " << PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) << endl;
		  // cout<< " SD " <<SD[ipad-1]<<endl;
		  // cout<<" mean "<<mean[ipad-1]<<endl;  
		  // cout<<"  error       "<<mean_err[ipad-1]<<endl;
		  // cout << " Chi2/DoF " << RChi2[ipad-1] << endl;
		  
		  ipad++;

		  //Again Use the peak to fit the SPE with a Gaussian to determine the mean
		  /*  Gauss2->SetRange(xpeaks[0]-10, xpeaks[0]+10);
		      Gauss2->SetParameter(1, Gauss2->GetParameter(1));
		      Gauss2->SetParameter(2, Gauss2->GetParameter(2) );
		      Gauss2->SetParameter(4, Gauss2->GetParameter(4));
		      Gauss2->SetParameter(5, Gauss2->GetParameter(5));
		      Gauss2->SetParLimits(0, 3., PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])));
		      Gauss2->SetParLimits(1, xpeaks[0]-3, xpeaks[0]+3);
		      Gauss2->SetParLimits(2, 0.5, 10.);
		      Gauss2->SetParLimits(3, 0., 500.);
		      Gauss2->SetParLimits(4, xpeaks[1]-3, xpeaks[1]+3);
		      Gauss2->SetParLimits(5, 0.5, 10.);
		      fFullShow ? PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQ") : PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQN");
		      if (fFullShow) PulseInt_quad[iquad][ipmt]->GetXaxis()->SetRangeUser(0,40);*/
		} 
	    }}}

    }
