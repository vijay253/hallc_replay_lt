
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
  ADC_max = 35;
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


  fTim1 = new TH1F("Timing_PMT1", "ADC TDC Diff PMT1 ; Time (ns) ;Counts", 200, -10.0, 50.0);
  GetOutputList()->Add(fTim1);
  fTim2 = new TH1F("Timing_PMT2", "ADC TDC Diff PMT2 ; Time (ns) ;Counts", 200, -10.0, 50.0);
  GetOutputList()->Add(fTim2);
  fTim3 = new TH1F("Timing_PMT3", "ADC TDC Diff PMT3 ; Time (ns) ;Counts", 200, -10.0, 50.0);
  GetOutputList()->Add(fTim3);
  fTim4 = new TH1F("Timing_PMT4", "ADC TDC Diff PMT4 ; Time (ns) ;Counts", 200, -10.0, 50.0);
  GetOutputList()->Add(fTim4);

  //Timing and Beta cut visualizations
  fBeta_Cut = new TH1F("Beta_Cut", "Beta cut used for 'good' hits;Beta;Counts", 100, -0.1, 1.5);
  GetOutputList()->Add(fBeta_Cut); 
  
  fBeta_Full = new TH1F("Beta_Full", "Full beta for events;Beta;Counts", 100, -0.1, 1.5);
  GetOutputList()->Add(fBeta_Full);

  // fTim = new TH1F("Tim", "Time ; Time ;Counts", 100, -50.0, 0.0);
  // GetOutputList()->Add(fTim);

  fTiming_Cut = new TH1F("Timing_Cut", "Timing cut used for 'good' hits; Time (ns);Counts", 200, -10, 50);
  GetOutputList()->Add(fTiming_Cut);

  //fTiming_Cut = new TH2F("Timing_Cut", "Timing cut used for good hits ; Time (ns); Counts", 500, -50, 50, 500, 0.0, 200);
  //GetOutputList()->Add(fTiming_Cut);

  fTiming_Full = new TH1F("Timing_Full", "Full timing information for events;Time (ns);Counts", 200, -40, 50);
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
  //if (*Ndata_P_tr_beta != 1) return kTRUE;
  
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
	  fTiming_Full->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);       

	    if(ipmt ==0){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >38 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 34) continue;                      //13 9

	    fTim1->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);

	  }

	  if(ipmt ==1){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >36 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 33) continue;                          //12 7

	    fTim2->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	  }
	  if(ipmt ==2){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >36 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 33) continue;                           //12 7

	    fTim3->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	  }
	  if(ipmt ==3){

	    if(P_hgcer_goodAdcTdcDiffTime[ipmt] >38 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 30) continue;                                  //12 8
 
	    fTim4->Fill(P_hgcer_goodAdcTdcDiffTime[ipmt]);
      
	    }
	  // cut modified by VK, 24/05/19

	  // if(P_hgcer_goodAdcTdcDiffTime[ipmt] >40 || P_hgcer_goodAdcTdcDiffTime[ipmt] < 30) continue;
	   
	   //fTiming_Cut->Fill(P_hgcer_xAtCer[ipmt],P_hgcer_yAtCer[ipmt]);

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
	      }//Marks end of no particle ID strategy */
	  	  
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
		}//Marks end of tracksfired strategy with no particle ID*/

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
		}//Marks end of tracksfired with electrons*/

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
		}//Marks end of tracksfired with electrons*/
	  
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
      Beta->SaveAs("Beta.pdf");

      //Canvas to show timing cut information
      TCanvas *Timing;
      Timing = new TCanvas("Timing", "Timing information for events");
      Timing->Divide(2,1);
      Timing->cd(1);
      fTiming_Full->Draw("Colz");
      Timing->cd(2);
      fTiming_Cut->Draw("Colz");
      Timing->SaveAs("Timing.pdf");
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
      
      /* TCanvas *pmt1_2;
      pmt1_2 = new TCanvas("pmt1_2","Bet. PMT1 &PMT2");
      pmt1_2->Divide(2,1);
      pmt1_2->cd(1);
      fPMT1_2->Draw("Colz");*/
   
    } 

  //Show the particle cuts performed in the histogram forming
  if (fCut)
    {
      TCanvas *cut_enorm = new TCanvas("cut_enorm", "Visualization of etotnorm");
      fCut_enorm->Draw();

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
  TF1 *Gauss3 = new TF1("Gauss3",gauss,20,9,9);
  Gauss3->SetParNames("Amplitude 1","Mean 1","Std. Dev. 1","Amplitude 2","Mean 2","Std. Dev. 2","Amplitude 3","Mean 3","Std. Dev. 3");

  //Poisson distribution to remove high NPE background
  TF1 *Poisson = new TF1("Poisson",poisson,0.0,5.0,2.0);
  Poisson->SetParNames("Mean", "Amplitude");

 //sum_gauss_poisson  distribution to remove high NPE background

  TF1 *sum_gauss_poisson = new TF1("sum_gauss_poisson", sum_gauss_poisson1, 0.0, 5.0,11.0);
  sum_gauss_poisson->SetParNames("Amplitude 1","Mean 1","Std. Dev. 1","Amplitude 2","Mean 2","Std. Dev. 2","Amplitude 3","Mean 3","Std. Dev. 3", "Mean 4", "Amplitude 4");


  //Note about Poisson background, the mean varies between detectors/operating conditions so this quantity may require user input
  Double_t Poisson_mean;
  Poisson_mean = 5.5;  

  //Linear function used to determine goodness-of-fit for NPE spacing
  TF1 *Linear = new TF1("Linear",linear,0,5,2);
  Linear->SetParNames("Slope", "Intercept");
      
  //An array is used to store the means for the SPE, and to determine NPE spacing
  Double_t mean[3];
  Double_t SD[3];
  Double_t mean_err[3];
  Double_t x_npe[3], y_npe[3], x_err[3], y_err[3];
  Double_t RChi2[3];
  Bool_t GoodFit[3];

  //Two more arrays are used to store the estimates for the calibration constants and another two to store goodness of calibration
  Double_t calibration_mk1[4], calibration_mk1Err[4], calibration_mk2[4], calibration_mk2Err[4], pmt_calib[4], pmt_calib_mk2[4];

  TPaveText *GoodFitText = new TPaveText (0.65, 0.15, 0.85, 0.2, "NDC");
  GoodFitText->SetTextColor(kGreen);
  GoodFitText->AddText("Good fit");
  TPaveText *BadFitText = new TPaveText (0.65, 0.15, 0.85, 0.2, "NDC");  
  BadFitText->SetTextColor(kRed);
  BadFitText->AddText("Bad fit");

  TString outputpdf = "PMT_Fits.pdf";                        //Name of the pdf file
  //Array to hold the Poisson character of the calibrations
  Double_t Pois_Chi[2];
  Pois_Chi[0] = 0.0, Pois_Chi[1] = 0.0;
  gStyle->SetOptFit(111);
  //Main loop for calibration
  for (Int_t ipmt=0; ipmt < (fhgc_pmts); ipmt++)
    {  
      //Initialize the various arrays (calibration arrays are explicitly filled)
      for (Int_t i=0; i<4; i++)
	{
	  mean[i] = 0.0;
	  x_npe[i] = 0, y_npe[i] = 0, x_err[i] = 0, y_err[i] = 0;
	  RChi2[i] = 0;
	} 
      //Begin strategy for quadrant cut calibration
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
		  //Perform search for the SPE and save the peak into the array xpeaks   //0.001
		   fFullShow ? s->Search(PulseInt_quad[iquad][ipmt], 2.0, "nobackground",0.05) : s->Search(PulseInt_quad[iquad][ipmt], 2.0, "nobackground&&nodraw",0.05);
                    
		    TList *functions = PulseInt_quad[iquad][ipmt]->GetListOfFunctions(); 
		    TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
		  
		  if ( pm == nullptr)               
		    {
		      cout << "pm is null!!!\n\n ";                                   
		      cout << "ipmt = " << ipmt << " and iquad = " << iquad <<endl;
		      continue;
		    }	 
		  	
		   Double_t k[3];
		 
		   Double_t * xpeaks = pm->GetX();   

		   cout<<" xpeaks[0] ="<<xpeaks[0]<<endl;
		   cout<<"xpeaks[1] = "<<xpeaks[1]<<endl; 
		   
		    if (xpeaks[1] < xpeaks[0]) 
 
		      { cout<<"we are here"<<endl;
                     k[0] =  xpeaks[0]; xpeaks[0] = xpeaks[1]; xpeaks[1] = k[0];}
		    cout<<" xpeaks[0] ="<<xpeaks[0]<<endl;
		    cout<<"xpeaks[1] = "<<xpeaks[1]<<endl;
		    
		   
		  //Use the peak to fit the SPE with a Gaussian to determine the mean
                  if(fFullShow) PulseInt_quad[iquad][ipmt]->Draw("E");
		  Gauss2->SetRange(0,18);
                  Gauss2->SetParameter(0, 200);
		  Gauss2->SetParameter(1, xpeaks[0]);
		  Gauss2->SetParameter(2, 5.25);
                  Gauss2->SetParameter(3, 40);
		  Gauss2->SetParameter(4, xpeaks[1]);	
		  Gauss2->SetParameter(5, 2.25);
		  Gauss2->SetParLimits(0, 0., PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])));
		  Gauss2->SetParLimits(1, xpeaks[0]-2, xpeaks[0]+2);
		  Gauss2->SetParLimits(2, 0.5, 10.0);
		  Gauss2->SetParLimits(3, 0., PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[1])));
		  Gauss2->SetParLimits(4, xpeaks[1]-0.5, xpeaks[1]+0.5);
		  Gauss2->SetParLimits(5, 0.5, 4.0);
		  fFullShow ? PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQN") : PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQN");		
         	  //  if (fFullShow) PulseInt_quad[iquad][ipmt]->GetXaxis()->SetRangeUser(0,50);
	
                 //Again Use the peak to fit the SPE with a Gaussian to determine the mean
                  Gauss2->SetParameter(0,Gauss2->GetParameter(0));
		  Gauss2->SetParameter(1,Gauss2->GetParameter(1));
		  Gauss2->SetParameter(2,Gauss2->GetParameter(2));
                  Gauss2->SetParameter(3,Gauss2->GetParameter(3));
		  Gauss2->SetParameter(4,Gauss2->GetParameter(4));	
		  Gauss2->SetParameter(5,Gauss2->GetParameter(5));
		  // Gauss2->SetLineColor(2);
	          Gauss2->SetLineWidth(2);
	       	  fFullShow ? PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQ"): PulseInt_quad[iquad][ipmt]->Fit("Gauss2","RQ");

		  // Draw the Gauss fuctions on first & second peaks using parameters from Gauss2 fit  
		  TF1 *g1 = new TF1("g1","gaus",0,35);
		  if (xpeaks[1] < xpeaks[0]) {
                  g1->SetParameter(0,Gauss2->GetParameter(3));
		  g1->SetParameter(1,Gauss2->GetParameter(4));
		  g1->SetParameter(2,Gauss2->GetParameter(5));
                  g1->SetLineColor(1);
                  g1->Draw("same");}
		  // cout<<"g1  = "<<g1->GetParameter(1);


		  else {g1->SetParameter(0,Gauss2->GetParameter(0));
		  g1->SetParameter(1,Gauss2->GetParameter(1));
		  g1->SetParameter(2,Gauss2->GetParameter(2));
                  g1->SetLineColor(1);
                  g1->Draw("same");}

	       	  TF1 *g2 = new TF1("g2","gaus",0,35);
                  if (xpeaks[1] < xpeaks[0]) {
                  g2->SetParameter(0,Gauss2->GetParameter(0));
		  g2->SetParameter(1,Gauss2->GetParameter(1));	
		  g2->SetParameter(2,Gauss2->GetParameter(2));
                  g2->SetParLimits(2, 0.5, 10.0);
                  g2->SetLineColor(3);	      	       
		  g2->Draw("same");}
    
		  else {
                  g2->SetParameter(0,Gauss2->GetParameter(3));
		  g2->SetParameter(1,Gauss2->GetParameter(4));	
		  g2->SetParameter(2,Gauss2->GetParameter(5));
                  g2->SetParLimits(2, 0.5, 10.0);
                  g2->SetLineColor(3);	      	       
		  g2->Draw("same");} 

		  if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) mean[ipad-1] = Gauss2->GetParameter(1); 
		  if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) SD[ipad-1] = Gauss2->GetParameter(2); 
		  if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) RChi2[ipad-1] = Gauss2->GetChisquare()/Gauss2->GetNDF(); 
		  if (xpeaks[0] > 2.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 90) mean_err[ipad-1] = Gauss2->GetParError(1);

		   // Set Boolean of whether fit is good or not here
		  
		   if (RChi2[ipad-1] < 0.5 || RChi2[ipad-1] > 5) {
		    GoodFit[ipad-1] = kFALSE; 
		    BadFitText->Draw("same");
		          // cout << "BAD FIT RCHI2 "  << RChi2[ipad-1] << endl;
		  } 
		    else if  (RChi2[ipad-1] > 0.5 && RChi2[ipad-1] < 5) {
		         GoodFit[ipad-1] = kTRUE;
		         GoodFitText->Draw("same");
			 //cout << "GOOD FIT RCHI2 "  << RChi2[ipad-1] << endl;
		    } 
		  
		  // cout << xpeaks[0] <<endl;
		  // cout<< " Amplitude " << PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) << endl;
		  // cout<< " SD " <<SD[ipad-1]<<endl;
		  // cout<<" mean "<<mean[ipad-1]<<endl;  
		  // cout<<"  error       "<<mean_err[ipad-1]<<endl;
		  // cout << " Chi2/DoF " << RChi2[ipad-1] << endl;
		  
		   ipad++;

                 } 
		}

	  
	  //Obtain the conversion from ADC to NPE by taking the average of the SPE means
	  Double_t xscale = 0.0;
	  Double_t xscaleErr = 0.0;
	  Double_t num_peaks = 0.0;
	  Double_t WeightAvgSum1 = 0.0;
	  Double_t WeightAvgSum2 = 0.0;
	  for (Int_t i=0; i<3; i++)
	    {
	      // if (mean[i] == 0.0) continue;
	      //xscale += mean[i]; // Old version of getting xscale which was simple mean of means
	      //num_peaks += 1.0;
	      if (GoodFit[i] == kFALSE) continue;                     // Take weighted avg 26/8/19 SK
	      if (mean_err[i] == 0) continue;                        
	      WeightAvgSum1 += mean[i]/(mean_err[i]*mean_err[i]);                
	      WeightAvgSum2 += 1/(mean_err[i]*mean_err[i]);
	    }
	  //cout << "Weighted Average = " << WeightAvgSum1/WeightAvgSum2 << " pm " << 1/TMath::Sqrt(WeightAvgSum2) << endl;
	  //if (num_peaks != 0.0) xscale = xscale/num_peaks;
	  xscale = WeightAvgSum1/WeightAvgSum2;

	  // cout<<" xscale =" <<xscale<<endl;
	  xscaleErr = 1/TMath::Sqrt(WeightAvgSum2);    

	  //Perform check if the statistics were too low to get a good estimate of the SPE mean
	  if (xscale < 1.0)
	    {
	      //Repeat the exact same procedure for the SPE of each quadrant, except now its for the full PMT spectra
	      if(fFullShow) low_stats_ipmt = new TCanvas(Form("low_stats_%d",ipmt),Form("Low stats analysis for PMT%d",ipmt+1));
	      if(fFullShow) low_stats_ipmt->cd(1);
	      PulseInt[ipmt]->GetXaxis()->SetRangeUser(0,20);
	      fFullShow ? s->Search(PulseInt[ipmt], 3.5, "nobackground", 0.001) : s->Search(PulseInt[ipmt], 2.0, "nobackground&&nodraw", 0.001);
	      TList *functions = PulseInt[ipmt]->GetListOfFunctions();
	      TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	      Double_t *xpeaks = pm->GetX();
	      Gauss1->SetRange(xpeaks[0]-1, xpeaks[0]+1);
	      Gauss1->SetParameter(1, xpeaks[0]);
	      Gauss1->SetParameter(2, 3);
	      Gauss1->SetParLimits(0, 0., 3000.);
	      //Gauss1->SetParLimits(1, xpeaks[0]-1, xpeaks[0]+1);
	      //Gauss1->SetParLimits(2, 0.5, 20.0);
	      PulseInt[ipmt]->GetXaxis()->SetRangeUser(0,200);
	      fFullShow ? PulseInt[ipmt]->Fit("Gauss1","RQ") : PulseInt[ipmt]->Fit("Gauss1","RQ");
	      xscale = Gauss1->GetParameter(1);
	      xscaleErr = Gauss1->GetParError(1);
	      }	  


	  //Scale full ADC spectra according to the mean of the SPE. This requires filling a new histogram with the same number of bins but scaled min/max
	  Int_t nbins;
	  nbins = (PulseInt[ipmt]->GetXaxis()->GetNbins());

	  //cout<<"nbins = "<<nbins<<endl;

	  //With the scale of ADC to NPE create a histogram that has the conversion applied
	  fscaled[ipmt] = new TH1F(Form("fscaled_PMT%d", ipmt+1), Form("Scaled ADC spectra for PMT%d; NPE; Normalized Counts",ipmt+1), nbins,0.0,5.0);
	  
	  //Fill this histogram bin by bin
	  for (Int_t ibin=0; ibin<nbins; ibin++)
	    {
	      Double_t y = PulseInt[ipmt]->GetBinContent(ibin);
	      Double_t x = PulseInt[ipmt]->GetXaxis()->GetBinCenter(ibin);
	      Double_t x_scaled = x/xscale;
	      Int_t bin_scaled = fscaled[ipmt]->GetXaxis()->FindBin(x_scaled); 
	      fscaled[ipmt]->SetBinContent(bin_scaled,y);
	    }

	  //Normalize the histogram for ease of fitting
	   fscaled[ipmt]->Scale(1.0/fscaled[ipmt]->Integral(), "width");               

	  //Begin the removal of the Poisson-like background
	  if (fFullShow) background_ipmt = new TCanvas(Form("backgrounf_pmt%d",ipmt), Form("NPE spectra for PMT%d with Poisson-like background",ipmt+1));
	    if (fFullShow) background_ipmt->cd(1);
	    // Poisson->SetParameter(0, 4);
	    //  Poisson->SetParameter(1, 0.85);
	     // Poisson->SetParLimits(0, 0.0, 4.0);
	     // Poisson->SetParLimits(1, 0.8, 0.9);
	     //fFullShow ? fscaled[ipmt]->Fit("Poisson","RQN") : fscaled[ipmt]->Fit("Poisson","RQN");
	     // Poisson->SetParameter(0, Poisson->GetParameter(0));
	     // Poisson->SetParameter(1, Poisson->GetParameter(1));
	     // Poisson->SetParLimits(0, 0.0, 4.0);
	     //Poisson->SetParLimits(1, 0.1, 1.2);
	   
	    //  fFullShow ? fscaled[ipmt]->Fit("Poisson","RQN") : fscaled[ipmt]->Fit("Poisson","RQN");
	    //  TF1 *g3 = new TF1("g3","Gauss3",0,2);
	          sum_gauss_poisson->SetRange(0,5.0);
                  sum_gauss_poisson->SetParameter(0,0.6);
		  sum_gauss_poisson->SetParameter(1,1.0);	
		  sum_gauss_poisson->SetParameter(2,0.39);
		  sum_gauss_poisson->SetParameter(3,0.15);
		  sum_gauss_poisson->SetParameter(4,2);
		  sum_gauss_poisson->SetParameter(5,0.24);
		  sum_gauss_poisson->SetParameter(6,0.12);
		  sum_gauss_poisson->SetParameter(7,3.0);
		  sum_gauss_poisson->SetParameter(8,0.075);
		  sum_gauss_poisson->SetParameter(9,4.0);
		  sum_gauss_poisson->SetParameter(10,0.15);

		  // cout<<" error in scale "<<xscaleErr<<endl;
		  // sum_gauss_poisson->SetParLimits(0, 0 ,1.2);
                  sum_gauss_poisson->SetParLimits(1, 1 - 2*xscaleErr, 1 + 2*xscaleErr);
		  // sum_gauss_poisson->SetParLimits(2, 0.1, 0.3);
                  //sum_gauss_poisson->SetParLimits(3, 0, 0.4);
                  sum_gauss_poisson->SetParLimits(4, 2 - 3*xscaleErr,2 +  3*xscaleErr);
		  // sum_gauss_poisson->SetParLimits(5, 0, 3);
		  // sum_gauss_poisson->SetParLimits(6, 0.0, 0.24);
                  sum_gauss_poisson->SetParLimits(7, 3 - 3*xscaleErr, 3 + 3*xscaleErr);
		  // sum_gauss_poisson->SetParLimits(8, 0.0, 3);
		  // sum_gauss_poisson->SetParLimits(9, 3.9, 4.2);
		  // sum_gauss_poisson->SetParLimits(10, 0.0,0.8);*/
	       	  // sum_gauss_poisson->SetLineColor(2);	 
	     fFullShow ? fscaled[ipmt]->Fit("sum_gauss_poisson","RQN") : fscaled[ipmt]->Fit("sum_gauss_poisson","RQN");
	    /*  cout<< "Gauss3 (0) = "<<Gauss3->GetParameter(0)<<endl;
	    cout<< "Gauss3 (1) = "<<Gauss3->GetParameter(1)<<endl;
	    cout<< "Gauss3 (2) = "<<Gauss3->GetParameter(2)<<endl;
	    cout<< "Gauss3 (3) = "<<Gauss3->GetParameter(3)<<endl;
	    cout<< "Gauss3 (4) = "<<Gauss3->GetParameter(4)<<endl;
	    cout<< "Gauss3 (5) = "<<Gauss3->GetParameter(5)<<endl;
	    cout<< "Gauss3 (6) = "<<Gauss3->GetParameter(6)<<endl;
	    cout<< "Gauss3 (7) = "<<Gauss3->GetParameter(7)<<endl;
	    cout<< "Gauss3 (8) = "<<Gauss3->GetParameter(8)<<endl;*/

	    sum_gauss_poisson->SetParameter(0,sum_gauss_poisson->GetParameter(0));
	    sum_gauss_poisson->SetParameter(1,sum_gauss_poisson->GetParameter(1));	
	    sum_gauss_poisson->SetParameter(2,sum_gauss_poisson->GetParameter(2));
	    sum_gauss_poisson->SetParameter(3,sum_gauss_poisson->GetParameter(3));
	    sum_gauss_poisson->SetParameter(4,sum_gauss_poisson->GetParameter(4));
	    sum_gauss_poisson->SetParameter(5,sum_gauss_poisson->GetParameter(5));
	    sum_gauss_poisson->SetParameter(6,sum_gauss_poisson->GetParameter(6));
	    sum_gauss_poisson->SetParameter(7,sum_gauss_poisson->GetParameter(7));
	    sum_gauss_poisson->SetParameter(8,sum_gauss_poisson->GetParameter(8));
	    sum_gauss_poisson->SetParameter(9,sum_gauss_poisson->GetParameter(9));
	    sum_gauss_poisson->SetParameter(10,sum_gauss_poisson->GetParameter(10));
	    sum_gauss_poisson->SetLineColor(5);

	      fFullShow ? fscaled[ipmt]->Fit("sum_gauss_poisson","RQ") : fscaled[ipmt]->Fit("sum_gauss_poisson","RQ");

            TF1 *g3 = new TF1("g3","gaus",0,5.0);
	    {  
                  g3->SetParameter(0,sum_gauss_poisson->GetParameter(0));
		  g3->SetParameter(1,sum_gauss_poisson->GetParameter(1));
		  g3->SetParameter(2,sum_gauss_poisson->GetParameter(2));
                  g3->SetLineColor(1);
                  g3->Draw("same");}
            TF1 *g4 = new TF1("g4","gaus",0,5.0);
	    {
                  g4->SetParameter(0,sum_gauss_poisson->GetParameter(3));
		  g4->SetParameter(1,sum_gauss_poisson->GetParameter(4));
		  g4->SetParameter(2,sum_gauss_poisson->GetParameter(5));
                  g4->SetLineColor(3);
                  g4->Draw("same");}
            TF1 *g5 = new TF1("g5","gaus",0,5.0);
	    {
                  g5->SetParameter(0,sum_gauss_poisson->GetParameter(6));
		  g5->SetParameter(1,sum_gauss_poisson->GetParameter(7));
		  g5->SetParameter(2,sum_gauss_poisson->GetParameter(8));
                  g5->SetLineColor(4);
                  g5->Draw("same");}

	    /*  Poisson->SetRange(3.1, 5);

            Poisson->SetParameter(0, sum_gauss_poisson->GetParameter(9));
	    Poisson->SetParameter(1, sum_gauss_poisson->GetParameter(10));
	    // Poisson->SetParLimits(0, 0.0, 4.0);
	    Poisson->SetParLimits(1, 0.0,0.9);
	    Poisson->Draw("same");*/

            TF1 *poiss = new TF1("poiss","[1]*TMath::Poisson(x,[0])",0.0,5.0);{
            poiss->SetParameter(0,sum_gauss_poisson->GetParameter(9));
            poiss->SetParameter(1,sum_gauss_poisson->GetParameter(10));
	    poiss->Draw("same"); }

	    /* {
                  g6->SetParameter(0,sum_gauss_poisson->GetParameter(9));
		  g6->SetParameter(1,sum_gauss_poisson->GetParameter(10));	       
                  g6->SetLineColor(5);
                  g6->Draw("same");}*/


	  //Make and fill histogram with the background removed
	    /*fscaled_nobackground[ipmt] = new TH1F(Form("fscaled_nobackground_pmt%d", ipmt+1), Form("NPE spectra background removed for PMT%d; NPE; Normalized Counts",ipmt+1), 200, 0, 15);

	  for (Int_t ibin=1; ibin<nbins; ibin++)
	    {
	      Double_t y = Poisson->Eval(fscaled[ipmt]->GetXaxis()->GetBinCenter(ibin));
	      Double_t y2 = fscaled[ipmt]->GetBinContent(ibin) - y;
	      //if (ipmt == 1) cout << "Poisson Value: " << y << "\nfscaled Value: " << fscaled[ipmt]->GetBinContent(ibin) << endl;
	      fscaled_nobackground[ipmt]->SetBinContent(ibin,y2);
	      }*/

	  if (fFullShow) final_spectra_ipmt = new TCanvas(Form("final_Spectra_%d",ipmt), Form("NPE spectra Background Removed for PMT%d",ipmt+1));
	  // if (fFullShow) final_spectra_ipmt->Divide(2,1);
	  //if (fFullShow) final_spectra_ipmt->cd(1);
          /*Gauss3->SetRange(0,3.5);
	  Gauss3->SetParameter(0,0.5);
	  Gauss3->SetParameter(1,1.5);
	  Gauss3->SetParameter(2,0.2);
	  Gauss3->SetParameter(3,0.2); 
	  Gauss3->SetParameter(4,2.0);
          Gauss3->SetParameter(5,0.3);
          Gauss3->SetParameter(6,0.05);
          Gauss3->SetParameter(7,3.0);
          Gauss3->SetParameter(8,0.05);
	  Gauss3->SetParLimits(1, 0.5, 1.5);
	  Gauss3->SetParLimits(2, 0.0, 2.5);
	  Gauss3->SetParLimits(3, 0.0, 0.5);
	  Gauss3->SetParLimits(4, 1.5, 2.5);
	  Gauss3->SetParLimits(5, 0.2, 0.6);
	  Gauss3->SetParLimits(6, 0.0, 0.2);
	  Gauss3->SetParLimits(7, 2.5, 3.5);
	  Gauss3->SetParLimits(8, 0.02, 0.5);
	  fFullShow ? fscaled_nobackground[ipmt]->Fit("Gauss3","RQ") : fscaled_nobackground[ipmt]->Fit("Gauss3","RQ");
	  if (fFullShow) fscaled_nobackground[ipmt]->GetXaxis()->SetRangeUser(0,4);
	  if (fFullShow) fscaled_nobackground[ipmt]->GetYaxis()->SetRangeUser(0,0.6);*/
	 
	  //Create a TGraphErrors to determine the spacing of the NPE
          gROOT->SetStyle("111111");

	  y_npe[0] = sum_gauss_poisson->GetParameter(1), y_npe[1] =  sum_gauss_poisson->GetParameter(4), y_npe[2] = sum_gauss_poisson->GetParameter(7);
	  y_err[0] =  sum_gauss_poisson->GetParError(1), y_err[1] = sum_gauss_poisson->GetParError(4), y_err[2] = sum_gauss_poisson->GetParError(7);
	  x_npe[0] = 1, x_npe[1] = 2, x_npe[2] = 3;
	  TGraphErrors *gr_npe = new TGraphErrors(3, x_npe, y_npe, x_err, y_err);
          gr_npe->GetXaxis()->SetLimits(0, 5);
          gr_npe->SetMinimum(0);
          gr_npe->SetMaximum(5);
	  gr_npe->SetTitle(Form("Linear Spacing of PE for PMT%d",ipmt+1));
	  gr_npe->GetXaxis()->SetTitle("NPE number");
	  gr_npe->GetYaxis()->SetTitle("Photoelectron peak (NPE)");

	  //Plot this graph with the NPE spectra
	  // if (fFullShow) final_spectra_ipmt->cd(1);
	  Linear->SetParameter(0.0, 1.0);
          Linear->SetParameter(1.0, 0.0);

	  if (fFullShow) gr_npe->Draw("A*"); 

	  //  fFullShow ? gr_npe->Fit("Linear","RQ+") : gr_npe->Fit("Linear","RQ+");
     
	  //Plot this graph with the NPE spectra
	  // if (fFullShow) final_spectra_ipmt->cd(1);
	  // Linear->SetParameter(0.0, 1.0);
	  // Linear->SetParameter(1.0, 0.0);

	  // To check the linearfit with linear function making intercept zero

	  TF1 *gl = new TF1("gl","[0]*x", 0.0,5.0);
	  //  gl->SetRange(-2.0,4.0);
 
	    gl->SetParNames("Slope");

	    //  gl->SetParameters(1.0 ,0.0 );
	    gl->SetParameter(0.0,1.0);
	    //  gl->SetParameter(1,0.0);
	    gl->SetLineColor(4);

	  fFullShow ? gr_npe->Fit("Linear","RQ") : gr_npe->Fit("Linear","RQ");
	  fFullShow ? gr_npe->Fit("gl","RQ+") : gr_npe->Fit("gl","RQ+");

	  Double_t sl, sl_err, ChiNDF;
          ChiNDF = gl->GetChisquare()/gl->GetNDF(); 
	  sl = gl->GetParameter(0);
	  sl_err = gl->GetParError(0);
	  TPaveText *t = new TPaveText(0.60, 0.62, 0.90, 0.75, "NDC");{
           t->SetTextColor(kBlack);
	   t->AddText(Form(" Chi/NDF = %0f", ChiNDF ));
	   t->AddText(Form(" Slope = %0f", sl ));
	   t->AddText(Form(" Slope_Error = %0f", sl_err ));
	   t->AddText(Form(" Intercept = 0 " ));
           t->Draw();}
	     

 	  calibration_mk1[ipmt] = xscale;
	  calibration_mk1Err[ipmt] = xscaleErr;
	  pmt_calib[ipmt] = abs(1.0 -sum_gauss_poisson->GetParameter(1));

	   //Initial calibration constant has been obtained. Now I multiply it by the slope of the spacing of the NPE (should be approx. 1) for a second estimate

	  Double_t xscale_mk2 = xscale * Linear->GetParameter(0);  

	  //  cout<<" slope of linear fit"<<Linear->GetParameter(0)<<endl;     
	  // Double_t xscale_mk2Err = Sqrt(Power(xscaleErr*Linear->GetParameter(0),2) +  Power(xscale*Linear->GetParError(0),2)); 
         Double_t xscale_mk2Err = (Sqrt(Power(xscaleErr,2) +  Power(Linear->GetParError(0),2))); 
	  //Take this new xscale and repeat the exact same procedure as before
	  fscaled_mk2[ipmt] = new TH1F(Form("fhgc_scaled_mk2_PMT%d", ipmt+1), Form("Scaled ADC spectra for PMT%d; NPE; Normalized Counts",ipmt+1), 200, 0, 5.0);
	 
	  //Fill this histogram bin by bin
	  for (Int_t ibin=0; ibin<nbins; ibin++)
	    {
	      Double_t y = PulseInt[ipmt]->GetBinContent(ibin);
	      Double_t x = PulseInt[ipmt]->GetXaxis()->GetBinCenter(ibin);
	      Double_t x_scaled = x/xscale_mk2;
	      Int_t bin_scaled = fscaled_mk2[ipmt]->GetXaxis()->FindBin(x_scaled); 
	      fscaled_mk2[ipmt]->SetBinContent(bin_scaled,y);
	    }
	 
	  //Normalize the histogram for ease of fitting
	  fscaled_mk2[ipmt]->Scale(1.0/fscaled_mk2[ipmt]->Integral(), "width");

	  //Begin the removal of the Poisson-like background
	  if (fFullShow) background_mk2_ipmt = new TCanvas(Form("background_mk2_pmt%d",ipmt), Form("NPE spectra for PMT%d with Poisson-like background",ipmt+1));
	  if (fFullShow) background_mk2_ipmt->cd(1);
	  /*  Poisson->SetParameter(0, 2);
	  Poisson->SetParameter(1, 2.0);
	  Poisson->SetParLimits(0, 4.0, 1.0);
	  fFullShow ? fscaled_mk2[ipmt]->Fit("Poisson","RQ"):fscaled_mk2[ipmt]->Fit("Poisson","RQ");*/
 if (fFullShow) background_ipmt->cd(1);
 // Poisson->SetParameter(0, 4);
 // Poisson->SetParameter(1, 0.85);
	     // Poisson->SetParLimits(0, 0.0, 4.0);
	     // Poisson->SetParLimits(1, 0.8, 0.9);
	     //fFullShow ? fscaled[ipmt]->Fit("Poisson","RQN") : fscaled[ipmt]->Fit("Poisson","RQN");
	     // Poisson->SetParameter(0, Poisson->GetParameter(0));
	     // Poisson->SetParameter(1, Poisson->GetParameter(1));
	     // Poisson->SetParLimits(0, 0.0, 4.0);
	     //Poisson->SetParLimits(1, 0.1, 1.2);
	   
	    //  fFullShow ? fscaled[ipmt]->Fit("Poisson","RQN") : fscaled[ipmt]->Fit("Poisson","RQN");
	    //  TF1 *g3 = new TF1("g3","Gauss3",0,2);
	          sum_gauss_poisson->SetRange(0,5.0);
                  sum_gauss_poisson->SetParameter(0,0.6);
		  sum_gauss_poisson->SetParameter(1,0.95);	
		  sum_gauss_poisson->SetParameter(2,0.39);
		  sum_gauss_poisson->SetParameter(3,0.2);
		  sum_gauss_poisson->SetParameter(4,2);
		  sum_gauss_poisson->SetParameter(5,0.24);
		  sum_gauss_poisson->SetParameter(6,0.12);
		  sum_gauss_poisson->SetParameter(7,3.0);
		  sum_gauss_poisson->SetParameter(8,0.075);
		  sum_gauss_poisson->SetParameter(9,4.0);
		  sum_gauss_poisson->SetParameter(10,0.15);

		  /* sum_gauss_poisson->SetParLimits(0, 0 ,1.2);
                  sum_gauss_poisson->SetParLimits(1, 0.8,1.1 );
                  sum_gauss_poisson->SetParLimits(2, 0.1, 0.3);
                  sum_gauss_poisson->SetParLimits(3, 0, 0.4);
                  sum_gauss_poisson->SetParLimits(4, 1.8,2.2);
                  sum_gauss_poisson->SetParLimits(5, 0, 3);
                  sum_gauss_poisson->SetParLimits(6, 0.0, 0.24);
                  sum_gauss_poisson->SetParLimits(7, 2.8, 3.2);
                  sum_gauss_poisson->SetParLimits(8, 0.0, 3);
		  sum_gauss_poisson->SetParLimits(9, 3.9, 4.2);
		  sum_gauss_poisson->SetParLimits(10, 0.0,0.8);*/
	       	  // sum_gauss_poisson->SetLineColor(2);	 
	     fFullShow ? fscaled_mk2[ipmt]->Fit("sum_gauss_poisson","RQN") : fscaled_mk2[ipmt]->Fit("sum_gauss_poisson","RQN");
	    /*  cout<< "Gauss3 (0) = "<<Gauss3->GetParameter(0)<<endl;
	    cout<< "Gauss3 (1) = "<<Gauss3->GetParameter(1)<<endl;
	    cout<< "Gauss3 (2) = "<<Gauss3->GetParameter(2)<<endl;
	    cout<< "Gauss3 (3) = "<<Gauss3->GetParameter(3)<<endl;
	    cout<< "Gauss3 (4) = "<<Gauss3->GetParameter(4)<<endl;
	    cout<< "Gauss3 (5) = "<<Gauss3->GetParameter(5)<<endl;
	    cout<< "Gauss3 (6) = "<<Gauss3->GetParameter(6)<<endl;
	    cout<< "Gauss3 (7) = "<<Gauss3->GetParameter(7)<<endl;
	    cout<< "Gauss3 (8) = "<<Gauss3->GetParameter(8)<<endl;*/

	    sum_gauss_poisson->SetParameter(0,sum_gauss_poisson->GetParameter(0));
	    sum_gauss_poisson->SetParameter(1,sum_gauss_poisson->GetParameter(1));	
	    sum_gauss_poisson->SetParameter(2,sum_gauss_poisson->GetParameter(2));
	    sum_gauss_poisson->SetParameter(3,sum_gauss_poisson->GetParameter(3));
	    sum_gauss_poisson->SetParameter(4,sum_gauss_poisson->GetParameter(4));
	    sum_gauss_poisson->SetParameter(5,sum_gauss_poisson->GetParameter(5));
	    sum_gauss_poisson->SetParameter(6,sum_gauss_poisson->GetParameter(6));
	    sum_gauss_poisson->SetParameter(7,sum_gauss_poisson->GetParameter(7));
	    sum_gauss_poisson->SetParameter(8,sum_gauss_poisson->GetParameter(8));
	    sum_gauss_poisson->SetParameter(9,sum_gauss_poisson->GetParameter(9));
	    sum_gauss_poisson->SetParameter(10,sum_gauss_poisson->GetParameter(10));
	    sum_gauss_poisson->SetLineColor(5);

	      fFullShow ? fscaled_mk2[ipmt]->Fit("sum_gauss_poisson","RQ") : fscaled_mk2[ipmt]->Fit("sum_gauss_poisson","RQ");

            TF1 *g6 = new TF1("g6","gaus",0,5.0);
	    {  
                  g6->SetParameter(0,sum_gauss_poisson->GetParameter(0));
		  g6->SetParameter(1,sum_gauss_poisson->GetParameter(1));
		  g6->SetParameter(2,sum_gauss_poisson->GetParameter(2));
                  g6->SetLineColor(1);
                  g6->Draw("same");}
            TF1 *g7 = new TF1("g7","gaus",0,5.0);
	    {
                  g7->SetParameter(0,sum_gauss_poisson->GetParameter(3));
		  g7->SetParameter(1,sum_gauss_poisson->GetParameter(4));
		  g7->SetParameter(2,sum_gauss_poisson->GetParameter(5));
                  g7->SetLineColor(3);
                  g7->Draw("same");}
            TF1 *g8 = new TF1("g8","gaus",0,5.0);
	    {
                  g8->SetParameter(0,sum_gauss_poisson->GetParameter(6));
		  g8->SetParameter(1,sum_gauss_poisson->GetParameter(7));
		  g8->SetParameter(2,sum_gauss_poisson->GetParameter(8));
                  g8->SetLineColor(4);
                  g8->Draw("same");}

	    /*  Poisson->SetRange(3.1, 5);

            Poisson->SetParameter(0, sum_gauss_poisson->GetParameter(9));
	    Poisson->SetParameter(1, sum_gauss_poisson->GetParameter(10));
	    // Poisson->SetParLimits(0, 0.0, 4.0);
	    Poisson->SetParLimits(1, 0.0,0.9);
	    Poisson->Draw("same");*/

            TF1 *poiss1 = new TF1("poiss1","[1]*TMath::Poisson(x,[0])",0.0,5.0);{
            poiss1->SetParameter(0,sum_gauss_poisson->GetParameter(9));
            poiss1->SetParameter(1,sum_gauss_poisson->GetParameter(10));
	    poiss1->Draw("same"); }

	  //Make and fill histogram with the background removed
	    /* fscaled_mk2_nobackground[ipmt] = new TH1F(Form("fscaled_mk2_nobackground_pmt%d", ipmt+1), Form("NPE spectra background removed for PMT%d; NPE; Normalized Counts",ipmt+1), 200, 0, 30);

	  for (Int_t ibin=0; ibin<nbins; ibin++)
	    {
	      Double_t y = Poisson->Eval(fscaled_mk2[ipmt]->GetXaxis()->GetBinCenter(ibin));
	      Double_t y2 = fscaled_mk2[ipmt]->GetBinContent(ibin) - y;
	      fscaled_mk2_nobackground[ipmt]->SetBinContent(ibin,y2);
	      }*/

	  if (fFullShow) final_spectra_mk2_ipmt = new TCanvas(Form("final_Spectra_mk2_%d",ipmt), Form("NPE spectra Background Removed for PMT%d",ipmt+1));
	  // if (fFullShow) final_spectra_mk2_ipmt->Divide(2,1);
	  // if (fFullShow) final_spectra_mk2_ipmt->cd(1);
	  /* Gauss3->SetParameters(0.08, 1.0, 0.22, 0.029, 2, 0.5, 0.15, 3, 0.26);
	  Gauss3->SetParLimits(1, 0.5, 1.5);
	  Gauss3->SetParLimits(2, 0.0, 1.0);
	  Gauss3->SetParLimits(3, 0.0, 1.0);
	  Gauss3->SetParLimits(4, 1.5, 2.5);
	  Gauss3->SetParLimits(5, 0.2, 0.6);
	  Gauss3->SetParLimits(6, 0.0, 1.0);
	  Gauss3->SetParLimits(7, 2.5, 3.5);
	  Gauss3->SetParLimits(8, 0.2, 0.5);
	  fFullShow ? fscaled_mk2_nobackground[ipmt]->Fit("Gauss3","RQ") : fscaled_mk2_nobackground[ipmt]->Fit("Gauss3","RQN");
	  if (fFullShow) fscaled_mk2_nobackground[ipmt]->GetXaxis()->SetRangeUser(0,5);
	  if (fFullShow) fscaled_mk2_nobackground[ipmt]->GetYaxis()->SetRangeUser(0,0.7);*/

	  //Create a TGraphErrors to determine the spacing of the NPE
	  y_npe[0] = sum_gauss_poisson->GetParameter(1), y_npe[1] = sum_gauss_poisson->GetParameter(4), y_npe[2] = sum_gauss_poisson->GetParameter(7);
	  y_err[0] = sum_gauss_poisson->GetParError(1), y_err[1] = sum_gauss_poisson->GetParError(4), y_err[2] = sum_gauss_poisson->GetParError(7);
	  x_npe[0] = 1, x_npe[1] = 2, x_npe[2] = 3;
	  TGraphErrors *gr_npe_mk2 = new TGraphErrors(3, x_npe, y_npe, x_err, y_err);
	  gr_npe_mk2->GetXaxis()->SetTitle("NPE number");
	  gr_npe_mk2->GetYaxis()->SetTitle("Photoelectron peak (NPE)");

	  //Plot this graph with the NPE spectra
	  if (fFullShow) final_spectra_mk2_ipmt->cd(1);
	  Linear->SetParameters(1.0, 0.0);
	  fFullShow ? gr_npe_mk2->Fit("Linear","RQ") : gr_npe_mk2->Fit("Linear","RQ");
	  if (fFullShow) gr_npe_mk2->Draw("A*");
	  calibration_mk2[ipmt] = xscale_mk2;
	  calibration_mk2Err[ipmt] = xscale_mk2Err;
	  pmt_calib_mk2[ipmt] = abs(1.0 - sum_gauss_poisson->GetParameter(1));
				   
	    } // This brace marks the end of the quadrant cut strategy
               

       //Begin the TrackFired cut calibration
      if (fTrack)
	{
	  //TSpectrum class is used to find the SPE peak using the search method
	  TSpectrum *s = new TSpectrum(2); 

	  //Create Canvas to show the search result for the SPE
	  if (fFullShow) quad_cuts_ipmt = new TCanvas(Form("quad_cuts_%d",ipmt), Form("First Photoelectron peaks PMT%d",ipmt+1));
	  if (fFullShow) quad_cuts_ipmt->Divide(2,2);
	  //if (fFullShow) quad_cuts_ipmt->cd(1);

	  for (Int_t iquad = 0; iquad < 4; iquad++)
	    {
	      //Perform search for the SPE and save the peak into the array xpeaks
	      if (fFullShow) quad_cuts_ipmt->cd(iquad+1);

	      PulseInt_quad[ipmt][ipmt]->GetXaxis()->SetRangeUser(0,50);
	      fFullShow ? s->Search(PulseInt_quad[iquad][ipmt], 1.0, "nobackground", 0.001) : s->Search(PulseInt_quad[iquad][ipmt], 1.5, "nobackground&&nodraw", 0.001);
	      TList *functions = PulseInt_quad[iquad][ipmt]->GetListOfFunctions();
	      TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
	      Double_t *xpeaks = pm->GetX();
	      PulseInt_quad[iquad][ipmt]->GetXaxis()->SetRangeUser(-1,200);

	      //Use the peak to fit the SPE with a Gaussian to determine the mean
	      Gauss1->SetRange(xpeaks[0]-3, xpeaks[0]+3);
	      Gauss1->SetParameter(1, xpeaks[0]);
	      Gauss1->SetParameter(2, 10.);
	      Gauss1->SetParLimits(0, 0., 2000.);
	      Gauss1->SetParLimits(1, xpeaks[0]-3, xpeaks[0]+3);
	      Gauss1->SetParLimits(2, 0.5, 10.);
	      fFullShow ? PulseInt_quad[iquad][ipmt]->Fit("Gauss1","RQ") : PulseInt_quad[iquad][ipmt]->Fit("Gauss1","RQN");
	      //Store the mean of the SPE in the mean array provided it is not zero, passes a loose statistical cut, and is above a minimum channel number
	      if (xpeaks[0] != 0.0 && PulseInt_quad[iquad][ipmt]->GetBinContent(PulseInt_quad[iquad][ipmt]->GetXaxis()->FindBin(xpeaks[0])) > 10 && ipmt != iquad){	
		mean[iquad] = Gauss1->GetParameter(1);
		mean_err[iquad] = Gauss1->GetParError(1);
	      }
	    } 
	  
	  Double_t xscale = 0.0;
	  Double_t xscaleErr = 0.0;
	  Double_t WeightAvgSum1 = 0.0;
	  Double_t WeightAvgSum2 = 0.0;
	  for (Int_t i=0; i<3; i++)
	    {
	      // Need to define Boolean for this loop, earlier it's based on X2 range of fit, do similar here. Case not used for now SK 27/8/19
	      //if (GoodFit[i] == kFALSE) continue;                     // Take weighted avg 26/8/19 SK
	      WeightAvgSum1 += mean[i]/(mean_err[i]*mean_err[i]);               
	      WeightAvgSum2 += 1/(mean_err[i]*mean_err[i]);
	    }
	  xscale = WeightAvgSum1/WeightAvgSum2;
	  xscaleErr = 1/TMath::Sqrt(WeightAvgSum2);
	  
	  calibration_mk1[ipmt] = xscale;
	  calibration_mk1Err[ipmt] = xscaleErr;

	  //Scale full ADC spectra according to the mean of the SPE. This requires filling a new histogram with the same number of bins but scaled min/max
	  Int_t nbins;
	  nbins = PulseInt_quad[ipmt][ipmt]->GetXaxis()->GetNbins();

	  fscaled[ipmt] = new TH1F(Form("fscaled_PMT%d", ipmt+1), Form("Scaled ADC spectra for PMT%d; NPE; Normalized Counts",ipmt+1), 210, -1, 20);

	  //Fill this histogram bin by bin
	  for (Int_t ibin=0; ibin<nbins; ibin++)
	    {
	      Double_t y = PulseInt_quad[ipmt][ipmt]->GetBinContent(ibin);
	      Double_t x = PulseInt_quad[ipmt][ipmt]->GetXaxis()->GetBinCenter(ibin);
	      Double_t x_scaled = x/calibration_mk1[ipmt];
	      Int_t bin_scaled = fscaled[ipmt]->GetXaxis()->FindBin(x_scaled); 
	      fscaled[ipmt]->SetBinContent(bin_scaled,y);
	    }

	  //Normalize the histogram for ease of fitting
	  fscaled[ipmt]->Scale(1.0/fscaled[ipmt]->Integral(), "width");
	  
	  if (fFullShow) final_spectra_ipmt = new TCanvas(Form("final_Spectra_%d",ipmt), Form("Calibrated spectra for PMT%d; NPE; Normalized Counts",ipmt+1));
	  if (fFullShow) final_spectra_ipmt->cd(1);

	  //Find the location of the SPE and subtract from 1.0 to determine accuracy of calibration
	  Gauss1->SetRange(0.50, 1.50);
	  Gauss1->SetParameter(0, 0.05);
	  Gauss1->SetParameter(1, 1.0);
	  Gauss1->SetParameter(2, 0.3);
	  Gauss1->SetParLimits(0, 0.0, 1.0);
	  Gauss1->SetParLimits(1, 0.5, 1.5);
	  Gauss1->SetParLimits(2, 0.1, 0.5);
	  fFullShow ? fscaled[ipmt]->Fit("Gauss1","RQ") : fscaled[ipmt]->Fit("Gauss1","RQN");
			   
	  //calibration_mk2[ipmt] = calibration_mk1[ipmt]*Gauss1->GetParameter(1);
	  calibration_mk2[ipmt] = xscale * Gauss1->GetParameter(1);

	  calibration_mk2Err[ipmt] = Sqrt(Power(xscaleErr*Gauss1->GetParameter(1),2) +  Power(xscale*Gauss1->GetParError(1),2)); 
	  pmt_calib[ipmt] = abs(1.0 - Gauss1->GetParameter(1));

	  //Scale full ADC spectra according to the mean of the SPE. This requires filling a new histogram with the same number of bins but scaled min/max
	  fscaled_mk2[ipmt] = new TH1F(Form("fscaled_mk2_PMT%d", ipmt+1), Form("Scaled ADC spectra for PMT%d; NPE; Normalized Counts",ipmt+1), 210, -1, 20);
	  
	  //Fill this histogram bin by bin
	  for (Int_t ibin=0; ibin<nbins; ibin++)
	    {
	      Double_t y = PulseInt_quad[ipmt][ipmt]->GetBinContent(ibin);
	      Double_t x = PulseInt_quad[ipmt][ipmt]->GetXaxis()->GetBinCenter(ibin);
	      Double_t x_scaled = x/calibration_mk2[ipmt];
	      Int_t bin_scaled = fscaled_mk2[ipmt]->GetXaxis()->FindBin(x_scaled); 
	      fscaled_mk2[ipmt]->SetBinContent(bin_scaled,y);
	    }

	  //Normalize the histogram for ease of fitting
	  fscaled_mk2[ipmt]->Scale(1.0/fscaled_mk2[ipmt]->Integral(), "width");
	  
	  if (fFullShow) final_spectra_mk2_ipmt = new TCanvas(Form("final_Spectra_mk2_%d",ipmt), Form("Calibrated spectra for PMT%d; NPE; Normalized Counts",ipmt+1));
	  if (fFullShow) final_spectra_mk2_ipmt->cd(1);
	 

	  //Find the location of the SPE and subtract from 1.0 to determine accuracy of calibration
	  Gauss1->SetRange(0.50, 1.50);
	  Gauss1->SetParameter(0, 0.05);
	  Gauss1->SetParameter(1, 1.0);
	  Gauss1->SetParameter(2, 0.3);
	  Gauss1->SetParLimits(0, 0.0, 0.1);
	  Gauss1->SetParLimits(1, 0.5, 1.5);
	  Gauss1->SetParLimits(2, 0.1, 0.5);
	  fFullShow ? fscaled_mk2[ipmt]->Fit("Gauss1","RQ") : fscaled_mk2[ipmt]->Fit("Gauss1","RQN");

	  pmt_calib_mk2[ipmt] = abs(1.0 - Gauss1->GetParameter(1));
	 

      	} //This brace marks the end of TracksFired strategy

      //Begin investigation of Poisson-like behaviour of calibrated spectra..only valid if particle ID is applied
      if (fCut)
	{
	  fscaled_combined[ipmt] = new TH1F(Form("fscaled_combined%d",ipmt+1), Form("Scaled ADC spectra for PMT %d", ipmt+1), 300, 0, 20);

	  fscaled_combined_mk2[ipmt] = new TH1F(Form("fscaled_combined_mk2%d",ipmt+1), Form("Scaled ADC spectra with Second Calibration for PMT %d", ipmt+1), 300, 0, 20);
  
	  Int_t nbins = PulseInt[ipmt]->GetXaxis()->GetNbins();
	  Double_t xmean = calibration_mk1[ipmt];
	  Double_t xmean_mk2 = calibration_mk2[ipmt];

	  fscaled_temp[ipmt] = new TH1F(Form("fscaled_temp_pmt%d",ipmt+1), Form("Scaled ADC spectra for PMT %d", ipmt+1), 300, 0, 20);
	  fscaled_temp_mk2[ipmt] = new TH1F(Form("fscaled_temp_mk2_pmt%d",ipmt+1), Form("Scaled ADC spectra for PMT %d", ipmt+1), 300, 0, 20);

	  //Fill this histogram bin by bin
	  for (Int_t ibin=0; ibin < nbins; ibin++)
	    {
	      Double_t y = PulseInt[ipmt]->GetBinContent(ibin);
	      Double_t x = PulseInt[ipmt]->GetXaxis()->GetBinCenter(ibin);
	      Double_t x_scaled_mk1 = x/xmean;
	      Double_t x_scaled_mk2 = x/xmean_mk2;
	      Int_t bin_scaled_mk1 = fscaled_temp[ipmt]->GetXaxis()->FindBin(x_scaled_mk1); 
	      Int_t bin_scaled_mk2 = fscaled_temp_mk2[ipmt]->GetXaxis()->FindBin(x_scaled_mk2);
	      fscaled_temp[ipmt]->SetBinContent(bin_scaled_mk1,y);
	      fscaled_temp_mk2[ipmt]->SetBinContent(bin_scaled_mk2,y);
	    }
	  fscaled_combined[ipmt]->Add(fscaled_temp[ipmt]);
	  fscaled_combined_mk2[ipmt]->Add(fscaled_temp_mk2[ipmt]);	

	  //Normalize the histogram for ease of fitting
	  fscaled_combined[ipmt]->Scale(1.0/fscaled_combined[ipmt]->Integral(), "width");
	  fscaled_combined_mk2[ipmt]->Scale(1.0/fscaled_combined[ipmt]->Integral(), "width");
	}     // This brace marks the end of the loop over PMTs

              // Write to PDF file

      /*   if (ipmt == 0){
	quad_cuts[ipmt]->Print(outputpdf + "[");
	quad_cuts[ipmt]->Print(outputpdf);
      }
       if (ipmt == 1 ){
	quad_cuts[ipmt]->Print(outputpdf);
      }
       if (ipmt == 2 ){
	quad_cuts[ipmt]->Print(outputpdf);
      }
       if (ipmt == 3){
	quad_cuts[ipmt]->Print(outputpdf);
	quad_cuts[0]->Print(outputpdf + "]");
	} */  
     
    }   // Combine each PMT into one final histogram

  if (fCut)
    {
      fscaled_total = new TH1F("fscaled_total", "Scaled ADC spectra for all PMTs;NPE;Normalized Counts", 300, 0, 20);
      fscaled_total_mk2 = new TH1F("fscaled_total_mk2", "Scaled ADC spectra for all PMTs;NPE;Normalized Counts", 300, 0, 20);
      for (Int_t i=0; i<4; i++)
	{
	  fscaled_total->Add(fscaled_combined[i]);
	  fscaled_total_mk2->Add(fscaled_combined_mk2[i]);
	}

      fscaled_total->Scale(1.0/fscaled_total->Integral(), "width");
      fscaled_total_mk2->Scale(1.0/fscaled_total_mk2->Integral(), "width");

      //Display the Poisson characteristics of the ADC spectra
      if (fFullShow) scaled_total = new TCanvas("scaled_total", "Scaled ADC of all PMTs showing Poisson Fit");
      if (fFullShow) scaled_total->Divide(2,1);
      if (fFullShow) scaled_total->cd(1);
      Poisson->SetRange(0, 20);
      Poisson->SetParameter(0, Poisson_mean);
      Poisson->SetParameter(1, 0.8);
      Poisson->SetParLimits(0, Poisson_mean - 1.0, Poisson_mean + 3.0);
      fFullShow ? fscaled_total->Fit("Poisson","RQ") : fscaled_total->Fit("Poisson","RQN");
      Pois_Chi[0] = Poisson->GetChisquare();
      if (fFullShow) scaled_total->cd(2);
      fFullShow ? fscaled_total_mk2->Fit("Poisson","RQ") : fscaled_total_mk2->Fit("Poisson","RQN");
      Pois_Chi[1] = Poisson->GetChisquare();
      } 

  printf("\n\n"); 

  //Output the actual calibration information
  //cout << "Calibration constants are (where the '*' indicates the better value)\nPMT#: First Guess  Second Guess\n" << endl;
  cout << "Calibration constants are \nPMT#:   First Guess" << setw(25) << "Second Guess\n" << endl;
  for (Int_t i=0; i<4; i++)
    {
      cout << Form("PMT%d:", i+1) << setw(8) << Form("%3.3f", calibration_mk1[i]) << " +/- " << Form("%3.3f", calibration_mk1Err[i]) << setw(13) << Form("%3.3f", calibration_mk2[i]) << " +/- " << Form("%3.3f", calibration_mk2Err[i]) <<  "\n";
    }

  printf("\n");

  if (fCut) cout << (Pois_Chi[0] < Pois_Chi[1] ? "First Guess":"Second Guess") << " better characterizes the full Poisson character" << endl;

  //Start the process of writing the calibration information to file
 
  ofstream calibration;
  calibration.open("calibration_temp.txt", ios::out);

  if (!calibration.is_open()) cout << "Problem saving calibration constants, may have to update constants manually!" << endl;

  else
    {
      calibration << "phgcer_adc_to_npe = "; 
      for (Int_t ipmt = 0; ipmt < 4; ipmt++)
	{
	  calibration << Form("1./%3.3f. ", (pmt_calib[ipmt] < pmt_calib_mk2[ipmt]) ? calibration_mk1[ipmt] : calibration_mk2[ipmt]);
	}

      calibration.close();
    }
  
}
