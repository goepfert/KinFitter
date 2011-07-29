#define TToyGentT_cxx

#include <iostream>

#include "KinFitter/TToyGentT.h"
#include "KinFitter/TFitConstraintM.h"
#include "KinFitter/TFitConstraintMGaus.h"
#include "KinFitter/TFitConstraintMBW.h"
#include "KinFitter/TAbsFitParticle.h"
#include "KinFitter/TKinFitter.h"

#include "TBenchmark.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TString.h"

ClassImp(TToyGentT)

TToyGentT::TToyGentT( const TAbsFitParticle* Whad_d1, 
		      const TAbsFitParticle* Whad_d2, 
		      const TAbsFitParticle* Tophad_d2, 
		      const TAbsFitParticle* Wlep_d1, 
		      const TAbsFitParticle* Wlep_d2, 
		      const TAbsFitParticle* Toplep_d2 ) 
  : TNamed("NoName", "NoTitle") {
  
  // Clone input particles
  _iniWhad_d1   = Whad_d1->clone( Whad_d1->GetName() + (TString) "INI" );
  _Whad_d1      = Whad_d1->clone( Whad_d1->GetName() + (TString) "SMEAR" );
  _iniWhad_d2   = Whad_d2->clone( Whad_d2->GetName() + (TString) "INI" );
  _Whad_d2      = Whad_d2->clone( Whad_d2->GetName() + (TString) "SMEAR" );
  _iniTophad_d2 = Tophad_d2->clone( Tophad_d2->GetName() + (TString) "INI" );
  _Tophad_d2    = Tophad_d2->clone( Tophad_d2->GetName() + (TString) "SMEAR" );

  _iniWlep_d1   = Wlep_d1->clone( Wlep_d1->GetName() + (TString) "INI" );
  _Wlep_d1      = Wlep_d1->clone( Wlep_d1->GetName() + (TString) "SMEAR" );
  _iniWlep_d2   = Wlep_d2->clone( Wlep_d2->GetName() + (TString) "INI" );
  _Wlep_d2      = Wlep_d2->clone( Wlep_d2->GetName() + (TString) "SMEAR" );
  _iniToplep_d2 = Toplep_d2->clone( Toplep_d2->GetName() + (TString) "INI" );
  _Toplep_d2    = Toplep_d2->clone( Toplep_d2->GetName() + (TString) "SMEAR" );
  
  _inimeasParticles.clear();
  _iniunmeasParticles.clear();
  _measParticles.clear();
  _unmeasParticles.clear();

  _measConstraints.clear();

  _printPartIni = false;
  _printConsIni = false;
  _printSmearedPartBefore = false;
  _printConsBefore = false;
  _printConsAfter = false;
  _printPartAfter = false;
  _doCheckConstraintsTruth = true;

  _verbosity=1;
}

TToyGentT::~TToyGentT() {

  delete _iniWhad_d1;
  delete _iniWhad_d2;
  delete _iniTophad_d2;
  delete _iniWlep_d1;
  delete _iniWlep_d2;
  delete _iniToplep_d2;

  delete _Whad_d1;
  delete _Whad_d2;
  delete _Tophad_d2;
  delete _Wlep_d1;
  delete _Wlep_d2;
  delete _Toplep_d2;
}

Bool_t TToyGentT::doToyExperiments( Int_t nbExperiments ) {

  // define fitter
  TKinFitter fitter;

  fitter.addMeasParticle(_Whad_d1);
  _inimeasParticles.push_back( _iniWhad_d1 );
  _measParticles.push_back( _Whad_d1 );
  fitter.addMeasParticle(_Whad_d2);
  _inimeasParticles.push_back( _iniWhad_d2 );
  _measParticles.push_back( _Whad_d2 );

  fitter.addMeasParticle(_Tophad_d2);
  _inimeasParticles.push_back( _iniTophad_d2 );
  _measParticles.push_back( _Tophad_d2 );

  fitter.addMeasParticle(_Wlep_d1);
  _inimeasParticles.push_back( _iniWlep_d1 );
  _measParticles.push_back( _Wlep_d1 );

  fitter.addMeasParticle(_Wlep_d2);
  /**
     set the right unmeasured parameter of the neutrino
     this depends on your parametrization
     e.g. par 1 for PtEtaPhi, par 2 for PxPyPz
   */
  fitter.setParamUnmeas( _Wlep_d2, 1);

  _inimeasParticles.push_back( _iniWlep_d2 );
  _measParticles.push_back( _Wlep_d2 );
  fitter.addMeasParticle(_Toplep_d2);
  _inimeasParticles.push_back( _iniToplep_d2 );
  _measParticles.push_back( _Toplep_d2 );

  double WhadM = ((*_Whad_d1->getIni4Vec()) + (*_Whad_d2->getIni4Vec())).M();
  double TophadM = ((*_Whad_d1->getIni4Vec()) + (*_Whad_d2->getIni4Vec()) + (*_Tophad_d2->getIni4Vec())).M();
  double WlepM = ((*_Wlep_d1->getIni4Vec()) + (*_Wlep_d2->getIni4Vec())).M();
  double ToplepM = ((*_Wlep_d1->getIni4Vec()) + (*_Wlep_d2->getIni4Vec()) + (*_Toplep_d2->getIni4Vec())).M();
  
  TFitConstraintM MConsWhad( "MassConstraint Whad", "Mass-Constraint Whad", 0, 0, WhadM);
  MConsWhad.addParticles1( _Whad_d1, _Whad_d2 );
 
  TFitConstraintM MConsTophad( "MassConstraint Tophad", "Mass-Constraint Tophad", 0, 0, TophadM);
  MConsTophad.addParticles1( _Whad_d1, _Whad_d2, _Tophad_d2 );

  TFitConstraintM MConsWlep( "MassConstraint Wlep", "Mass-Constraint Wlep", 0, 0, WlepM);
  MConsWlep.addParticles1( _Wlep_d1, _Wlep_d2 );

  TFitConstraintM MConsToplep( "MassConstraint Toplep", "Mass-Constraint Toplep", 0, 0, ToplepM);
  MConsToplep.addParticles1( _Wlep_d1, _Wlep_d2, _Toplep_d2 );

  _measConstraints.push_back( &MConsWhad );
  _measConstraints.push_back( &MConsWlep );
  _measConstraints.push_back( &MConsTophad );
  _measConstraints.push_back( &MConsToplep );

  fitter.addConstraint(&MConsWhad);
  fitter.addConstraint(&MConsWlep);
  fitter.addConstraint(&MConsTophad);
  fitter.addConstraint(&MConsToplep);

  /** equal mass constraint of leptonic and hadronic tops */
  //  TFitConstraintM MConsTops( "MassConstraint Toplep=Tophad", "Mass-Constraint Toplep=Tophad", 0, 0, 0);
  //  MConsTops.addParticles1( _Whad_d1, _Whad_d2, _Tophad_d2 );
  //  MConsTops.addParticles2( _Wlep_d1, _Wlep_d2, _Toplep_d2 );
  //  _measConstraints.push_back( &MConsTops );
  //  fitter.addConstraint(&MConsTops);

  fitter.setMaxNbIter( 50 );
  fitter.setMaxDeltaS( 5e-5 );
  fitter.setMaxF( 1e-4 );
  fitter.setVerbosity( _verbosity );

  if( _printPartIni ) {
    cout << endl
   	 << "----------------------------------" << endl;
    cout << "--- PRINTING INITIAL PARTICLES ---" << endl;
    cout << "----------------------------------" << endl ;
    _iniWhad_d1->print();
    _iniWhad_d2->print();
    _iniTophad_d2->print();
    _iniWlep_d1->print();
    _iniWlep_d2->print();
    _iniToplep_d2->print();
    cout << endl << endl;
  }
 
  if( _printConsIni ) {
    cout << endl
	 << "-------------------------------------------------" << endl;
    cout << "INITIAL CONSTRAINTS " << endl ;
    cout << "-------------------------------------------------" << endl;
    cout << "     M Whad: "   << MConsWhad.getCurrentValue()
	 << "     M Tophad: " << MConsTophad.getCurrentValue()
	 << "     M Wlep: "   << MConsWlep.getCurrentValue()
	 << "     M Tophad: " << MConsToplep.getCurrentValue() << endl << endl;
  }
  
  // Check initial constraints
  if (  _doCheckConstraintsTruth ) {
    if (fitter.getF() > fitter.getMaxF()) {
      cout << "Initial constraints are not fulfilled." << endl;
      //return false;
    }
  }
  
  // create histograms
  createHists(&fitter);

  // Start Benchmark
  gBenchmark->Start("ToyExperiments");

  // perform pseudo experiments
  for (int i = 0; i < nbExperiments; i++) {

    smearParticles();
    
    if( _printSmearedPartBefore ) {
      cout <<  endl  
	   << "-------------------------------------------------------" << endl ;
      cout << "--- PRINTING SMEARED PARTICLES BEFORE FIT FOR experiment # " <<i+1 << endl;
      cout << "-------------------------------------------------------" << endl;
      _Whad_d1->print();
      _Whad_d2->print();
      _Tophad_d2->print();
      _Wlep_d1->print();
      _Wlep_d2->print();
      _Toplep_d2->print();
    }
    
    if( _printConsBefore ) {
      cout << endl
	   << "-------------------------------------------------" << endl;
      cout << "INITIAL (SMEARED) CONSTRAINTS FOR experiment # "<< i+1 << endl ;
      cout << "-------------------------------------------------" << endl;
      cout << "     M Whad: "   << MConsWhad.getCurrentValue()
	   << "     M Tophad: " << MConsTophad.getCurrentValue()
	   << "     M Wlep: "   << MConsWlep.getCurrentValue()
	   << "     M Tophad: " << MConsToplep.getCurrentValue() << endl << endl;
    }
    
    fitter.fit();
    
    if( _printConsAfter) {
      cout << endl
	   << "-------------------------------------------------" << endl;
      cout << " CONSTRAINTS AFTER FIT FOR experiment # "<< i+1 << endl ;
      cout << "-------------------------------------------------" << endl;
      cout << "     M Whad: "   << MConsWhad.getCurrentValue()
	   << "     M Tophad: " << MConsTophad.getCurrentValue()
	   << "     M Wlep: "   << MConsWlep.getCurrentValue()
	   << "     M Tophad: " << MConsToplep.getCurrentValue() << endl << endl;
    }
    
    if( _printPartAfter ) {
      cout <<  endl  
	   << "--------------------------------------------------------" << endl ;
      cout << "--- PRINTING PARTICLES AFTER FIT FOR experiment # "<< i+1 << endl;
      cout << "--------------------------------------------------------" << endl;
      _Whad_d1->print();
      _Whad_d2->print();
      _Tophad_d2->print();
      _Wlep_d1->print();
      _Wlep_d2->print();
      _Toplep_d2->print();
    }
    
    _histStatus->Fill( fitter.getStatus() );
    _histNIter->Fill( fitter.getNbIter() );
    if ( fitter.getStatus() == 0 ) {
      _histPChi2->Fill( TMath::Prob( fitter.getS(), fitter.getNDF() ) );
      _histChi2->Fill( fitter.getS());
      fillPull1();
      fillPull2();
      fillPar();
      fillM();
      fillConstraints();
    }

    if (i % 176 == 0) {
      cout << "\r";  
      cout <<" ------ "<< (Double_t) i/nbExperiments*100. << " % PROCESSED ------";
      cout.flush();
    }    
  }

  // Stop Benchmark
  cout << endl;
  gBenchmark->Stop("ToyExperiments");
  gBenchmark->Show("ToyExperiments");
  
  return true;
}


void TToyGentT::smearParticles() {
  // Smear measured particles

  for (UInt_t p = 0; p < _measParticles.size(); p++) {
    TAbsFitParticle* particle =  _measParticles[p];
    TAbsFitParticle* iniParticle =  _inimeasParticles[p];
    TMatrixD parIni( *(iniParticle->getParIni()) );
    const TMatrixD* covM = iniParticle->getCovMatrix();
    for (int m = 0; m < iniParticle->getNPar(); m++) {
      parIni(m, 0) += gRandom->Gaus(0., TMath::Sqrt( (*covM)(m,m) ) );
    }
    TLorentzVector* ini4Vec = iniParticle->calc4Vec( &parIni );
    particle->setIni4Vec( ini4Vec );
    delete ini4Vec;
  }
}


void TToyGentT::fillConstraints(){

  int histIdx=0;
  for (UInt_t c = 0; c <  _measConstraints.size(); c++) {
    //TAbsFitConstraint* constraints = _measConstraints.at(c);
    
    const TMatrixD* CovPar = _measConstraints.at(c)->getCovMatrixDeltaAlpha();
    
    if(CovPar){
      TMatrixD pardiff( *(_measConstraints.at(c)->getParCurr()) );
      pardiff -= *(_measConstraints.at(c)->getParIni());

      for (int i = 0; i < pardiff.GetNrows(); i++) {
	( (TH1D*) _histsDiff3[histIdx])->Fill( pardiff(i, 0) );
	pardiff(i, 0) /= TMath::Sqrt( (*CovPar)(i, i) );
	( (TH1D*) _histsPull3[histIdx])->Fill( pardiff(i, 0) );
	( (TH1D*) _histsError3[histIdx])->Fill( TMath::Sqrt( (*CovPar)(i, i) ) );
	
	histIdx++;
      }
    }
    
  }//end constraints
  
}


void TToyGentT::fillPull1() {

  Int_t histindex = 0;
  for (UInt_t p = 0; p < _measParticles.size(); p++) {

    TLorentzVector vectrue( *_inimeasParticles[p]->getIni4Vec() );
    TMatrixD* partrue = _measParticles[p]->transform( vectrue );

    const TMatrixD* parfit = _measParticles[p]->getParCurr();
    TMatrixD parpull( *parfit );
    parpull -= (*partrue);
    const TMatrixD* covMatrixFit = _measParticles[p]->getCovMatrixFit();
    
    for (int i = 0; i < parpull.GetNrows(); i++) {
      if( !TMath::IsNaN(parpull(i, 0))         )((TH1D*) _histsDiff1[histindex])->Fill( parpull(i, 0) );
      parpull(i, 0) /= TMath::Sqrt( (*covMatrixFit)(i, i) );
      if( !TMath::IsNaN(parpull(i, 0))         )((TH1D*) _histsPull1[histindex])->Fill( parpull(i, 0) );
      if( !TMath::IsNaN((*covMatrixFit)(i, i)) )((TH1D*) _histsError1[histindex])->Fill( TMath::Sqrt( (*covMatrixFit)(i, i) ) );
      histindex++;
    }
  }

}


void TToyGentT::fillPull2() {

  Int_t histindex = 0;
  for (UInt_t p = 0; p <  _measParticles.size(); p++) {

    const TMatrixD* pull = _measParticles[p]->getPull();
    const TMatrixD* VDeltaY = _measParticles[p]->getCovMatrixDeltaY();
    TMatrixD pardiff( *(_measParticles[p]->getParCurr()) );
    pardiff -= *(_measParticles[p]->getParIni());
    for (int i = 0; i < pull->GetNrows(); i++) {
      if( !TMath::IsNaN((*pull)(i, 0))     )( (TH1D*) _histsPull2[histindex])->Fill( (*pull)(i, 0) );
      if( !TMath::IsNaN((*VDeltaY)(i, i) ) )( (TH1D*) _histsError2[histindex])->Fill( TMath::Sqrt( (*VDeltaY)(i, i) ) );
      if( !TMath::IsNaN(pardiff(i, 0))     )( (TH1D*) _histsDiff2[histindex])->Fill( pardiff(i, 0) );
      histindex++;
    }
  }
  
}


void TToyGentT::fillPar() {

  Int_t histindex = 0;
  for (UInt_t p = 0; p <  _measParticles.size(); p++) {

    const TMatrixD* partrue  = _inimeasParticles[p]->getParIni();
    const TMatrixD* parsmear = _measParticles[p]->getParIni();
    const TMatrixD* parfit   = _measParticles[p]->getParCurr();
    for (int i = 0; i < partrue->GetNrows(); i++) {
      if( !TMath::IsNaN((*partrue)(i, 0))  )( (TH1D*) _histsParTrue[histindex])->Fill( (*partrue)(i, 0) );
      if( !TMath::IsNaN((*parsmear)(i, 0)) )( (TH1D*) _histsParSmear[histindex])->Fill( (*parsmear)(i, 0) );
      if( !TMath::IsNaN((*parfit)(i, 0))   )( (TH1D*) _histsParFit[histindex])->Fill( (*parfit)(i, 0) );
      histindex++;
    }
  }

}


void TToyGentT::fillM() {

  _histMWhadTrue->Fill( ( *(_iniWhad_d1->getIni4Vec()) + *(_iniWhad_d2->getIni4Vec()) ).M() );
  _histMWhadSmear->Fill( ( *(_Whad_d1->getIni4Vec()) + *(_Whad_d2->getIni4Vec()) ).M() );
  _histMWhadFit->Fill( ( *(_Whad_d1->getCurr4Vec()) + *(_Whad_d2->getCurr4Vec()) ).M() );
  
  _histMTophadTrue->Fill( ( *(_iniWhad_d1->getIni4Vec()) + *(_iniWhad_d2->getIni4Vec()) + *(_iniTophad_d2->getIni4Vec()) ).M() );
  _histMTophadSmear->Fill( ( *(_Whad_d1->getIni4Vec()) + *(_Whad_d2->getIni4Vec()) + *(_Tophad_d2->getIni4Vec()) ).M() );
  _histMTophadFit->Fill( ( *(_Whad_d1->getCurr4Vec()) +  *(_Whad_d2->getCurr4Vec()) + *(_Tophad_d2->getCurr4Vec()) ).M() );
  
  _histMWlepTrue->Fill( ( *(_iniWlep_d1->getIni4Vec()) + *(_iniWlep_d2->getIni4Vec()) ).M() );
  _histMWlepSmear->Fill( ( *(_Wlep_d1->getIni4Vec()) + *(_Wlep_d2->getIni4Vec()) ).M() );
  _histMWlepFit->Fill( ( *(_Wlep_d1->getCurr4Vec()) + *(_Wlep_d2->getCurr4Vec()) ).M() );
  
  _histMToplepTrue->Fill( ( *(_iniWlep_d1->getIni4Vec()) + *(_iniWlep_d2->getIni4Vec()) + *(_iniToplep_d2->getIni4Vec()) ).M() );
  _histMToplepSmear->Fill( ( *(_Wlep_d1->getIni4Vec()) + *(_Wlep_d2->getIni4Vec()) + *(_Toplep_d2->getIni4Vec()) ).M() );
  _histMToplepFit->Fill( ( *(_Wlep_d1->getCurr4Vec()) +  *(_Wlep_d2->getCurr4Vec()) + *(_Toplep_d2->getCurr4Vec()) ).M() );

}


void  TToyGentT::createHists(TKinFitter* fitter) {

  // fitter
  _histStatus = new TH1D( "hStatus", "Status of the Fit; fit status", 16, -1, 15);
  _histNIter = new TH1D( "hNIter", "Number of iterations; # iterations", 100, 0, 100);
  _histPChi2 = new TH1D( "hPChi2", "Chi2 probability; #chi^{2} probability", 100, 0., 1.);
  _histChi2 = new TH1D(  "hChi2", "Chi2; #chi^{2} ", 200, 0., 20.);

  // masses of composites before and after fit
  _histMWhadTrue  = new TH1D( "histMWhadTrue",  "histMWhadTrue; m(jj) [GeV]",  200, 0., 180.);
  _histMWhadSmear = new TH1D( "histMWhadSmear", "histMWhadSmear; m(jj) [GeV]", 200, 40., 140.);
  _histMWhadFit   = new TH1D( "histMWhadFit",   "histMWhadFit; m(jj) [GeV]",   200, 40., 140.);
  _histMTophadTrue  = new TH1D( "histMTophadTrue",  "histMTophadTrue; m(jjj) [GeV]",  200, 80., 280.);
  _histMTophadSmear = new TH1D( "histMTophadSmear", "histMTophadSmear; m(jjj) [GeV]", 200, 80., 280.);
  _histMTophadFit   = new TH1D( "histMTophadFit",   "histMTophadFit; m(jjj) [GeV]",   200, 80., 280.);

  _histMWlepTrue  = new TH1D( "histMWlepTrue",  "histMWlepTrue; m(lnu) [GeV]",  200, 0., 180.);
  _histMWlepSmear = new TH1D( "histMWlepSmear", "histMWlepSmear; m(lnu) [GeV]", 200, 0., 180.);
  _histMWlepFit   = new TH1D( "histMWlepFit",   "histMWlepFit; m(lnu) [GeV]",   200, 0., 180.);
  _histMToplepTrue  = new TH1D( "histMToplepTrue",  "histMToplepTrue; m(lnuj) [GeV]",  200, 80., 280.);
  _histMToplepSmear = new TH1D( "histMToplepSmear", "histMToplepSmear; m(lnuj) [GeV]", 200, 80., 280.);
  _histMToplepFit   = new TH1D( "histMToplepFit",   "histMToplepFit; m(lnuj) [GeV]",   200, 80., 280.);


  _histsParTrue.Clear();
  _histsParSmear.Clear();
  _histsParFit.Clear();

  _histsPull1.Clear();
  _histsError1.Clear();
  _histsDiff1.Clear();
  _histsPull2.Clear();
  _histsError2.Clear();
  _histsDiff2.Clear();

  // TObjArray of TObjArrays
  TObjArray histarrays;
  histarrays.Add( &_histsParTrue );
  histarrays.Add( &_histsParSmear );
  histarrays.Add( &_histsParFit );
  histarrays.Add( &_histsPull1 );
  histarrays.Add( &_histsError1 );
  histarrays.Add( &_histsDiff1 );
  histarrays.Add( &_histsPull2 );
  histarrays.Add( &_histsError2 );
  histarrays.Add( &_histsDiff2 );

  TString histnames[]  = {"hParTrue", "hParSmear", "hParFit", "hPull1", "hError1", "hDiff1", "hPull2", "hError2", "hDiff2" };
  TString histtitles[] = {"True", "Smear", "Fit", "Pull1", "Error1", "Diff1", "Pull2", "Error2", "Diff2" };

  // Histogram bounds
  TArrayD arrmin( histarrays.GetEntries() );
  TArrayD arrmax( histarrays.GetEntries() );
  arrmin[0] = 0.;   arrmax[0] = 2.;   // Parameters, true
  arrmin[1] = 0.;   arrmax[1] = 2.;   // Parameters, smear
  arrmin[2] = 0.;   arrmax[2] = 2.;   // Parameters, fit

  arrmin[3] = -3.;  arrmax[3] = 3.;   // Pull1
  arrmin[4] = 0.;   arrmax[4] = 10.;  // Error1
  arrmin[5] = -2.5; arrmax[5] = 2.5;  // Diff1
  arrmin[6] = -3.;  arrmax[6] = 3.;   // Pull2
  arrmin[7] = 10.;  arrmax[7] = 10.;  // Error2
  arrmin[8] = -2.5; arrmax[8] = 2.5;  // Diff2

  for (UInt_t p = 0; p <  _measParticles.size(); p++) {

    TAbsFitParticle* particle = _measParticles[p];
    const TMatrixD* covMatrix = particle->getCovMatrix();

    for (int i = 0; i < particle->getNPar(); i++) {
      for (int h = 0; h < histarrays.GetEntries(); h++ ) {
	TString name = histnames[h] + "_" + (TString) particle->GetName()  + Form("_%i",i);
	TString title = (TString) particle->GetName() + "_" + histtitles[h] + Form("_Par_%i",i);
	title.ReplaceAll("SMEAR","");
	if ( h < 3) {
	  const TMatrixD* parfit = _measParticles[p]->getParCurr();
	  if( (fitter->getParamMeas(particle,i)) ){
	    arrmin[h] = (*parfit)(i,0)-4*TMath::Sqrt((*covMatrix)(i,i));
	    arrmax[h] = (*parfit)(i,0)+4*TMath::Sqrt((*covMatrix)(i,i));	  
	  } else {
	    arrmin[h] = -0.1 ;
	    arrmax[h] = 0.1;
	  }
	}
	TH1D* newhisto =  new TH1D( name, title, 100, arrmin[h], arrmax[h]) ; 
	((TObjArray*) histarrays[h])->Add( newhisto );
      }
    }
  }//end measPaticles


  _histsPull3.Clear();
  _histsError3.Clear();
  _histsDiff3.Clear();

  // TObjArray of TObjArrays
  TObjArray histsConstraints;
  histsConstraints.Add( &_histsPull3 );
  histsConstraints.Add( &_histsError3 );
  histsConstraints.Add( &_histsDiff3 );
  TString histnames3[]  = {"Pull3", "Error3", "Diff3"};

  for (UInt_t c = 0; c < _measConstraints.size(); c++) {

    TAbsFitConstraint* constraints = _measConstraints.at(c);
    for (int i = 0; i < constraints->getNPar(); i++) {
      for (int h = 0; h < histsConstraints.GetEntries(); h++ ) {
	TString name = histnames3[h] + "_" + (TString) constraints->GetName()+ Form("_%i",i);
	
	TH1D* newhisto =  new TH1D( name, name, 100, -3, 3) ; 
	((TObjArray*) histsConstraints[h])->Add( newhisto );
      }
    }
  }//end constraints
  
}
