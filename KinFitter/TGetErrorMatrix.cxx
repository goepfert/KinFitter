
using namespace std;

#include "KinFitter/TGetErrorMatrix.h"
#include "KinFitter/TAbsFitParticle.h"


#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TStyle.h"
#include "TPaveText.h"


ClassImp(TGetErrorMatrix)
  
  
TGetErrorMatrix::TGetErrorMatrix( const TAbsFitParticle* theParticle, 
				  TString VarX, Int_t NbinsX, Double_t xlow, Double_t xup,
				  TString VarY, Int_t NbinsY, Double_t ylow, Double_t yup,
				  Int_t histResolNBins, Double_t histResolMin, Double_t histResolMax ):
  TNamed( "NoName", "NoTitle")
//   _VarX(VarX),
//   _VarY(VarY),
//   _xlow(xlow),
//   _xup(xup),
//   _NbinsX(NbinsX),
//   _ylow(ylow),
//   _yup(yup),
//   _NbinsY(NbinsY),
//   _BinWidthX((_xup - _xlow) / _NbinsX),
//   _BinWidthY((_yup - _ylow) / _NbinsY),
//   _xAxisTitle(VarX),
//   _yAxisTitle(VarY)
{ 

  // Constructor takes the following parameters:
  //
  // theParticle:         particle with four-vector to calculate uncertainties for
  // VarX:                name of the X-axis 
  // NbinsX, xlow, xup:   number of bins, lower and upper X-axis boundaries
  // VarY:                name of the Y-axis 
  // NbinsY, ylow, yup:   number of bins, lower and upper X-axis boundaries
  // histResolNBins:      number of bins if resolution histograms
  // histResolMin, histResolMax: minimum and maximum of resolution histograms

  _verbose = false;

  _theParticle = theParticle->clone( theParticle->GetName() );
  _nPar = _theParticle->getNPar();
  _xAxis = new TAxis( NbinsX, xlow, xup );
  _xAxis->SetNameTitle( VarX, VarX );
  _yAxis = new TAxis( NbinsY, ylow, yup );
  _yAxis->SetNameTitle( VarY, VarY );
  _resolAxis = new TAxis( histResolNBins, histResolMin, histResolMax );
  _resolAxis->SetNameTitle( "Resol", "(meas. - true)" );

  CreateHistos();  
  DefineFitFuncs();
} 

TGetErrorMatrix::TGetErrorMatrix(const TString &name, const TString &title,
				 const TAbsFitParticle* theParticle,
				 TString VarX, Int_t NbinsX, Double_t xlow, Double_t xup,
				 TString VarY, Int_t NbinsY, Double_t ylow, Double_t yup,
				 Int_t histResolNBins, Double_t histResolMin, Double_t histResolMax ):
  TNamed( name, title)
{ 
  // Constructor takes the following parameters:
  //
  // name, title:         name and title of object to create
  // theParticle:         particle with four-vector to calculate uncertainties for
  // VarX:                name of the X-axis 
  // NbinsX, xlow, xup:   number of bins, lower and upper X-axis boundaries
  // VarY:                name of the Y-axis 
  // NbinsY, ylow, yup:   number of bins, lower and upper X-axis boundaries
  // histResolNBins:      number of bins if resolution histograms
  // histResolMin, histResolMax: minimum and maximum of resolution histograms

  _verbose = false;
  _theParticle = theParticle->clone( theParticle->GetName() );
  _nPar = _theParticle->getNPar();
  _xAxis = new TAxis( NbinsX, xlow, xup );
  _xAxis->SetNameTitle( VarX, VarX );
  _yAxis = new TAxis( NbinsY, ylow, yup );
  _yAxis->SetNameTitle( VarY, VarY );
  _resolAxis = new TAxis( histResolNBins, histResolMin, histResolMax );
  _resolAxis->SetNameTitle( "Resol", "(meas. - true)" );

  CreateHistos();  
  DefineFitFuncs();  
}

TGetErrorMatrix::TGetErrorMatrix( const TString &name, const TString &title,
				  const TAbsFitParticle* theParticle,
				  TAxis* xAxis, TAxis* yAxis, TAxis* resolAxis ):
  TNamed( name, title)
{ 
  // Constructor takes the following parameters:
  //
  // name, title:         name and title of object to create
  // theParticle:         particle with four-vector to calculate uncertainties for
  // VarX:                name of the X-axis 
  // NbinsX, xlow, xup:   number of bins, lower and upper X-axis boundaries
  // VarY:                name of the Y-axis 
  // NbinsY, ylow, yup:   number of bins, lower and upper X-axis boundaries
  // histResolNBins:      number of bins if resolution histograms
  // histResolMin, histResolMax: minimum and maximum of resolution histograms

  _verbose = false;
  _theParticle = theParticle->clone( theParticle->GetName() );
  _nPar = _theParticle->getNPar();
  _xAxis = (TAxis*) xAxis->Clone();
  _yAxis = (TAxis*) yAxis->Clone();
  _resolAxis = (TAxis*) resolAxis->Clone();

  CreateHistos();  
  DefineFitFuncs();  
}

TGetErrorMatrix::TGetErrorMatrix():
  TNamed( "NoName", "NoTitle" )
{
  _nPar = 0;
  _verbose = false;
  _xAxis = new TAxis( 1, 0., 1. );
  _xAxis->SetNameTitle( "getErrorMatrixX", "getErrorMatrixX" );
  _yAxis = new TAxis( 1, 0., 1. );
  _yAxis->SetNameTitle( "getErrorMatrixY", "getErrorMatrixY" );
  _resolAxis = new TAxis( 120, -3., 3. );
  _resolAxis->SetNameTitle( "Resol", "(meas. - true)" );

  DefineFitFuncs();

}

TGetErrorMatrix::~TGetErrorMatrix() {

}

void TGetErrorMatrix::DefineFitFuncs() {

  TF2* pol3 = new TF2("mypol3"," [0]*x + [1]*x*x + [2]*x*x*x + [3]*y + [4]*y*y + [5]*y*y*y ", 
		      _xAxis->GetBinLowEdge(1), _xAxis->GetBinUpEdge(_xAxis->GetNbins()), 
		      _yAxis->GetBinLowEdge(1), _yAxis->GetBinUpEdge(_xAxis->GetNbins()));
  TF2* pol4 = new TF2("mypol4"," [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*y + [5]*y*y + [6]*y*y*y + [7]*y*y*y*y ",
		      _xAxis->GetBinLowEdge(1), _xAxis->GetBinUpEdge(_xAxis->GetNbins()), 
		      _yAxis->GetBinLowEdge(1), _yAxis->GetBinUpEdge(_xAxis->GetNbins()));
  TF2* pol5 = new TF2("mypol5",
		      "[0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*x*x*x*x*x +  [5]*y + [6]*y*y + [7]*y*y*y + [8]*y*y*y*y + [9]*y*y*y*y*y ", 
		      _xAxis->GetBinLowEdge(1), _xAxis->GetBinUpEdge(_xAxis->GetNbins()), 
		      _yAxis->GetBinLowEdge(1), _yAxis->GetBinUpEdge(_xAxis->GetNbins()));
  TF2* pol9 = new TF2("mypol9", "xpol8 + ypol8(9)", 
		      _xAxis->GetBinLowEdge(1), _xAxis->GetBinUpEdge(_xAxis->GetNbins()), 
		      _yAxis->GetBinLowEdge(1), _yAxis->GetBinUpEdge(_xAxis->GetNbins()));
  
}

void TGetErrorMatrix::CreateHistos() {
  // Creates histograms

  _histarrays.Add( &_histsPar1 );
  _histarrays.Add( &_histsPar2 );
  _histarrays.Add( &_histsPar3 );
  _histarrays.Add( &_histsPar4 );
  
  _resulthistarray.Add( &_histErrorMean);
  _resulthistarray.Add( &_histErrorSigma);
  _resulthistarray.Add( &_histErrorChi2NDF);

  _histnames[0] = "Par1";
  _histnames[1] = "Par2";
  _histnames[2] = "Par3";
  _histnames[3] = "Par4";
  _results[0] = "Mean";
  _results[1] = "Sigma";
  _results[2] = "Chi2NDF";

  if (_verbose) cout << "---- The following Histograms are defined -----" << endl;

  //1D histograms
  TArrayD resolAxisBinning( _resolAxis->GetNbins()+1 );
  copyBinning( resolAxisBinning, *_resolAxis );

  for (int h = 0; h < _histarrays.GetEntries(); h++ ) {
    for (int Xindex = 0; Xindex <= _xAxis->GetNbins() + 1; Xindex++) {
      for (int Yindex = 0; Yindex <= _yAxis->GetNbins() + 1; Yindex++) {
	TString name = GetName() + (TString) "_" + _histnames[h] + _xAxis->GetName(); 
	name += Xindex; 
	name += _yAxis->GetName();
	name += Yindex;
	if ( Xindex == 0 || Yindex == 0 ) name += "_UNDERFLOW";
	if ( Xindex == _xAxis->GetNbins()+1 || Yindex == _yAxis->GetNbins()+1 ) name += "_OVERFLOW";

	TH1D* newhisto =  new TH1D( name, name, _resolAxis->GetNbins(), resolAxisBinning.GetArray() );
	newhisto->SetDirectory(0);
	((TObjArray*) _histarrays[h])->Add( newhisto );
	if (_verbose) cout << name << " , ";
      }
      if (_verbose) cout << endl;
    }
    if (_verbose) cout << endl;
  }
  if (_verbose) cout << endl;

  // 2D histograms
  TArrayD xAxisBinning( _xAxis->GetNbins()+1 );
  copyBinning( xAxisBinning, *_xAxis );
  TArrayD yAxisBinning( _yAxis->GetNbins()+1 );
  copyBinning( yAxisBinning, *_yAxis );

  for (int r = 0; r < _resulthistarray.GetEntries(); r++ ) {
    for (int h = 0; h < _nPar; h++ ) {	
      TString name = GetName() + (TString) "_" + _histnames[h] + _results[r];

      TH2D*  newhisto = new TH2D( name, name, 
				  _xAxis->GetNbins(), xAxisBinning.GetArray(), 
				  _yAxis->GetNbins(), yAxisBinning.GetArray() ) ;
      newhisto->SetDirectory(0);
      newhisto->GetXaxis()->SetTitle( _xAxis->GetTitle() );
      newhisto->GetYaxis()->SetTitle( _yAxis->GetTitle() );
      ((TObjArray*) _resulthistarray[r])->Add( newhisto );      
    }
    if (_verbose) cout << endl;
  } 
  if (_verbose)  {
    cout << endl << _xAxis->GetName() <<" defined from " << _xAxis->GetBinLowEdge(1) << " to " 
	 << _xAxis->GetBinUpEdge(_xAxis->GetNbins()) << " in " << _xAxis->GetNbins()
	 << " Bins ==> _BinWidth == "<< _xAxis->GetBinWidth(1) << endl;     
    cout << _yAxis->GetName() <<" defined from " << _yAxis->GetBinLowEdge(1) << " to " 
	 << _yAxis->GetBinUpEdge(_yAxis->GetNbins()) << " in " << _yAxis->GetNbins()
	 << " Bins ==> _BinWidth == "<< _yAxis->GetBinWidth(1) << endl;     
    cout <<"---- CreateHistos() finished ----" << endl;
  }

  // Correlation plots
  for ( Int_t x = 1; x <= 4; x++ ) {
    for ( Int_t y = 1; y <= 4; y++ ) {

      Int_t index = x + (y-1)*4 - 1;

      Double_t min = 3.*_resolAxis->GetBinLowEdge(1);
      Double_t max = 3.*_resolAxis->GetBinUpEdge(_resolAxis->GetNbins());
      Int_t nbins = 360;

      TString histname = GetName() + (TString) "_CorrMeas_"; 
      histname += x;
      histname += "_";
      histname += y;
      _histsCorrMeas[index] = new TH2D( histname, histname, nbins, min, max, nbins, min, max );
      _histsCorrMeas[index]->SetDirectory(0);
      _histsCorrMeas[index]->GetXaxis()->SetTitle( _histnames[x-1] );
      _histsCorrMeas[index]->GetYaxis()->SetTitle( _histnames[y-1] );

      histname = GetName() + (TString) "_CorrTrue_"; 
      histname += x;
      histname += "_";
      histname += y;
      _histsCorrTrue[index] = new TH2D( histname, histname, nbins, min, max, nbins, min, max );
      _histsCorrTrue[index]->SetDirectory(0);
      _histsCorrTrue[index]->GetXaxis()->SetTitle( _histnames[x-1] );
      _histsCorrTrue[index]->GetYaxis()->SetTitle( _histnames[y-1] );

    }
  }

}

void TGetErrorMatrix::Fill( TLorentzVector* TrueVec, TLorentzVector* MeasVec, Double_t XTrue,  Double_t XMeas,
			    Double_t YTrue,  Double_t YMeas) {
  // Fill true and measured parameter resolution in bins XTrue/YTrue amd XMeas/YMeas


  // update the particle's 4vec
  _theParticle->setIni4Vec(MeasVec);

  // get the true and meas parameters
  TMatrixD* TrueParams = _theParticle->transform(*TrueVec);
  const TMatrixD* MeasParams = _theParticle->getParIni();

//   // choose which histogram to fill
//   Int_t HistX = ( (Int_t) ((XMeas - _xlow) / _BinWidthX) ) + 1;
//   Int_t HistY = ( (Int_t) ((YMeas - _ylow) / _BinWidthY) ) + 1;
//   if ( XMeas == _xup ) HistX = _NbinsX;
//   if ( YMeas == _yup ) HistY = _NbinsY;

//   // overflow
//   if ( XMeas > _xup ) HistX = _NbinsX + 1;
//   if ( YMeas > _yup ) HistY = _NbinsY + 1;
  
//   // underflow
//   if ( XMeas < _xlow ) HistX = 0;
//   if ( YMeas < _ylow ) HistY = 0;


  // choose which histogram to fill
  Int_t HistX = _xAxis->FindBin( XMeas );
  Int_t HistY = _yAxis->FindBin( YMeas );

  for (int ParIndex = 0; ParIndex < TrueParams->GetNrows(); ParIndex++) {


    Double_t ParDiff =  (*MeasParams)(ParIndex,0) - (*TrueParams)(ParIndex,0);
    TH1D* HistToFill = Find1DHist(ParIndex, HistX, HistY);

//     if (ParIndex == 0 && ParDiff > 0.8) { 

//       cout << "PT Truth = " << (*TrueParams)(ParIndex,0) << ",  PT Reco = " << (*MeasParams)(ParIndex,0) << endl;

//     }


    //    cout << "Histogramm to fill: (" << HistX << "/" << HistY << "): " << HistToFill->GetName() << endl;
    if (HistToFill != NULL) {
      HistToFill -> Fill(ParDiff);
    } else {
      TString name = _histnames[ParIndex] + _xAxis->GetName() ;
      name += HistX;
      name += _yAxis->GetName();
      name += HistY;    
      cout << "HistoFill " << name <<" is NULL !!!!!!!!" << endl;
    }
    
    


    // Correlations
    for (int ParIndexY = 0; ParIndexY < TrueParams->GetNrows(); ParIndexY++) {
      TH2D* histcorrmeas = FindCorrMeasHist( ParIndex+1, ParIndexY+1 );
      histcorrmeas->Fill( (*MeasParams)(ParIndex,0), (*MeasParams)(ParIndexY,0) );
      histcorrmeas = FindCorrMeasHist( ParIndexY+1, ParIndex+1 );
      histcorrmeas->Fill( (*MeasParams)(ParIndexY,0), (*MeasParams)(ParIndex,0) );

      TH2D* histcorrtrue = FindCorrTrueHist( ParIndex+1, ParIndexY+1 );
      histcorrtrue->Fill( (*TrueParams)(ParIndex,0), (*TrueParams)(ParIndexY,0) );
      histcorrtrue = FindCorrTrueHist( ParIndexY+1, ParIndex+1 );
      histcorrtrue->Fill( (*TrueParams)(ParIndexY,0), (*TrueParams)(ParIndex,0) );
    }

  }
  delete TrueParams;

}

void TGetErrorMatrix::ResetFitParameters() {
  // reset the fit parameters

  for (int index=0; index < 3; index++ ) {
    _params[index] = -999.;
    _fitErrors[index] = -999.;
  }
  _chi2 = -999.;
  _ndf =  1;
}


void TGetErrorMatrix::Fit1DHisto(TH1D* HistToFit, Int_t MinBinEntry, const char* formula) {
  // Fit 1-dimensional histograms containing differences of 
  // true and measured parameters

  if (HistToFit->GetEntries() >= MinBinEntry) {
    
    // Do the fit
    TF1* func = (TF1*) gROOT->GetFunction(formula);
    func->SetParameters(HistToFit->GetMaximum(), HistToFit->GetMean(), HistToFit->GetRMS());
    HistToFit->Fit(func,"0","",HistToFit->GetMean() - 4.*(HistToFit->GetRMS()),
		   HistToFit->GetMean() + 4.*(HistToFit->GetRMS()) );
    TList* funclist =  HistToFit->GetListOfFunctions();
    func = (TF1*)funclist->At(0);	
    const Int_t kNotDraw = 1<<9;	    
    func ->ResetBit(kNotDraw);
    
    // Check whether fit makes sense
    if ( (func->GetParameter(1) > HistToFit->GetMean() - 8.*HistToFit->GetRMS()) &&
	 (func->GetParameter(1) < HistToFit->GetMean() + 8.*HistToFit->GetRMS()) ) {
      func -> GetParameters( &_params[0] );
      for (int i=0; i< func->GetNpar(); i++) {
	_fitErrors[i] = func -> GetParError(i);
      }
      _chi2 = func->GetChisquare();
      _ndf = func->GetNDF();
    }
    
  }
}  


void TGetErrorMatrix::Fill2DHist(TH2D* TheHist, int whichresult, Double_t X, Double_t Y) {
  // Fill 2-dimensional histograms containiing result mean(0), sigma(1), and chi2(2)

  Int_t bin = TheHist->FindBin( X, Y );
  if(whichresult ==_resulthistarray.GetEntries()-1) {
    TheHist->SetBinContent( bin, _chi2/_ndf); 
    TheHist->SetBinError(bin, 0);
  }
  else {
    TheHist->SetBinContent( bin, _params[whichresult+1]); 
    TheHist->SetBinError(bin, _fitErrors[whichresult+1]);
  }

}

void TGetErrorMatrix::Fit1DandFill2D() {
  // Fit all 1-dimensional histograms and fill results
  // into 2-dimensional histograms

  for (int h = 0; h < _nPar; h++ ) {
    
    //loop over hists in x
    for (int Xindex = 1; Xindex <= _xAxis->GetNbins(); Xindex++) {

      //loop over hists in y
      for (int Yindex = 1; Yindex <= _yAxis->GetNbins(); Yindex++) {

	ResetFitParameters();

	//get the histogram
	TH1D* HistToFit = Find1DHist(h, Xindex, Yindex);	
	if (HistToFit != NULL) {
	  Double_t X = _xAxis->GetBinCenter( Xindex );
	  Double_t Y = _yAxis->GetBinCenter( Yindex );

// 	  cout << "Fitting histogram (" << X << "," << Y << ")" << endl;
	  Fit1DHisto(HistToFit, 10, "gaus");
  
	  // Fill fit results mean, sigma chi2/ndf
	  if ( _chi2 >= 0. ) {
	    for (int res=0; res < _resulthistarray.GetEntries() ; res++) {
	      TH2D* TheHist = Find2DHist(h,res);
	      Fill2DHist(TheHist, res,  X, Y);
	    }
	  }
	} else {
	  TString name = _histnames[h] + _xAxis->GetName() ;
	  name += Xindex;
	  name += _yAxis->GetName();
	  name += Yindex;  
	  cout << "HistoFit " << name << " is NULL !!!!!" << endl;
	}
      }
    }
  }  
}

void TGetErrorMatrix::Fit2DResults(const char* formula) {

  for (int h = 0; h < _nPar; h++ ) {
    for (int res=0; res < _resulthistarray.GetEntries() ; res++) {
    
      TH2D* TheHist = Find2DHist(h,res);
      TheHist->Fit(formula);

    }
  }
}

void TGetErrorMatrix::setAxisTitles( TString& xTitle, TString& yTitle ) { 

  _xAxis->SetTitle( xTitle );
  _yAxis->SetTitle( yTitle );

  for (int h = 0; h < _histarrays.GetEntries(); h++ ) {
    for (int res=0; res < _resulthistarray.GetEntries() ; res++) {
      TH2D* TheHist = Find2DHist(h,res);
      TheHist->SetXTitle( xTitle );
      TheHist->SetYTitle( yTitle );
      TheHist->SetZTitle( _results[res] );
    }
  }

}

void TGetErrorMatrix::Plot(Bool_t ParamDistributions, Bool_t Results, 
			   const char* formula, Option_t* FuncDrawoption, Option_t* HistDrawOption ) {

  if (ParamDistributions) {   
    for (int h = 0; h < _nPar; h++ ) {
      const Float_t space = 0;  
      TString CanName = GetName() + (TString) "_" + _histnames[h] + "Canvas";
      TCanvas* Canvas = new TCanvas(CanName, CanName, 0, 0, 800, 1000);
      Canvas->Divide( _yAxis->GetNbins()+2, _xAxis->GetNbins()+2, space, space);
      Int_t canindex = 0;
      for (int Xindex = 0; Xindex < _xAxis->GetNbins()+2; Xindex++) {
	for (int Yindex = 0; Yindex < _yAxis->GetNbins()+2; Yindex++) {

	  canindex++;

	  //get histogram
	  TH1D* HistToPlot = Find1DHist(h, Xindex, Yindex); 

	  //plot
	  Canvas->cd(canindex);
	  HistToPlot->Draw();

	  // write binnnig
	  Double_t VarXLow = _xAxis->GetBinLowEdge(Xindex);
	  Double_t VarXUp  = _xAxis->GetBinUpEdge(Xindex);
	  Double_t VarYLow = _yAxis->GetBinLowEdge(Yindex);
	  Double_t VarYUp  = _yAxis->GetBinUpEdge(Yindex);
	  if (Xindex > 0 && Xindex <= _xAxis->GetNbins())
	    write( Form( "%3.1f <= %s < %3.1f", VarXLow, (const char*) _xAxis->GetTitle(), VarXUp ), 0.29, 0.9, 0.6, 1. );
	  else if (Xindex == 0)
	    write( Form( "%s < %3.1f", (const char*) _xAxis->GetTitle(), VarXUp ), 0.29, 0.9, 0.6, 1. );
	  else if (Xindex == _xAxis->GetNbins()+1)
	    write( Form( "%3.1f <= %s", VarXLow, (const char*)_xAxis->GetTitle() ), 0.29, 0.9, 0.6, 1. );
	  
	  if (Yindex > 0 && Yindex <= _yAxis->GetNbins())
	    write( Form( "%4.3f <= %s < %4.3f", VarYLow, (const char*) _yAxis->GetTitle(), VarYUp ) , 0.29, 0.82, 0.6, 0.92 );
	  else  if (Yindex == 0)
	    write( Form( "%s < %4.3f", (const char*) _yAxis->GetTitle(), VarYUp ) , 0.29, 0.82, 0.6, 0.92 );
	  else if (Yindex == _yAxis->GetNbins()+1)
	    write( Form( "%4.3f <= %s", VarYLow, (const char*) _yAxis->GetTitle() ), 0.29, 0.82, 0.6, 0.92 );
	  
	}
      }
    }
  }
  
  if (Results) {
    for (int h = 0; h < _nPar; h++ ) {
      TString CanName = GetName() + (TString) "_" + _histnames[h] + "ResultCanvas";
      TCanvas* ResultCanvas = new TCanvas(CanName, CanName,0,0,800,800);
      ResultCanvas -> Divide(2,2);
      for (int res=0; res < _resulthistarray.GetEntries() ; res++) {
	ResultCanvas->cd(res+1);
 	TH2D* TheHist = Find2DHist(h,res);
	TString formString =  formula;

	// Draw the histogram
	TheHist->Draw(HistDrawOption);
	TheHist->SetFillColor(38);	  
	TheHist->SetXTitle( _xAxis->GetTitle() );
	TheHist->SetYTitle( _yAxis->GetTitle() );
	TheHist->SetZTitle( _results[res] );

	// Draw the formula if specified
 	if (formString.Length() != 0) {
    	  TF2* fitpol3 = (TF2*) TheHist->FindObject(formula);
   	  TString Opt = "same " + (TString) FuncDrawoption;
	  fitpol3->Draw(Opt);
	}
      }
    }
  }

}

void TGetErrorMatrix::WriteTXT(TString Filename) {

  //  ofstream outfile( Filename, ios_base::app );  
  ofstream outfile( Filename );  
  outfile << Filename << endl;
  for (int h = 0; h < _histarrays.GetEntries(); h++ ) {

    outfile << endl 
	    <<"----------------------------"<< endl;
    outfile <<"   Parameter " << h+1 << " : " << endl;
    outfile <<"----------------------------"<< endl 
	    << endl;
    
    outfile << endl << "-- values for    :"<< _results[0] << "  " << _results[1] << "  " << _results[2] 
	    << " --" << endl << endl;  
    TH2D* TheHist1 = Find2DHist(h,0);
    TH2D* TheHist2 = Find2DHist(h,1);
    TH2D* TheHist3 = Find2DHist(h,2);	
    //LOOP OVER BINs
    for (int xbin = 0; xbin < TheHist1->GetNbinsX(); xbin++) {
      for (int ybin = 0; ybin < TheHist1->GetNbinsY(); ybin++) {
	outfile << _xAxis->GetName() << xbin << "  " << _yAxis->GetName() << ybin << " : " 
		<< TheHist1->GetBinContent(xbin,ybin)
		<< "  " <<  TheHist2->GetBinContent(xbin,ybin)
		<< "  " <<  TheHist3->GetBinContent(xbin,ybin)
		<< endl;
      }
    }	
  }

  outfile.close();
}

void TGetErrorMatrix::WriteROOT(TString Filename, TString opt) {
  // Write all histograms to root file

  TDirectory *dir = gDirectory;
  TFile *g = new TFile(Filename, opt);
  dir->GetList()->Write();    
  g->Close();

}

void TGetErrorMatrix::WritePS(Bool_t Gauss, Bool_t Results, TString Filename) {
  // Write plots to ps-file. You need to call Plot() before
  

  for (int h = 0; h < _nPar; h++ ) {
    if (Gauss) {
      TString CanName = GetName() + (TString) "_" + _histnames[h] + "Canvas";
      TCanvas* TheCanvas = (TCanvas*)gROOT->FindObject(CanName);
      if (h==0) TheCanvas->Print(Filename+"(");
      else TheCanvas->Print(Filename);
    }
    if (Results) {
      TString CanName = GetName() + (TString) "_" + _histnames[h] + "ResultCanvas";
      TCanvas* TheCanvas = (TCanvas*)gROOT->FindObject(CanName);
      if (h==_histarrays.GetEntries()-1)TheCanvas->Print(Filename+")");
      else TheCanvas->Print(Filename);
    }
  }

}

void TGetErrorMatrix::write( char* text, Double_t xlo,  Double_t ylo, Double_t xup, Double_t yup) {

  TPaveText* WriteSignal = new TPaveText(xlo, ylo, xup, yup, "NDC");
  WriteSignal->AddText(text);
  WriteSignal->SetFillColor(10);
  WriteSignal->SetBorderSize(0);
  WriteSignal->SetTextColor(kBlue);
  WriteSignal->SetTextSize(0.08);
  WriteSignal->Draw();

}

TH1D* TGetErrorMatrix::Find1DHist(Int_t ParIndex, Int_t X, Int_t Y) {
  // Returns the histogram corresponding to parameter 'ParIndex',
  // and bins 'X, Y'
  
  Int_t Hindex = X*(_yAxis->GetNbins()+2) + Y;
  TH1D* TheHist = (TH1D*)(((TObjArray*) _histarrays[ParIndex])->At(Hindex));
  return TheHist;

}

TH2D* TGetErrorMatrix::Find2DHist(Int_t ParIndex, Int_t result) {
  // Returns the histogram corresponding to parameter 'ParIndex'. 
  // 'result' means the histogram containing Mean (0), Sigma(1), or Chi2(2)

  TH2D* TheHist = (TH2D*)(((TObjArray*) _resulthistarray[result])->At(ParIndex));
  return TheHist;
}

TMatrixD* TGetErrorMatrix::getCovMatrix( Double_t X, Double_t Y ) {
  // Returns the covariance matrix corresponding to bin 'X, Y'

  TMatrixD* m = new TMatrixD( _nPar, _nPar );

  for (int parIndex = 0; parIndex < _nPar; parIndex++) {

    TH2D* hist = (TH2D*) ((TObjArray*) _resulthistarray[Sigma])->At(parIndex);
    Int_t bin = hist->FindBin( X, Y );
    Double_t sigma = hist->GetBinContent( bin );
    (*m)(parIndex, parIndex) = sigma*sigma;
    if ( sigma < 0. ) {
      if (_verbose) cout << "No proper error defined for parameter " << X << "/" << Y << "(index = " << parIndex << "): " << sigma << endl;
      delete m;
      return 0;
    }
  }
  
  return m;

}

TMatrixD* TGetErrorMatrix::getBiasCorrection( Double_t X, Double_t Y ) {
  // Returns (nPar, 1) matrix containing negative mean values of resolution distributions.
  // The returned matrix can be added to the measured particle parameters to correct the
  // bias.

  TMatrixD* m = new TMatrixD( _nPar, 1 );

  for (int parIndex = 0; parIndex < _nPar; parIndex++) {

    TH2D* hist = (TH2D*) ((TObjArray*) _resulthistarray[Mean])->At(parIndex);
    Int_t bin = hist->FindBin( X, Y );
    Double_t mean = hist->GetBinContent( bin );
    (*m)(parIndex, 0) = -1. * mean;
  }

  return m;

}

Double_t TGetErrorMatrix::getCorrMeas(  Int_t par1, Int_t par2 ) {
  // returns the correlation factor of parameter "par1" and "par2"
  // of measured quantities

  TH2D* hist = FindCorrMeasHist( par1, par2 );
  return (Double_t) hist->GetCorrelationFactor();
  
}

Double_t TGetErrorMatrix::getCorrTrue(  Int_t par1, Int_t par2 ) {
  // returns the correlation factor of parameter "par1" and "par2"
  // of true quantities

  TH2D* hist = FindCorrTrueHist( par1, par2 );
  return (Double_t) hist->GetCorrelationFactor();

}

TH2D* TGetErrorMatrix::FindCorrMeasHist( Int_t par1, Int_t par2 ) { 
  // returns 2D-plot par1 vs. par2 filled with true values

  if ( par1 >= 1 && par1 <= _nPar && par2 >= 1 && par2 <= _nPar ) {
    Int_t index = par1 + (par2-1)*4 - 1; 
    return _histsCorrMeas[index]; 
  } else {
    return 0;
  }

}

TH2D* TGetErrorMatrix::FindCorrTrueHist( Int_t par1, Int_t par2 ) {
  // returns 2D-plot par1 vs. par2 filled with measured values

  if ( par1 >= 1 && par1 <= _nPar && par2 >= 1 && par2 <= _nPar ) {
    Int_t index = par1 + (par2-1)*4 - 1;
    return _histsCorrTrue[ index ];
  } else {
    return 0;
  }

}

void TGetErrorMatrix::copyBinning( TArrayD& array, TAxis& axis ) {

  if ( array.GetSize() != axis.GetNbins() + 1 ) {
    cout << "Array has wrong size " << array.GetSize() << endl;
    return;
  }

  for ( Int_t index = 1; index <= axis.GetNbins()+1; index++) {
    array[index-1] = axis.GetBinLowEdge( index );
  }

}
