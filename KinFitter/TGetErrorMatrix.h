
using namespace std;

#ifndef TGETERRORMATRIX_H
#define TGETERRORMATRIX_H


#include "TMatrixD.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TLorentzVector.h"


class TAbsFitParticle;

class TGetErrorMatrix : public TNamed {

public :

  enum ResultType {
    Mean = 0,
    Sigma = 1,
    Chi2 = 2
  };

  TGetErrorMatrix(const TString &name, const TString &title,
		  const TAbsFitParticle* theParticle, TString VarX, Int_t NbinsX, Double_t xlow, Double_t xup,
		  TString VarY, Int_t NbinsY, Double_t ylow, Double_t yup,
		  Int_t histResolNBins = 120, Double_t histResolMin = -3., Double_t histResolMax = 3. );
  TGetErrorMatrix(const TAbsFitParticle* theParticle, TString VarX, Int_t NbinsX, Double_t xlow, Double_t xup,
		  TString VarY, Int_t NbinsY, Double_t ylow, Double_t yup,
		  Int_t histResolNBins = 120, Double_t histResolMin = -3., Double_t histResolMax = 3. );
  TGetErrorMatrix( const TString &name, const TString &title,
		   const TAbsFitParticle* theParticle,
		   TAxis* xAxis, TAxis* yAxis, TAxis* resolAxis );
  TGetErrorMatrix();		   
  virtual ~TGetErrorMatrix();
  
  void Fill( TLorentzVector* TrueVec, TLorentzVector* MeasVec, Double_t XTrue,  Double_t XMeas,
	     Double_t YTrue,  Double_t YMeas);
  TMatrixD* getCovMatrix( Double_t X, Double_t Y );
  TMatrixD* getBiasCorrection( Double_t X, Double_t Y );

  void setAxisTitles( TString& xTitle, TString& yTitle );
  void Plot(Bool_t ParamDistributions=false, Bool_t Results=true, const char* formula="",
	    Option_t* FuncDrawoption="surf1", Option_t* HistDrawOption="p0");

  void Fit1DHisto(TH1D* HistToFit, Int_t MinBinEntry, const char* formula);
  void Fill2DHist(TH2D* TheHist, int whichresult, Double_t X, Double_t Y);
  void Fit1DandFill2D(); 
  void Fit2DResults(const char* formula);
  void WriteTXT(TString Filename);
  void WriteROOT(TString Filename, TString opt = "UPDATE");
  void WritePS(Bool_t Gauss, Bool_t Results, TString Filename);
  void write( char* text, Double_t xlo,  Double_t ylo, Double_t xup, Double_t yup);

  TH1D* Find1DHist(Int_t ParIndex, Int_t X, Int_t Y);
  TH2D* Find2DHist(Int_t ParIndex, Int_t result);

  TH2D* FindCorrMeasHist(Int_t par1, Int_t par2);
  TH2D* FindCorrTrueHist(Int_t par1, Int_t par2);

  Double_t getCorrMeas(Int_t par1, Int_t par2);
  Double_t getCorrTrue(Int_t par1, Int_t par2);

  const TAbsFitParticle* getParticle() const { return _theParticle; }
  Int_t getNPar() { return _nPar; }

  void CreateHistos();
  void DefineFitFuncs();
  void ResetFitParameters();

  const TAxis* getXaxis() { return _xAxis; }
  const TAxis* getYaxis() { return _yAxis; }

  // Resolution histogram binning
//   void setHistResolBinning( Int_t nBins, Double_t min, Double_t max ) { _histResolMin = min; _histResolMax = max; _histResolNBins = nBins; }
//   Double_t getHistResolMin() { return _histResolMin; }
//   Double_t getHistResolMax() { return _histResolMax; }
//   Int_t getHistResolNBins() { return _histResolNBins; }

protected:

  void copyBinning( TArrayD& array, TAxis& axis );     // Copies axis binning into array

private :

  TAbsFitParticle* _theParticle;   // Particle with 4vector parametrization (not saved in root-file due to problems with root3, root4 should work)
  Int_t _nPar;                     // Number of parameters of chosen parametrization
  
  TObjArray _histarrays;           // Contains ObjArrays _histsPar1 ... _histsPar4 with resolution histograms
  TObjArray _histsPar1;            // Contains resolution histograms for parameter 1
  TObjArray _histsPar2;            // Contains resolution histograms for parameter 2
  TObjArray _histsPar3;            // Contains resolution histograms for parameter 3
  TObjArray _histsPar4;            // Contains resolution histograms for parameter 4

  TH2D* _histsCorrMeas[16];        // Contains correlation plots of measured parameters 1..4
  TH2D* _histsCorrTrue[16];        // Contains correlation plots of true parameters 1..4
  
  TObjArray _resulthistarray;      // Contains results mean (0), sigma (1), and chi2 (2) of the resolutions distributions
  TObjArray _histErrorMean;        // Contains mean values of resolution distributions (TH2D) for the parameters 1 ... 4
  TObjArray _histErrorSigma;       // Contains sigma values of resolution distributions (TH2D) for the parameters 1 ... 4
  TObjArray _histErrorChi2NDF;     // Contains Chi2 values of resolution distributions (TH2D) for the parameters 1 ... 4

  TAxis* _xAxis;                   // Axis defining x-binning
  TAxis* _yAxis;                   // Axus defining y-binning
  TAxis* _resolAxis;               // Axis for resolution plots

//   TString _VarX;                   // Name of x-axis binning in 2D histograms
//   TString _VarY;                   // Name of y-axis binning in 2D histograms

//   Double_t _xlow;                  // Lower boundary of x-axis binning in 2D histograms
//   Double_t _xup;                   // Upper boundary of x-axis binning in 2D histograms
//   Int_t _NbinsX;                   // Number of bins of x-axis binning in 2D histograms
//   Double_t _BinWidthX;             // Bin width of x-axis binning in 2D histograms

//   Double_t _ylow;                  // Lower boundary of y-axis binning in 2D histograms
//   Double_t _yup;                   // Upper boundary of y-axis binning in 2D histograms
//   Int_t _NbinsY;                   // Number of bins of y-axis binning in 2D histograms
//   Double_t _BinWidthY;             // Bin width of x-axis binning in 2D histograms

//   TString _xAxisTitle;            // Title of x-axis in 2D plots
//   TString _yAxisTitle;            // Title of y-axis in 2D plots

//   Double_t _histResolMin;         // Lower boundary of resolution histograms
//   Double_t _histResolMax;         // Upper boundary of resolution histograms
//   Int_t _histResolNBins;          // Number of bins of resolution histograms

  Double_t _params[4];
  Double_t _fitErrors[4];

  Double_t _chi2;
  Double_t _ndf;

  TString _histnames[4];
  TString _results[3];   

  Bool_t _verbose;                   // Toggle verbose output

  ClassDef(TGetErrorMatrix, 1)
};

#endif
