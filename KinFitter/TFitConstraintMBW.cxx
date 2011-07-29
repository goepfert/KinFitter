// //Classname: TFitConstraintMBW
// Author: Thomas Goepfert (TU Dresden), Jan E. Sundermann (U Freiburg)


//________________________________________________________________
// 
// TFitConstraintMBW::
// --------------------
//
// Fit constraint: mass conservation ( m_i - m_j - alpha == 0 )
//                 alpha: Breit-Wigner
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitConstraintMBW.h"
#include "KinFitter/TFitConstraintM.h"

#include "TLorentzVector.h"
#include "TClass.h"

ClassImp(TFitConstraintMBW)

//----------------
// Constructor --
//----------------
TFitConstraintMBW::TFitConstraintMBW()
  : TFitConstraintM() 
{

  init();

}

TFitConstraintMBW::TFitConstraintMBW(vector<TAbsFitParticle*>* ParList1,
				     vector<TAbsFitParticle*>* ParList2, 
				     Double_t Mass,
				     Double_t Width)
  : TFitConstraintM(ParList1, ParList2, Mass ) 
{
  
  init();
  setMassConstraint( Mass, Width );

}

TFitConstraintMBW::TFitConstraintMBW(const TString &name, const TString &title,
				     vector<TAbsFitParticle*>* ParList1,
				     vector<TAbsFitParticle*>* ParList2, 
				     Double_t Mass,
				     Double_t Width)
  : TFitConstraintM( name, title, ParList1, ParList2, Mass )
{

  init();
  setMassConstraint( Mass, Width );
  
}

void
TFitConstraintMBW::init() {

  _nPar = 1;
  _iniparameters.ResizeTo(1,1);
  _iniparameters(0,0) = 1.;
  _parameters.ResizeTo(1,1);
  _parameters = _iniparameters;
  setCovMatrix( 0 );
  _covMatrix(0,0) = 1.;

}


//--------------
// Destructor --
//--------------
TFitConstraintMBW::~TFitConstraintMBW() {

}

//--------------
// Operations --
//--------------

void TFitConstraintMBW::setMassConstraint(Double_t Mass, Double_t Width) { 
  
  // Set errors
  _width = Width;
  _sigma = _width / 2.; //arbitrary choice
  setCovMatrix( 0 );
  _covMatrix(0,0) = _sigma*_sigma;


  // Initialize mass
  _TheMassConstraint = Mass;
//   _iniparameters(0,0) = Mass;
  _iniparameters(0,0) = transformBWToGaus( Mass );
  _parameters(0,0) = _iniparameters(0,0);

}

Double_t TFitConstraintMBW::transformBWToGaus( Double_t M ) {

  Double_t M0 = _TheMassConstraint;

  Double_t mTransformed = M0 + TMath::Sqrt(2.) * _sigma * 
    TMath::ErfInverse( 2. / TMath::Pi() * TMath::ATan( (M - M0) / _sigma ) );
    
  return mTransformed;

}

Double_t TFitConstraintMBW::transformGausToBW( Double_t M ) {

  Double_t M0 = _TheMassConstraint;

  Double_t mTransformed = M0 + 
    _width/2. * TMath::Tan( TMath::Pi()/2. * TMath::Erf( (M-M0)/TMath::Sqrt(2.)/_sigma ) );
  
  return mTransformed;

}

Double_t TFitConstraintMBW::getInitValue() {
  // Get initial value of constraint (before the fit)

  Double_t InitValue = 
    CalcMass( &_ParList1, true ) - 
    CalcMass( &_ParList2, true ) - 
//     transformBWToGaus( _iniparameters(0,0) );
    transformGausToBW( _iniparameters(0,0) );

  return InitValue;

}

Double_t TFitConstraintMBW::getCurrentValue() {
  // Get value of constraint after the fit

  Double_t CurrentValue =
    CalcMass(&_ParList1,false) - 
    CalcMass(&_ParList2,false) - 
//     transformBWToGaus( _parameters(0,0) );
    transformGausToBW( _parameters(0,0) );

  return CurrentValue;

}

TMatrixD* TFitConstraintMBW::getDerivativeAlpha() {
  // Calculate dF/dAlpha

  TMatrixD* DerivativeMatrix = new TMatrixD(1,1);
  DerivativeMatrix->Zero();

  Double_t M0 = _TheMassConstraint;
  Double_t M = _parameters(0,0);
  Double_t aDerivative = 0.;

//   // BWToGaus
//   Double_t z = (M - M0) / _sigma;
//   Double_t InvErf = TMath::ErfInverse( 2. / TMath::Pi() * TMath::ATan(z) );
//   aDerivative = 
//     -1. * 
//     TMath::Sqrt( 2. / TMath::Pi() ) * 
//     TMath::Exp( InvErf*InvErf ) / 
//     1 / ( z*z + 1 );

  // GausToBW
  Double_t a = (M - M0) / TMath::Sqrt(2.) / _sigma;
  Double_t b = TMath::Pi()/2. * TMath::Erf(a);
  aDerivative = -1. *
    _width / 2. / _sigma *
    1. / (TMath::Cos(b)*TMath::Cos(b)) *
    TMath::Sqrt(TMath::Pi()/2.) *
    TMath::Exp(-1.*a*a);

  (*DerivativeMatrix)(0,0) = aDerivative;

  return DerivativeMatrix;

}

void TFitConstraintMBW::print() {

  cout << "__________________________" << endl << endl;
  cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

  cout << "initial value: " << getInitValue() << endl;
  cout << "current value: " << getCurrentValue() << endl;
  cout << "mean mass: " << _TheMassConstraint << endl;
  cout << "width: " << _width << endl;
//   cout << "initial mass: " << _iniparameters(0,0)*_TheMassConstraint  << endl;
//   cout << "current mass: " << _parameters(0,0)*_TheMassConstraint  << endl;
//   cout << "initial mass: " << _iniparameters(0,0)  << endl;
//   cout << "current mass: " << _parameters(0,0)  << endl;
  cout << "initial mass: " << transformGausToBW( _iniparameters(0,0) ) << endl;
  cout << "current mass: " << transformGausToBW( _parameters(0,0) ) << endl;


}
