// //Classname: TAbsFitConstraint
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)      


//________________________________________________________________
// 
// TAbsFitConstraint::
// --------------------
//
// Abstract base class for fit constraints
//

using namespace std;

#include "KinFitter/TAbsFitConstraint.h"
#include "TNamed.h"

#include <iostream>
#include "TClass.h"


ClassImp(TAbsFitConstraint)


TAbsFitConstraint::TAbsFitConstraint() 
  : TNamed("NoName","NoTitle")
  ,_covMatrix()
  ,_covMatrixFit()
  ,_covMatrixDeltaAlpha()
  ,_iniparameters()
  ,_parameters()

{
  _nPar = 0;
}

TAbsFitConstraint::TAbsFitConstraint(const TString &name, const TString &title) 
  : TNamed(name, title)
  ,_covMatrix()
  ,_covMatrixFit()
  ,_covMatrixDeltaAlpha()
  ,_iniparameters()
  ,_parameters()

{
  _nPar = 0;
}

TAbsFitConstraint::~TAbsFitConstraint() {

}

void TAbsFitConstraint::reset() {
  // Reset parameters to initial values

  if ( _nPar > 0 )
    _parameters = _iniparameters;  

  setCovMatrixFit( 0 );

}

void TAbsFitConstraint::setCovMatrix(const TMatrixD* theCovMatrix) {
  // Set measured alpha covariance matrix

  if ( _nPar == 0 )
    return;

  _covMatrix.ResizeTo(_nPar, _nPar);
  if( theCovMatrix == 0 ) {
    _covMatrix.Zero();
  } else if (theCovMatrix->GetNcols() ==_nPar && theCovMatrix->GetNrows() ==_nPar) {
    _covMatrix = (*theCovMatrix);
  } else {
    cout << GetName() 
	 << "::setCovMatrix()  Something is wrong with the definition of the covariance matrix" 
	 << endl;
  }

}

void TAbsFitConstraint::setCovMatrixFit(const TMatrixD* theCovMatrixFit) {
  // Set the fitted covariance matrix. Is called by the
  // kinematic fitter.

  if ( _nPar == 0 )
    return;

  _covMatrixFit.ResizeTo(_nPar, _nPar);
  if( theCovMatrixFit == 0 || _nPar == 0 ) {
    _covMatrixFit.Zero();
  } else if (theCovMatrixFit->GetNcols() ==_nPar && theCovMatrixFit->GetNrows() ==_nPar) {
    _covMatrixFit = (*theCovMatrixFit);
  } else {
    cout << GetName()
	 << "::setCovMatrixFit()  Something is wrong with the definition of the fit covariance matrix"
	 << endl;
  }

}

void TAbsFitConstraint::calcCovMatrixDeltaAlpha() {
  // Calculates V(deltaAlpha) ==  V(alpha_meas) - V(alpha_fit)

  if ( _nPar > 0 ) {
    _covMatrixDeltaAlpha.ResizeTo( _nPar, _nPar );
    _covMatrixDeltaAlpha = _covMatrix;
    if(_covMatrixFit.GetNrows() == _nPar && _covMatrixFit.GetNcols() == _nPar)
      _covMatrixDeltaAlpha -= _covMatrixFit;
    else 
      cout << GetName()
	   << "calcCovMatrixDeltaAlpha()::  _covMatrixFit probably not set" 
	   << endl;  
  }

}

void TAbsFitConstraint::applyDeltaAlpha(TMatrixD* corrMatrix) {
  // Apply corrections to the parameters wrt. to the
  // initial parameters alpha* = alpha + delta(alpha)

  if ( _nPar > 0 ) {
    _parameters = _iniparameters;
    _parameters += (*corrMatrix);
  }

}

void TAbsFitConstraint::setParIni(const TMatrixD* parini) {
  // Set initial parameter values (before the fit)

  if (parini == 0 || _nPar == 0) return;
  else if( parini->GetNrows() == _iniparameters.GetNrows() &&
	   parini->GetNcols() == _iniparameters.GetNcols() )
    _iniparameters = (*parini) ;
  else {
    cout << GetName() << "::setParIni()  Matrices don't fit" 
	 << endl;
    return;
  }

}

const TMatrixD* TAbsFitConstraint::getCovMatrix() const { 
  // Retrieve covariance matrix. Returns NULL if number of
  // parameters is 0.

  if ( _nPar > 0 ) {
    return &_covMatrix; 
  } else {
    return 0;
  }

}

const TMatrixD* TAbsFitConstraint::getCovMatrixFit() const { 
  // Retrieve covariance matrix of fit. Returns NULL if number of
  // parameters is 0.

  if ( _nPar > 0 ) {
    return &_covMatrixFit; 
  } else {
    return 0;
  }

}

const TMatrixD* TAbsFitConstraint::getCovMatrixDeltaAlpha() {
  // Returns covariance matrix 
  // V(deltaAlpha) ==  V(alpha_meas) - V(alpha_fit)

  calcCovMatrixDeltaAlpha();
  if (_nPar > 0) {
    return &_covMatrixDeltaAlpha;
  } else {
    return 0;
  }

}

const TMatrixD* TAbsFitConstraint::getParIni() { 
  // Retrieve initial parameters. Returns NULL if number of
  // parameters is 0.

  if (_nPar > 0) {
    return &_iniparameters; 
  } else {
    return 0;
  }

}

const TMatrixD* TAbsFitConstraint::getParCurr() { 
  // Retrieve current parameters. Returns NULL if number of
  // parameters is 0.

  if (_nPar > 0) {
    return &_parameters; 
  } else {
    return 0;
  }
    
}

void TAbsFitConstraint::print() {
  // Print information about constraint

  cout << "__________________________" << endl << endl;
  cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

  cout << "initial value: " << getInitValue() << endl;
  cout << "current value: " << getCurrentValue() << endl;

}

