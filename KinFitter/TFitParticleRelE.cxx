// //Classname: TFitParticleRelE
// Author: Thomas Goepfert (TU Dresden), Jan E. Sundermann (U Freiburg)

//________________________________________________________________
// 
// TFitParticleRelE::
// --------------------
//
// Particle with the relative energy change (a) as only free parameter. 
// The mass ist kept constant and the momentum scaled accordingly.
//
// p(fit) = SQRT[a**2 * E(fit)**2 - m**2 ]
// E(fit) = a * E(meas)
// Theta(fit) = Theta(meas)
// Phi(fit) = Phi(meas)
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticleRelE.h"

#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(TFitParticleRelE)

//----------------
// Constructor --
//----------------
TFitParticleRelE::TFitParticleRelE()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticleRelE::TFitParticleRelE( const TFitParticleRelE& fitParticle )
  :TAbsFitParticle( fitParticle.GetName(), fitParticle.GetTitle() )
{

  _nPar = fitParticle._nPar;
  _u1 = fitParticle._u1;
  _u2 = fitParticle._u2;
  _u3 = fitParticle._u3;
  _covMatrix.ResizeTo(  fitParticle._covMatrix );
  _covMatrix = fitParticle._covMatrix;
  _iniparameters.ResizeTo( fitParticle._iniparameters );
  _iniparameters = fitParticle._iniparameters;
  _parameters.ResizeTo( fitParticle._parameters );
  _parameters = fitParticle._parameters;
  _pini = fitParticle._pini;
  _pcurr = fitParticle._pcurr;

}

TFitParticleRelE::TFitParticleRelE(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()  
{
  init(pini, theCovMatrix);
}

TFitParticleRelE::TFitParticleRelE(const TString &name, const TString &title, 
			                       TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)  
{
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticleRelE::clone( TString newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticleRelE( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleRelE::~TFitParticleRelE() {

}

//--------------
// Operations --
//--------------
void TFitParticleRelE::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 1;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticleRelE::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    cout << "Parameter matrix has wrong size." << endl;
    return 0;
  }

  Double_t E = (*params)(0,0) * _pini.E();
  Double_t m = _pini.M();
  Double_t r = TMath::Sqrt( E*E - m*m );
  Double_t theta = _pini.Theta();
  Double_t phi = _pini.Phi();

  Double_t X = r*TMath::Cos(phi)*TMath::Sin(theta);
  Double_t Y = r*TMath::Sin(phi)*TMath::Sin(theta);
  Double_t Z = r*TMath::Cos(theta);
		
  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticleRelE::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values 

  if (pini == 0) {

    _u1.SetXYZ(0., 0., 0.);
    _u3.SetXYZ(0., 0., 0.);
    _u2.SetXYZ(0., 0., 0.);
    _pini.SetXYZT(0., 0., 0., 0.);
    _pcurr = _pini;

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 1.;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;
    
  } else {
    
//     Double_t r = pini->P();
    Double_t theta = pini->Theta();
    Double_t phi = pini->Phi();
//     Double_t E = pini->E();
    
    _pini = (*pini);
    _pcurr = _pini;
    
    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 1.;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _u1.SetXYZ( TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta) );
    _u2.SetXYZ( TMath::Cos(phi)*TMath::Cos(theta), TMath::Sin(phi)*TMath::Cos(theta), -1.*TMath::Sin(theta) );
    _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );

  }

}

TMatrixD* TFitParticleRelE::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(a) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/da).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,1);
  (*DerivativeMatrix) *= 0.;

  Double_t E = _parameters(0,0) * _pini.E();
  Double_t m = _pini.M();
  Double_t r = TMath::Sqrt( E*E - m*m );
  Double_t theta = _pini.Theta();
  Double_t phi = _pini.Phi();

  //1st column: dP/da
  (*DerivativeMatrix)(0,0) = _parameters(0,0) * _pini.E() * _pini.E() / r * TMath::Cos(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(1,0) = _parameters(0,0) * _pini.E() * _pini.E() / r * TMath::Sin(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(2,0) = _parameters(0,0) * _pini.E() * _pini.E() / r * TMath::Cos(theta);
  (*DerivativeMatrix)(3,0) = _pini.E();

  return DerivativeMatrix;

}

TMatrixD* TFitParticleRelE::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.E() / _pini.E();

  return tparams;

}
