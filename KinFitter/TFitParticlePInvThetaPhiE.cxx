// //Classname: TFitParticlePInvThetaPhiE
// Author: Verena Klose, Jan E. Sundermann (TU Dresden)

//________________________________________________________________
// 
// TFitParticlePInvThetaPhiE::
// --------------------
//
// Particle with 1/p, theta, phi, and E parametrization of the momentum 4vector
// (4 free parameters). The parametrization is chosen as follows:
//
// p_vec = (1/p, theta, phi)
// E(fit) =  E(meas)
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticlePInvThetaPhiE.h"

#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(TFitParticlePInvThetaPhiE)

//----------------
// Constructor --
//----------------
TFitParticlePInvThetaPhiE::TFitParticlePInvThetaPhiE()
  :TAbsFitParticle()
{ 
  init( 0, 0);
}

TFitParticlePInvThetaPhiE::TFitParticlePInvThetaPhiE( const TFitParticlePInvThetaPhiE& fitParticle )
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

TFitParticlePInvThetaPhiE::TFitParticlePInvThetaPhiE(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()
{ 
  init(pini, theCovMatrix);
}

TFitParticlePInvThetaPhiE::TFitParticlePInvThetaPhiE(const TString &name, const TString &title,
					       TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)
{ 
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticlePInvThetaPhiE::clone( TString newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticlePInvThetaPhiE( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticlePInvThetaPhiE::~TFitParticlePInvThetaPhiE() {

}


//--------------
// Operations --
//--------------
void TFitParticlePInvThetaPhiE::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);
  
}
 
TLorentzVector* TFitParticlePInvThetaPhiE::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    cout << "Parameter matrix has wrong size." << endl;
    return 0;
  }
  
  Double_t r = (*params)(0,0);
  Double_t theta = (*params)(1,0);
  Double_t phi = (*params)(2,0);
  Double_t E = (*params)(3,0);

  Double_t X = 1./r*TMath::Cos(phi)*TMath::Sin(theta);
  Double_t Y = 1./r*TMath::Sin(phi)*TMath::Sin(theta);
  Double_t Z = 1./r*TMath::Cos(theta);

  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticlePInvThetaPhiE::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values 

  if ( pini == 0 ) {

    _u1.SetXYZ(0., 0., 0.);
    _u3.SetXYZ(0., 0., 0.);
    _u2.SetXYZ(0., 0., 0.);
    _pini.SetXYZT(0., 0., 0., 0.);
    _pcurr = _pini;

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 0.;
    _iniparameters(1,0) = 0.;
    _iniparameters(2,0) = 0.;
    _iniparameters(3,0) = 0.;

    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

  } else {

    _pini = (*pini);
    _pcurr = _pini;

    Double_t r = 1/_pini.P();
    Double_t theta = _pini.Theta();
    Double_t phi = _pini.Phi();
    Double_t E = _pini.E();

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = r;
    _iniparameters(1,0) = theta;
    _iniparameters(2,0) = phi;
    _iniparameters(3,0) = E;

    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _u1.SetXYZ( TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta) );
    _u2.SetXYZ( TMath::Cos(phi)*TMath::Cos(theta), TMath::Sin(phi)*TMath::Cos(theta), -1.*TMath::Sin(theta) );
    _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );
  
  }

}

TMatrixD* TFitParticlePInvThetaPhiE::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, E) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,4);
  (*DerivativeMatrix) *= 0.;

  Double_t r = _parameters(0,0);
  Double_t p = 1./r;
  Double_t theta = _parameters(1,0);
  Double_t phi = _parameters(2,0);

  //1st column: dP/d(1/p)
  (*DerivativeMatrix)(0,0) = -1.*p*p*TMath::Cos(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(1,0) = -1.*p*p*TMath::Sin(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(2,0) = -1.*p*p*TMath::Cos(theta);
  (*DerivativeMatrix)(3,0) = 0.;

  //2nd column: dP/dtheta
  (*DerivativeMatrix)(0,1) = p*TMath::Cos(phi)*TMath::Cos(theta);
  (*DerivativeMatrix)(1,1) = p*TMath::Sin(phi)*TMath::Cos(theta);
  (*DerivativeMatrix)(2,1) = -1.*p*TMath::Sin(theta);
  (*DerivativeMatrix)(3,1) = 0.;

   //3rd column: dP/dphi
  (*DerivativeMatrix)(0,2) = -1.*p*TMath::Sin(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(1,2) = p*TMath::Cos(phi)*TMath::Sin(theta);;
  (*DerivativeMatrix)(2,2) = 0.;
  (*DerivativeMatrix)(3,2) = 0.;

   //4rd column: dP/dE
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = 1.;

  return DerivativeMatrix;  

}

TMatrixD* TFitParticlePInvThetaPhiE::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = 1./vec.P();
  (*tparams)(1,0) = vec.Theta();
  (*tparams)(2,0) = vec.Phi();
  (*tparams)(3,0) = vec.E();

  return tparams;

}
