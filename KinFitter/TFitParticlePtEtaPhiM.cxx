// //Classname: TFitParticlePtEtaPhiM
// Author: Jan E. Sundermann, Verena Klose, Thomas Goepfert (TU Dresden)

//________________________________________________________________
// 
// TFitParticlePtEtaPhiM::
// --------------------
//
// Particle with spherical  parametrization of the momentum 4vector and
// free mass (4 free parameters). The parametrization is chosen as follows:
//
// p = (pt, eta, phi)
// E(fit) =  sqrt( |p|^2 + d^2 )
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticlePtEtaPhiM.h"
#include "TLorentzVector.h"

#include "TMath.h"

ClassImp(TFitParticlePtEtaPhiM)

//----------------
// Constructor --
//----------------
TFitParticlePtEtaPhiM::TFitParticlePtEtaPhiM()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticlePtEtaPhiM::TFitParticlePtEtaPhiM( const TFitParticlePtEtaPhiM& fitParticle )
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

TFitParticlePtEtaPhiM::TFitParticlePtEtaPhiM(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()
{ 
  init(pini, theCovMatrix);
}

TFitParticlePtEtaPhiM::TFitParticlePtEtaPhiM(const TString &name, const TString &title,
						   TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)
{ 
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticlePtEtaPhiM::clone( TString newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticlePtEtaPhiM( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticlePtEtaPhiM::~TFitParticlePtEtaPhiM() {

}

//--------------
// Operations --
//--------------
void TFitParticlePtEtaPhiM::init(TLorentzVector* pini, const TMatrixD* theCovMatrix) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticlePtEtaPhiM::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    cout << "Parameter matrix has wrong size." << endl;
    return 0;
  }
  if (_pini.Pt() == 0){
    cout << GetName() << endl;
    cout << _pini.Pt() << endl;
    cout << _pini.Eta() << endl;
    cout << _pini.Phi() << endl;
    cout << _pini.M() << endl;
  }

  Double_t pt = (*params)(0,0);
  Double_t eta = (*params)(1,0);
  Double_t phi = (*params)(2,0);
  Double_t m = (*params)(3,0);

  TLorentzVector* vec = new TLorentzVector( 0., 0., 0., 0. );
  vec->SetPtEtaPhiM( pt, eta, phi, m);
  return vec;

}


void TFitParticlePtEtaPhiM::setIni4Vec(const TLorentzVector* pini) {
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
    _iniparameters(1,0) = 1.;
    _iniparameters(2,0) = 1.;
    _iniparameters(3,0) = 1.;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters(0,0) = 1.;
    _parameters(1,0) = 1.;
    _parameters(2,0) = 1.;   
    _parameters(3,0) = 1.;   
    
  } else {

    Double_t pt = pini->Pt();
    Double_t eta = pini->Eta();
    Double_t phi = pini->Phi();
    Double_t m = pini->M();
    
    _pini.SetPtEtaPhiM(pt, eta, phi, m);
    _pcurr = _pini;
    
    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = pt;
    _iniparameters(1,0) = eta;
    _iniparameters(2,0) = phi;
    _iniparameters(3,0) = m;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

//     _u1.SetXYZ( TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta) );
//     _u2.SetXYZ( TMath::Cos(phi)*TMath::Cos(theta), TMath::Sin(phi)*TMath::Cos(theta), -1.*TMath::Sin(theta) );
//     _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );

  }

}

TMatrixD* TFitParticlePtEtaPhiM::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(pt, eta, phi, d) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,_nPar);
  (*DerivativeMatrix) *= 0.;

  if(_pini.Pt() == 0){
    cout << GetName() << endl;
    cout << _pini.Pt() << endl;
    cout << _pini.Eta() << endl;
    cout << _pini.Phi() << endl;
    cout << _pini.M() << endl;
  }

  Double_t pt = _parameters(0,0);
  Double_t eta = _parameters(1,0);
  Double_t phi = _parameters(2,0);
  Double_t m = _parameters(3,0);

  Double_t px = pt * TMath::Cos(phi);
  Double_t py = pt * TMath::Sin(phi);
  Double_t pz = pt * TMath::SinH(eta);

  Double_t p2 = pt*pt + pz*pz;
  Double_t e2 = p2 + m*m;
  Double_t e = TMath::Sqrt(e2);

  //1st column: dP/dpt
  (*DerivativeMatrix)(0,0) = px / pt;
  (*DerivativeMatrix)(1,0) = py / pt;
  (*DerivativeMatrix)(2,0) = pz / pt;
  (*DerivativeMatrix)(3,0) = pt * TMath::CosH(eta) * TMath::CosH(eta) / e;

  //2nd column: dP/deta
  (*DerivativeMatrix)(0,1) = 0.;
  (*DerivativeMatrix)(1,1) = 0.;
  (*DerivativeMatrix)(2,1) = pt * TMath::CosH(eta);
  (*DerivativeMatrix)(3,1) = pt * pt * TMath::SinH(eta) * TMath::CosH(eta) / e;

   //3rd column: dP/dphi
  (*DerivativeMatrix)(0,2) = -1. * pt * TMath::Sin(phi);
  (*DerivativeMatrix)(1,2) = pt * TMath::Cos(phi);
  (*DerivativeMatrix)(2,2) = 0.;
  (*DerivativeMatrix)(3,2) = 0.;

   //4th column: dP/d(d)
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = m / e;


  return DerivativeMatrix;

}

TMatrixD* TFitParticlePtEtaPhiM::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.Pt();
  (*tparams)(1,0) = vec.Eta();
  (*tparams)(2,0) = vec.Phi();
  (*tparams)(3,0) = vec.M();

  return tparams;

}
