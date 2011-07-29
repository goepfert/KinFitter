// //Classname: TFitParticleRelPtEtaPhiM
// Author: Jan E. Sundermann, Verena Klose, Thomas Goepfert (TU Dresden)      

//________________________________________________________________
// 
// TFitParticleRelPtEtaPhiM::
// --------------------
//
// Particle with spherical  parametrization of the momentum 4vector and
// free mass (4 free parameters). The parametrization is chosen as follows:
//
// p = (a*pt, b*eta, c*phi)
// E(fit) =  sqrt( |p|^2 + (d*m)^2 )
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticleRelPtEtaPhiM.h"
#include "TLorentzVector.h"

#include "TMath.h"

ClassImp(TFitParticleRelPtEtaPhiM)

//----------------
// Constructor --
//----------------
TFitParticleRelPtEtaPhiM::TFitParticleRelPtEtaPhiM()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticleRelPtEtaPhiM::TFitParticleRelPtEtaPhiM( const TFitParticleRelPtEtaPhiM& fitParticle )
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

TFitParticleRelPtEtaPhiM::TFitParticleRelPtEtaPhiM(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()
{ 
  init(pini, theCovMatrix);
}

TFitParticleRelPtEtaPhiM::TFitParticleRelPtEtaPhiM(const TString &name, const TString &title,
						   TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)
{ 
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticleRelPtEtaPhiM::clone( TString newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticleRelPtEtaPhiM( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleRelPtEtaPhiM::~TFitParticleRelPtEtaPhiM() {

}

//--------------
// Operations --
//--------------
void TFitParticleRelPtEtaPhiM::init(TLorentzVector* pini, const TMatrixD* theCovMatrix) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticleRelPtEtaPhiM::calc4Vec( const TMatrixD* params ) {
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

  Double_t pt = (*params)(0,0) * _pini.Pt();
  Double_t eta = (*params)(1,0) * _pini.Eta();
  Double_t phi = (*params)(2,0) * _pini.Phi();
  Double_t m = (*params)(3,0) * _pini.M();

  TLorentzVector* vec = new TLorentzVector( 0., 0., 0., 0. );
  vec->SetPtEtaPhiM( pt, eta, phi, m);
  return vec;

}


void TFitParticleRelPtEtaPhiM::setIni4Vec(const TLorentzVector* pini) {
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
    _iniparameters(0,0) = 1.;
    _iniparameters(1,0) = 1.;
    _iniparameters(2,0) = 1.;
    _iniparameters(2,0) = 1.;
    
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

//     _u1.SetXYZ( TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta) );
//     _u2.SetXYZ( TMath::Cos(phi)*TMath::Cos(theta), TMath::Sin(phi)*TMath::Cos(theta), -1.*TMath::Sin(theta) );
//     _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );

  }

}

TMatrixD* TFitParticleRelPtEtaPhiM::getDerivative() {
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
    cout << _pini.E() << endl;
  }

  Double_t pt = _pini.Pt();
  Double_t a = _parameters(0,0);
  Double_t ptcur = a * pt;

  Double_t eta = _pini.Eta();
  Double_t b = _parameters(1,0);
  Double_t etacur = b * eta;

  Double_t phi = _pini.Phi();
  Double_t c = _parameters(2,0);
  Double_t phicur = c * phi;

  Double_t m = _pini.M();
  Double_t d = _parameters(3,0);
  Double_t mcur = d * m;

  Double_t pxcur = ptcur * TMath::Cos(phicur);
  Double_t pycur = ptcur * TMath::Sin(phicur);
  Double_t pzcur = ptcur * TMath::SinH(etacur);

  Double_t pcur2 = ptcur*ptcur + pzcur*pzcur;
  Double_t ecur2 = pcur2 + mcur*mcur;
  Double_t ecur = TMath::Sqrt(ecur2);

  //1st column: dP/da
  (*DerivativeMatrix)(0,0) = pxcur / a;
  (*DerivativeMatrix)(1,0) = pycur / a;
  (*DerivativeMatrix)(2,0) = pzcur / a;
  (*DerivativeMatrix)(3,0) = ptcur * pt * TMath::CosH(etacur) * TMath::CosH(etacur) / ecur;

  //2nd column: dP/db
  (*DerivativeMatrix)(0,1) = 0.;
  (*DerivativeMatrix)(1,1) = 0.;
  (*DerivativeMatrix)(2,1) = ptcur * eta * TMath::CosH(etacur);
  (*DerivativeMatrix)(3,1) = ptcur * ptcur * eta * TMath::SinH(etacur) * TMath::CosH(etacur) / ecur;

   //3rd column: dP/dc
  (*DerivativeMatrix)(0,2) = -1. * ptcur * phi * TMath::Sin(phicur);
  (*DerivativeMatrix)(1,2) = ptcur * phi * TMath::Cos(phicur);
  (*DerivativeMatrix)(2,2) = 0.;
  (*DerivativeMatrix)(3,2) = 0.;

   //4th column: dP/d(d)
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = mcur * m / ecur;

  return DerivativeMatrix;

}

TMatrixD* TFitParticleRelPtEtaPhiM::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.Pt()/_pini.Pt();
  (*tparams)(1,0) = vec.Eta()/_pini.Eta();
  (*tparams)(2,0) = vec.Phi()/_pini.Phi();
  (*tparams)(3,0) = vec.M()/_pini.M();

  return tparams;

}
