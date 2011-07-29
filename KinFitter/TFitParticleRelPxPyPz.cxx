// //Classname: TFitParticleRelPxPyPz
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)      

//________________________________________________________________
// 
// TFitParticleRelPxPyPz::
// --------------------
//
// Particle with cartesian 4vector parametrization and constrained mass
// [three free parameters (px, py, pz) with initial values (1, 1, 1)]
//
// p = px*u1 + py*u2 + pz*u3
// E = Sqrt( |p|^2 + m^2 )
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticleRelPxPyPz.h"

#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(TFitParticleRelPxPyPz)

//----------------
// Constructor --
//----------------
TFitParticleRelPxPyPz::TFitParticleRelPxPyPz()
  :TAbsFitParticle()
{ 
  init( 0, 0., 0);
}

TFitParticleRelPxPyPz::TFitParticleRelPxPyPz( const TFitParticleRelPxPyPz& fitParticle )
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

TFitParticleRelPxPyPz::TFitParticleRelPxPyPz(TVector3* p, Double_t M, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()
{ 
  init(p, M, theCovMatrix);
}

TFitParticleRelPxPyPz::TFitParticleRelPxPyPz(const TString &name, const TString &title,
					     TVector3* p, Double_t M, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)
{ 
  init(p, M, theCovMatrix);
}

TAbsFitParticle* TFitParticleRelPxPyPz::clone( TString newname ) const {
  // Returns a copy of itself

  TAbsFitParticle* myclone = new TFitParticleRelPxPyPz( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleRelPxPyPz::~TFitParticleRelPxPyPz() {

}


//--------------
// Operations --
//--------------
void TFitParticleRelPxPyPz::init(TVector3* p, Double_t M, const TMatrixD* theCovMatrix) {

  _nPar = 3;
  setIni4Vec(p, M);
  setCovMatrix(theCovMatrix);
  
}


TLorentzVector* TFitParticleRelPxPyPz::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    cout << "Parameter matrix has wrong size." << endl;
    return 0;
  }

  Double_t X = (*params)(0,0) * _pini.Px();
  Double_t Y = (*params)(1,0) * _pini.Py();
  Double_t Z = (*params)(2,0) * _pini.Pz();
  Double_t E = TMath::Sqrt(  X*X + Y*Y + Z*Z + _pini.M2() );

  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticleRelPxPyPz::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values 

  TVector3 vec( pini->Vect() );
  setIni4Vec( &vec, pini->M() );

}

void TFitParticleRelPxPyPz::setIni4Vec(const TVector3* p, Double_t M) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values

  if ( p == 0 ) {

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 0.;
    _iniparameters(1,0) = 0.;
    _iniparameters(2,0) = 0.;
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _pini.SetXYZM( 0., 0., 0., M);
    _pcurr = _pini;

  } else {

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 1.;
    _iniparameters(1,0) = 1.;
    _iniparameters(2,0) = 1.;
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _pini.SetXYZM( p->x(), p->y(), p->z(), M);
    _pcurr = _pini;

    _u1.SetXYZ( 1., 0., 0. );
    _u2.SetXYZ( 0., 1., 0. );
    _u3.SetXYZ( 0., 0., 1. );
  
  }

}

TMatrixD* TFitParticleRelPxPyPz::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(px, py, pz) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dpx, dP/dpy, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,3);
  (*DerivativeMatrix) *= 0.;

  Double_t px = _pini.Px();
  Double_t a = _parameters(0,0);
  Double_t pxcur = a * px;

  Double_t py = _pini.Py();
  Double_t b = _parameters(1,0);
  Double_t pycur = b * py;

  Double_t pz = _pini.Pz();
  Double_t c = _parameters(2,0);
  Double_t pzcur = c * pz;

  Double_t m = _pini.M();

  Double_t ecur2 = pxcur*pxcur + pycur*pycur + pzcur*pzcur + m*m;
  Double_t ecur = TMath::Sqrt(ecur2);

  //1st column: dP/da
  (*DerivativeMatrix)(0,0) = px;
  (*DerivativeMatrix)(1,0) = 0.;
  (*DerivativeMatrix)(2,0) = 0.;
  (*DerivativeMatrix)(3,0) = pxcur * px / ecur;

  //2nd column: dP/db
  (*DerivativeMatrix)(0,1) = 0.;
  (*DerivativeMatrix)(1,1) = py;
  (*DerivativeMatrix)(2,1) = 0.;
  (*DerivativeMatrix)(3,1) = pycur * py / ecur;

   //3rd column: dP/dc
  (*DerivativeMatrix)(0,2) = 0.;
  (*DerivativeMatrix)(1,2) = 0.;
  (*DerivativeMatrix)(2,2) = pz;
  (*DerivativeMatrix)(3,2) = pzcur * pz / ecur;

  return DerivativeMatrix;  

}

TMatrixD* TFitParticleRelPxPyPz::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  TVector3 vec3( vec.Vect() );

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec3.x()/_pini.Px();
  (*tparams)(1,0) = vec3.y()/_pini.Py();
  (*tparams)(2,0) = vec3.z()/_pini.Pz();

  return tparams;

}
