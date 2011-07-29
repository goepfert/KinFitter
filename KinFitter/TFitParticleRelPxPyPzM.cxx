// //Classname: TFitParticleRelPxPyPzM
// Author: Jan E. Sundermann, Verena Klose, Thomas Goepfert (TU Dresden)      

//________________________________________________________________
// 
// TFitParticleRelPxPyPzM::
// --------------------
//
// Particle with cartesian 4vector parametrization and constant mass
// [four free parameters (px, py, pz, d) with initial values (1, 1, 1, 1)]
//
// p = px*u1 + py*u2 + pz*u3
// E = Sqrt( |p|^2 + (d*m)^2 )
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticleRelPxPyPzM.h"

#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(TFitParticleRelPxPyPzM)

//----------------
// Constructor --
//----------------
TFitParticleRelPxPyPzM::TFitParticleRelPxPyPzM()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticleRelPxPyPzM::TFitParticleRelPxPyPzM( const TFitParticleRelPxPyPzM& fitParticle )
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

TFitParticleRelPxPyPzM::TFitParticleRelPxPyPzM(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()  
{
  init(pini, theCovMatrix);
}

TFitParticleRelPxPyPzM::TFitParticleRelPxPyPzM(const TString &name, const TString &title, 
			   TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)  
{
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticleRelPxPyPzM::clone( TString newname ) const {
  // Returns a copy of itself

  TAbsFitParticle* myclone = new TFitParticleRelPxPyPzM( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleRelPxPyPzM::~TFitParticleRelPxPyPzM() {

}

//--------------
// Operations --
//--------------
void TFitParticleRelPxPyPzM::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticleRelPxPyPzM::calc4Vec( const TMatrixD* params ) {
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
  Double_t M = (*params)(3,0) * _pini.M();

  Double_t E =  TMath::Sqrt( X*X + Y*Y + Z*Z + M*M );

  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticleRelPxPyPzM::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values

  if (pini == 0) {

    _u1.SetXYZ(0., 0., 0.);
    _u3.SetXYZ(0., 0., 0.);
    _u2.SetXYZ(0., 0., 0.);
    _pini.SetXYZT(0., 0., 0., 0.);
    _pcurr = _pini;

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0)=0;
    _iniparameters(1,0)=0;
    _iniparameters(2,0)=0;
    _iniparameters(3,0)=0.;

    _parameters.ResizeTo(_nPar,1);
    _parameters(0,0)=0.;
    _parameters(1,0)=0.;
    _parameters(2,0)=0.;   
    _parameters(3,0)=0.;   
    
  } else {
    
    _pini = (*pini);
    _pcurr = _pini;

    _u1.SetXYZ( 1., 0., 0. );
    _u2.SetXYZ( 0., 1., 0. );
    _u3.SetXYZ( 0., 0., 1. );
  
    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 1.;
    _iniparameters(1,0) = 1.;
    _iniparameters(2,0) = 1.;
    _iniparameters(3,0) = 1.;
    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

  }

}

TMatrixD* TFitParticleRelPxPyPzM::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(a, b, c, d) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dpx, dP/dpy, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,4);
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
  Double_t d = _parameters(3,0);
  Double_t mcur = d * m;

  Double_t ecur2 = pxcur*pxcur + pycur*pycur + pzcur*pzcur + mcur*mcur;
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

   //4th column: dP/dd
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = mcur * m / ecur;

  return DerivativeMatrix;  
}

TMatrixD* TFitParticleRelPxPyPzM::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.X()/_pini.Px();
  (*tparams)(1,0) = vec.Y()/_pini.Py();
  (*tparams)(2,0) = vec.Z()/_pini.Pz();
  (*tparams)(3,0) = vec.M()/_pini.M();

  return tparams;

}
