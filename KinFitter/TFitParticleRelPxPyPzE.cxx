// //Classname: TFitParticleRelPxPyPzE
// Author: Verena Klose (TU Dresden), Jan E. Sundermann

//________________________________________________________________
// 
// TFitParticleRelPxPyPzE::
// --------------------
//
// Particle with cartesian 4vector parametrization and free energy
// [four free parameters (a, b, c, d) with initial values (1, 1, 1, 1)]
//
// p_fit = a*px_meas*u_x + b*py_meas*u_y + c*pz_meas*u_z
// E_fit = d*E_meas
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticleRelPxPyPzE.h"

#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(TFitParticleRelPxPyPzE)

//----------------
// Constructor --
//----------------
TFitParticleRelPxPyPzE::TFitParticleRelPxPyPzE()
  :TAbsFitParticle()  
{
  init(0, 0);
}

TFitParticleRelPxPyPzE::TFitParticleRelPxPyPzE( const TFitParticleRelPxPyPzE& fitParticle )
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

TFitParticleRelPxPyPzE::TFitParticleRelPxPyPzE(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()  
{
  init(pini, theCovMatrix);
}

TFitParticleRelPxPyPzE::TFitParticleRelPxPyPzE(const TString &name, const TString &title, 
			   TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)  
{
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticleRelPxPyPzE::clone( TString newname ) const {
  // Returns a copy of itself

  TAbsFitParticle* myclone = new TFitParticleRelPxPyPzE( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;

}

//--------------
// Destructor --
//--------------
TFitParticleRelPxPyPzE::~TFitParticleRelPxPyPzE() {

}

//--------------
// Operations --
//--------------
void TFitParticleRelPxPyPzE::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);

}

TLorentzVector* TFitParticleRelPxPyPzE::calc4Vec( const TMatrixD* params ) {
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
  Double_t E = (*params)(3,0) * _pini.E();

  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticleRelPxPyPzE::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values

  if (pini == 0) {
    _pini.SetXYZT(0., 0., 0., 0.);
  } else {
    _pini = (*pini);
  }

  _pcurr = _pini;

  _u1.SetXYZ( 1., 0., 0. );
  _u2.SetXYZ( 0., 1., 0. );
  _u3.SetXYZ( 0., 0., 1. );

  _iniparameters.ResizeTo(_nPar,1);
  _parameters.ResizeTo(_nPar,1);

  _iniparameters(0,0) = 1.;
  _iniparameters(1,0) = 1.;
  _iniparameters(2,0) = 1.;
  _iniparameters(3,0) = 1.;
  _parameters = _iniparameters;

}

TMatrixD* TFitParticleRelPxPyPzE::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(a, b, c, d) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/da, dP/db, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,4);
  (*DerivativeMatrix) *= 0.;

  //1st column: dP/dx
  (*DerivativeMatrix)(0,0) = _pini.Px();
  (*DerivativeMatrix)(1,0) = 0.;
  (*DerivativeMatrix)(2,0) = 0.;
  (*DerivativeMatrix)(3,0) = 0.;

  //2nd column: dP/dy
  (*DerivativeMatrix)(0,1) = 0.;
  (*DerivativeMatrix)(1,1) = _pini.Py();
  (*DerivativeMatrix)(2,1) = 0.;
  (*DerivativeMatrix)(3,1) = 0.;

   //3rd column: dP/dz
  (*DerivativeMatrix)(0,2) = 0.;
  (*DerivativeMatrix)(1,2) = 0.;
  (*DerivativeMatrix)(2,2) = _pini.Pz();
  (*DerivativeMatrix)(3,2) = 0.;

   //4th column: dP/dd
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = _pini.E();

  return DerivativeMatrix;  
}

TMatrixD* TFitParticleRelPxPyPzE::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = vec.Px()/_pini.Px();
  (*tparams)(1,0) = vec.Py()/_pini.Px();
  (*tparams)(2,0) = vec.Pz()/_pini.Px();
  (*tparams)(3,0) = vec.E()/_pini.E();

  return tparams;

}
