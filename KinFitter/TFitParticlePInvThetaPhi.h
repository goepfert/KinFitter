#ifndef TFITPARTICLEPINVTHETAPHI_H
#define TFITPARTICLEPINVTHETAPHI_H

#include "TMatrixD.h"
#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFitParticlePInvThetaPhi: public TAbsFitParticle {

public :

  TFitParticlePInvThetaPhi();
  TFitParticlePInvThetaPhi( const TFitParticlePInvThetaPhi& fitParticle );
  TFitParticlePInvThetaPhi(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  TFitParticlePInvThetaPhi(const TString &name, const TString &title, 
			   TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePInvThetaPhi();
  virtual TAbsFitParticle* clone( TString newname = "" ) const;
  
  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, ...) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  void setIni4Vec(const TVector3* p, Double_t M);
  virtual TLorentzVector* calc4Vec( const TMatrixD* params );
  
 protected :
  
  void init(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  
  ClassDef(TFitParticlePInvThetaPhi, 1)    // Particle with spherical 4vector parametrization (1/r, theta, phi) and constrained mass
  
 private:
};

#endif
