#ifndef TFITPARTICLEPINVTHETAPHIE_H
#define TFITPARTICLEPINVTHETAPHIE_H

#include "TMatrixD.h"
#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFitParticlePInvThetaPhiE: public TAbsFitParticle {

public :

  TFitParticlePInvThetaPhiE();
  TFitParticlePInvThetaPhiE( const TFitParticlePInvThetaPhiE& fitParticle );
  TFitParticlePInvThetaPhiE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePInvThetaPhiE(const TString &name, const TString &title, 
			    TLorentzVector* pini, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePInvThetaPhiE();
  virtual TAbsFitParticle* clone( TString newname = "" ) const;
  
  // returns derivative dP/dy with P=(p,E) and y=(1/p, theta, phi, ...) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  virtual TLorentzVector* calc4Vec( const TMatrixD* params );
  
 protected :

  void init(TLorentzVector* pini, const TMatrixD* theCovMatrix );
  
  ClassDef(TFitParticlePInvThetaPhiE, 1)    // Particle with spherical 4vector parametrization (1/p, theta, phi) and free energy
  
private:
};

#endif
