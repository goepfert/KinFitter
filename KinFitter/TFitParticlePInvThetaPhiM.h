#ifndef TFITPARTICLEPINVTHETAPHIM_H
#define TFITPARTICLEPINVTHETAPHIM_H

#include "TMatrixD.h"
#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFitParticlePInvThetaPhiM: public TAbsFitParticle {
  
 public :

  TFitParticlePInvThetaPhiM();
  TFitParticlePInvThetaPhiM( const TFitParticlePInvThetaPhiM& fitParticle );
  TFitParticlePInvThetaPhiM(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePInvThetaPhiM(const TString &name, const TString &title, 
			    TLorentzVector* pini, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePInvThetaPhiM();
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
  
  ClassDef(TFitParticlePInvThetaPhiM, 1)    // Particle with spherical 4vector parametrization (1/p, theta, phi) and free energy
  
 private:
};

#endif
