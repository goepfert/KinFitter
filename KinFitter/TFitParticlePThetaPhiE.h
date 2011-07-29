#ifndef TFITPARTICLEPTHETAPHIE_H
#define TFITPARTICLEPTHETAPHIE_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticlePThetaPhiE: public TAbsFitParticle {

public :

  TFitParticlePThetaPhiE();
  TFitParticlePThetaPhiE( const TFitParticlePThetaPhiE& fitParticle );
  TFitParticlePThetaPhiE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePThetaPhiE(const TString &name, const TString &title, 
			 TLorentzVector* pini,
			 const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePThetaPhiE();
  virtual TAbsFitParticle* clone( TString newname = "" ) const;

  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, ...) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  virtual TLorentzVector* calc4Vec( const TMatrixD* params );

protected :

  void init(TLorentzVector* pini, const TMatrixD* theCovMatrix);

  ClassDef(TFitParticlePThetaPhiE, 1)    // Particle with spherical 4vector parametrization (p, theta, phi) and free energy

private:
  
};

#endif
