#ifndef TFITPARTICLEPTETAPHIE_H
#define TFITPARTICLEPTETAPHIE_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticlePtEtaPhiE: public TAbsFitParticle {

public :

  TFitParticlePtEtaPhiE();
  TFitParticlePtEtaPhiE( const TFitParticlePtEtaPhiE& fitParticle );
  TFitParticlePtEtaPhiE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePtEtaPhiE(const TString &name, const TString &title, 
			TLorentzVector* pini, const TMatrixD* theCovMatrix);

  virtual ~TFitParticlePtEtaPhiE();
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

  ClassDef(TFitParticlePtEtaPhiE, 1)    // Particle with spherical 4vector parametrization and free mass

private:
  
};

#endif
