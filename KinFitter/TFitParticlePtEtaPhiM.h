#ifndef TFITPARTICLEPTETAPHIM_H
#define TFITPARTICLEPTETAPHIM_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticlePtEtaPhiM: public TAbsFitParticle {

public :

  TFitParticlePtEtaPhiM();
  TFitParticlePtEtaPhiM( const TFitParticlePtEtaPhiM& fitParticle );
  TFitParticlePtEtaPhiM(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePtEtaPhiM(const TString &name, const TString &title, 
			TLorentzVector* pini, const TMatrixD* theCovMatrix);
  
  virtual ~TFitParticlePtEtaPhiM();
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
  
  ClassDef(TFitParticlePtEtaPhiM, 1)    // Particle with spherical 4vector parametrization and free mass
  
private:
  
};

#endif
