#ifndef TFITPARTICLEPXPYPZE_H
#define TFITPARTICLEPXPYPZE_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticlePxPyPzE: public TAbsFitParticle {

public :

  TFitParticlePxPyPzE();
  TFitParticlePxPyPzE( const TFitParticlePxPyPzE& fitParticle );
  TFitParticlePxPyPzE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticlePxPyPzE(const TString &name, const TString &title, 
	              TLorentzVector* pini, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePxPyPzE();
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

  ClassDef(TFitParticlePxPyPzE, 1)  // Particle with cartesian 4vector parametrization and free energy

private:
  
};

#endif
