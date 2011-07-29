#ifndef TFITPARTICLERELPXPYPZE_H
#define TFITPARTICLERELPXPYPZE_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticleRelPxPyPzE: public TAbsFitParticle {

public :

  TFitParticleRelPxPyPzE();
  TFitParticleRelPxPyPzE( const TFitParticleRelPxPyPzE& fitParticle );
  TFitParticleRelPxPyPzE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticleRelPxPyPzE(const TString &name, const TString &title, 
			 TLorentzVector* pini, const TMatrixD* theCovMatrix);
  virtual ~TFitParticleRelPxPyPzE();
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

  ClassDef(TFitParticleRelPxPyPzE, 1)  // Particle with cartesian 4vector parametrization and free energy

private:
  
};

#endif
