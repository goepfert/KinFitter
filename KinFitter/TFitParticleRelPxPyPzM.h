#ifndef TFITPARTICLERELPXPYPZM_H
#define TFITPARTICLERELPXPYPZM_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticleRelPxPyPzM: public TAbsFitParticle {

public :

  TFitParticleRelPxPyPzM();
  TFitParticleRelPxPyPzM( const TFitParticleRelPxPyPzM& fitParticle );
  TFitParticleRelPxPyPzM(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticleRelPxPyPzM(const TString &name, const TString &title, 
			 TLorentzVector* pini, const TMatrixD* theCovMatrix);
  virtual ~TFitParticleRelPxPyPzM();
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

  ClassDef(TFitParticleRelPxPyPzM, 1)  // Particle with cartesian 4vector parametrization

private:
  
};

#endif
