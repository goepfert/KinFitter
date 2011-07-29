#ifndef TFITPARTICLERELE_H
#define TFITPARTICLERELE_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticleRelE: public TAbsFitParticle {

public :

  TFitParticleRelE();
  TFitParticleRelE( const TFitParticleRelE& fitParticle );
  TFitParticleRelE(TLorentzVector* pini, const TMatrixD* theCovMatrix);
  TFitParticleRelE(const TString &name, const TString &title, 
			 TLorentzVector* pini,
			 const TMatrixD* theCovMatrix);
  virtual ~TFitParticleRelE();
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

  ClassDef(TFitParticleRelE, 1)

private:
  
};

#endif
