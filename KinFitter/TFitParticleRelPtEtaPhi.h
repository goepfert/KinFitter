#ifndef TFITPARTICLERELPTETAPHI_H
#define TFITPARTICLERELPTETAPHI_H


#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"


class TFitParticleRelPtEtaPhi: public TAbsFitParticle {

public :

  TFitParticleRelPtEtaPhi();
  TFitParticleRelPtEtaPhi( const TFitParticleRelPtEtaPhi& fitParticle );
  TFitParticleRelPtEtaPhi(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  TFitParticleRelPtEtaPhi(const TString &name, const TString &title, 
			 TVector3* p, Double_t M, const TMatrixD* theCovMatrix);

  virtual ~TFitParticleRelPtEtaPhi();
  virtual TAbsFitParticle* clone( TString newname = "" ) const;

  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, ...) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  void setIni4Vec(const TVector3* pini, Double_t M);


  virtual TLorentzVector* calc4Vec( const TMatrixD* params );

protected :

  void init(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);

  ClassDef(TFitParticleRelPtEtaPhi, 1)    // Particle with spherical 4vector parametrization and free mass

private:
  
};

#endif
