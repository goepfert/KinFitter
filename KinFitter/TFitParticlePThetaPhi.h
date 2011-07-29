#ifndef TFITPARTICLEPTHETAPHI_H
#define TFITPARTICLEPTHETAPHI_H

#include "TMatrixD.h"
#include "KinFitter/TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFitParticlePThetaPhi: public TAbsFitParticle {

public :

  TFitParticlePThetaPhi();
  TFitParticlePThetaPhi( const TFitParticlePThetaPhi& fitParticle );
  TFitParticlePThetaPhi(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  TFitParticlePThetaPhi(const TString &name, const TString &title, 
			TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePThetaPhi();
  virtual TAbsFitParticle* clone( TString newname = "" ) const;

  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, ...) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).
  virtual TMatrixD* getDerivative();
  virtual TMatrixD* transform(const TLorentzVector& vec);
  virtual void setIni4Vec(const TLorentzVector* pini);
  void setIni4Vec(const TVector3* p, Double_t M);
  virtual TLorentzVector* calc4Vec( const TMatrixD* params );

protected :

  void init(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);

  ClassDef(TFitParticlePThetaPhi, 1)    // Particle with spherical 4vector parametrization and constrained mass
  
private:

};

#endif
