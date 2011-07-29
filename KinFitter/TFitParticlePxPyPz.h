#ifndef TFITPARTICLEPXPYPZ_H
#define TFITPARTICLEPXPYPZ_H

#include "TMatrixD.h"
#include "TAbsFitParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFitParticlePxPyPz: public TAbsFitParticle {

public :

  TFitParticlePxPyPz();
  TFitParticlePxPyPz( const TFitParticlePxPyPz& fitParticle );
  TFitParticlePxPyPz(TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  TFitParticlePxPyPz(const TString &name, const TString &title, 
		     TVector3* p, Double_t M, const TMatrixD* theCovMatrix);
  virtual ~TFitParticlePxPyPz();
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

  ClassDef(TFitParticlePxPyPz, 1)    // Particle with cartesian 4vector parametrization and constrained mass
  
private:
};

#endif
