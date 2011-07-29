#ifndef TABSFITCONSTRAINT_H
#define TABSFITCONSTRAINT_H

#include <vector>
#include "TAbsFitParticle.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TString.h"

class TAbsFitConstraint : public TNamed {

public :

  TAbsFitConstraint();
  TAbsFitConstraint(const TString &name, const TString &title);
  virtual ~TAbsFitConstraint();

  // returns derivative df/dP with P=(p,E) and f the constraint f=0.
  // The matrix contains one row (df/dp, df/dE).
  virtual TMatrixD* getDerivative( TAbsFitParticle* particle ) = 0 ;
  virtual Double_t getInitValue() = 0;
  virtual Double_t getCurrentValue() = 0;

  // Methods needed for parameters in constraint
  Int_t getNPar() { return _nPar; } 
  virtual TMatrixD* getDerivativeAlpha() { return 0; }
  virtual const TMatrixD* getCovMatrix() const;
  virtual void setCovMatrix(const TMatrixD* theCovMatrix);
  virtual const TMatrixD* getCovMatrixFit() const;
  virtual void setCovMatrixFit(const TMatrixD* theCovMatrixFit);
  virtual const TMatrixD* getCovMatrixDeltaAlpha();
  const TMatrixD* getParIni();
  void  setParIni(const TMatrixD* parini);
  virtual void applyDeltaAlpha(TMatrixD* corrMatrix);
  const TMatrixD* getParCurr();

  virtual void print(); 
  virtual void reset();

protected :

  void calcCovMatrixDeltaAlpha();

  Int_t _nPar;                    // Number of free parameters 

  TMatrixD _covMatrix;            // covariance matrix
  TMatrixD _covMatrixFit;         // fitted covariance matrix
  TMatrixD _covMatrixDeltaAlpha;  // V(deltaAlpha) == V(alpha_meas) - V(alpha_fit)
  TMatrixD _iniparameters;        // initialized parameters (parameter values before the fit)
  TMatrixD _parameters;           // fitted parameters

  ClassDef(TAbsFitConstraint, 1)  // Abstract base class for fit constraints

};

#endif
