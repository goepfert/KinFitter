
using namespace std;

#ifndef TFITCONSTRAINTPT_H
#define TFITCONSTRAINTPT_H

#include "KinFitter/TAbsFitConstraint.h"
#include "TAbsFitParticle.h"
#include "TMatrixD.h"
#include <vector>


class TFitConstraintPt: public TAbsFitConstraint {

public :

  TFitConstraintPt( );

  TFitConstraintPt(  vector<TAbsFitParticle*>* particles, 
		     Double_t constraint = 0.);

  TFitConstraintPt(  const TString &name, const TString &title,
		     vector<TAbsFitParticle*>* particles, 
		     Double_t constraint = 0.);
  virtual ~TFitConstraintPt();

  void addParticle( TAbsFitParticle* particle );
  void addParticles( TAbsFitParticle* p1, TAbsFitParticle* p2 = 0, TAbsFitParticle* p3 = 0, TAbsFitParticle* p4 = 0,
		     TAbsFitParticle* p5 = 0, TAbsFitParticle* p6 = 0, TAbsFitParticle* p7 = 0, TAbsFitParticle* p8 = 0,
		     TAbsFitParticle* p9 = 0, TAbsFitParticle* p10 = 0);
    
  // returns derivative df/dP with P=(p,E) and f the constraint f=0.
  // The matrix contains one row (df/dp, df/dE).
  virtual TMatrixD* getDerivative( TAbsFitParticle* particle );
  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();
  void setConstraint( Double_t constraint ) { _constraint = constraint; }

  virtual void print(); 
 
protected :


private:
  vector<TAbsFitParticle*> _particles;    // Vector containing constrained particles
  Double_t _constraint;                   // Value of constraint

  ClassDef(TFitConstraintPt, 1)   // Fit constraint: energy and momentum conservation
};

#endif
