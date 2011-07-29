
using namespace std;

#ifndef TFITCONSTRAINTEPGAUS_H
#define TFITCONSTRAINTEPGAUS_H

#include "TFitConstraintEp.h"
#include "TAbsFitParticle.h"
#include "TMatrixD.h"
#include <vector>


class TFitConstraintEpGaus: public TFitConstraintEp {

public :

  TFitConstraintEpGaus( );

  TFitConstraintEpGaus( vector<TAbsFitParticle*>* particles1, 
			vector<TAbsFitParticle*>* particles2,
			TFitConstraintEpGaus::component thecomponent, 
			Double_t constraint = 0.,
			Double_t width = 0.);

  TFitConstraintEpGaus(  vector<TAbsFitParticle*>* particles, 
			 TFitConstraintEpGaus::component thecomponent, 
			 Double_t constraint = 0.,
			 Double_t width = 0.);

  TFitConstraintEpGaus(  const TString &name, const TString &title,
			 vector<TAbsFitParticle*>* particles1,
			 vector<TAbsFitParticle*>* particles2, 
			 TFitConstraintEpGaus::component thecomponent, 
			 Double_t constraint = 0.,
			 Double_t width = 0.);

  TFitConstraintEpGaus(  const TString &name, const TString &title,
			 vector<TAbsFitParticle*>* particles,
			 TFitConstraintEpGaus::component thecomponent, 
			 Double_t constraint = 0.,
			 Double_t width = 0.);

  virtual ~TFitConstraintEpGaus();

  
  void init();


  // returns derivative df/dP with P=(p,E) and f the constraint f=0.
  // The matrix contains one row (df/dp, df/dE).
  virtual TMatrixD* getDerivativeAlpha();
  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();



  void setConstraint( Double_t constraint, Double_t width );

  virtual void print(); 
 
protected :

  Double_t _width;                   // Momentum/Energy reolution

private:

  ClassDef(TFitConstraintEpGaus, 1)   // Fit constraint: energy and momentum conservation

};

#endif
