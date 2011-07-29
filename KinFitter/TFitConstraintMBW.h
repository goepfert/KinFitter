using namespace std;

#ifndef TFITCONSTRAINTMBW_H
#define TFITCONSTRAINTMBW_H

#include "TFitConstraintM.h"

#include <vector>

#if ROOT_VERSION_CODE >= 329728
#include "TMatrixDfwd.h"
#else
#include "TMatrixD.h"
#endif

class TAbsFitParticle;

class TFitConstraintMBW: public TFitConstraintM {

public :

  TFitConstraintMBW();
  TFitConstraintMBW(vector<TAbsFitParticle*>* ParList1,
		    vector<TAbsFitParticle*>* ParList2,
		    Double_t Mass = 0, Double_t Width = 0);
  TFitConstraintMBW(const TString &name, const TString &title,
		    vector<TAbsFitParticle*>* ParList1,
		    vector<TAbsFitParticle*>* ParList2,
		    Double_t Mass = 0, Double_t Width = 0);
  
  virtual ~TFitConstraintMBW();
  
  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();
  virtual TMatrixD* getDerivativeAlpha();

  void setMassConstraint(Double_t Mass, Double_t Width);
  
  virtual void print();
  Double_t transformBWToGaus( Double_t M );
  Double_t transformGausToBW( Double_t M );
  
  
protected :
  
  Double_t _width;
  Double_t _sigma;
  
  void init();
  
  ClassDef(TFitConstraintMBW, 1)   // Fit constraint: mass conservation
  
};

#endif

