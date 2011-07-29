using namespace std;

#ifndef TFITCONSTRAINTMGAUS_H
#define TFITCONSTRAINTMGAUS_H

#include "TFitConstraintM.h"

#include <vector>

#if ROOT_VERSION_CODE >= 329728
#include "TMatrixDfwd.h"
#else
#include "TMatrixD.h"
#endif

class TAbsFitParticle;

class TFitConstraintMGaus: public TFitConstraintM {

public :

  TFitConstraintMGaus();
  TFitConstraintMGaus(vector<TAbsFitParticle*>* ParList1,
		      vector<TAbsFitParticle*>* ParList2,
		      Double_t Mass = 0, Double_t Width = 0);
  TFitConstraintMGaus(const TString &name, const TString &title,
		      vector<TAbsFitParticle*>* ParList1,
		      vector<TAbsFitParticle*>* ParList2,
		      Double_t Mass = 0, Double_t Width = 0);

  virtual ~TFitConstraintMGaus();

  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();
  virtual TMatrixD* getDerivativeAlpha();

  void setMassConstraint(Double_t Mass, Double_t Width);

  virtual void print();

protected :

  Double_t _width;

  void init();

  ClassDef(TFitConstraintMGaus, 1)   // Fit constraint: mass conservation

};

#endif

