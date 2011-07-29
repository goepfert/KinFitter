using namespace std;

#ifndef TFITCONSTRAINTM_H
#define TFITCONSTRAINTM_H


#include "KinFitter/TAbsFitConstraint.h"
#include <vector>

#if ROOT_VERSION_CODE >= 329728
#include "TMatrixDfwd.h"
#else
#include "TMatrixD.h"
#endif

class TAbsFitParticle;

class TFitConstraintM: public TAbsFitConstraint {

public :

  TFitConstraintM();
  TFitConstraintM(vector<TAbsFitParticle*>* ParList1,
		  vector<TAbsFitParticle*>* ParList2,
		  Double_t Mass = 0);
  TFitConstraintM(const TString &name, const TString &title,
		  vector<TAbsFitParticle*>* ParList1,
		  vector<TAbsFitParticle*>* ParList2,
		  Double_t Mass = 0);

  virtual ~TFitConstraintM();

  void addParticle1( TAbsFitParticle* particle );
  void addParticle2( TAbsFitParticle* particle );
  void addParticles1( TAbsFitParticle* p1, TAbsFitParticle* p2 = 0, TAbsFitParticle* p3 = 0, TAbsFitParticle* p4 = 0,
		      TAbsFitParticle* p5 = 0, TAbsFitParticle* p6 = 0, TAbsFitParticle* p7 = 0, TAbsFitParticle* p8 = 0,
		      TAbsFitParticle* p9 = 0, TAbsFitParticle* p10 = 0);
  void addParticles2( TAbsFitParticle* p1, TAbsFitParticle* p2 = 0, TAbsFitParticle* p3 = 0, TAbsFitParticle* p4 = 0,
		      TAbsFitParticle* p5 = 0, TAbsFitParticle* p6 = 0, TAbsFitParticle* p7 = 0, TAbsFitParticle* p8 = 0,
		      TAbsFitParticle* p9 = 0, TAbsFitParticle* p10 = 0);
  void setMassConstraint(Double_t Mass) { _TheMassConstraint = Mass; }

  // returns derivative df/dP with P=(p,E) and f the constraint f=0 for
  // one particle. The matrix contains one row (df/dp, df/dE).
  virtual TMatrixD* getDerivative( TAbsFitParticle* particle );
  virtual Double_t getInitValue();
  virtual Double_t getCurrentValue();

  Bool_t OnList(vector<TAbsFitParticle*>* List, TAbsFitParticle* particle);
  Double_t CalcMass(vector<TAbsFitParticle*>* List, Bool_t IniVal);

  virtual void print();

protected :

  vector<TAbsFitParticle*> _ParList1;   // Vector containing first list of constrained particles ( sum[ m_i ] - sum[ m_j ] == 0 )
  vector<TAbsFitParticle*> _ParList2;   // Vector containing second list of constrained particles ( sum[ m_i ] - sum[ m_j ] == 0 )
  Double_t _TheMassConstraint;

  ClassDef(TFitConstraintM, 1)   // Fit constraint: mass conservation
};

#endif
