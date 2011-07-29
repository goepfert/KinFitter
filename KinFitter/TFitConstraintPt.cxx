// //Classname: TFitConstraintPt
// Author: Jan E. Sundermann, Verena Klose, Thomas Goepfert (TU Dresden)      

//________________________________________________________________
// 
// TFitConstraintPt::
// --------------------
//
// Fit constraint: Pt conservation ( pt_i - constraint == 0)
//                 Pt fixed to certain value
//


using namespace std;

#include "KinFitter/TFitConstraintPt.h"

#include <iostream>

#include "TClass.h"

ClassImp(TFitConstraintPt)

//----------------
// Constructor --
//----------------
TFitConstraintPt::TFitConstraintPt()
  :TAbsFitConstraint()
  ,_constraint(0.)
  ,_particles(0)
{}

TFitConstraintPt::TFitConstraintPt(vector<TAbsFitParticle*>* particles, 
				   Double_t constraint)
  :TAbsFitConstraint()
  ,_constraint(constraint)
  ,_particles(0)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Pt constraint will be calculated for those particles.

  if (particles) {
    _particles = (*particles);
  }
}

TFitConstraintPt::TFitConstraintPt(const TString &name, const TString &title,
				   vector<TAbsFitParticle*>* particles, 
				   Double_t constraint)
  :TAbsFitConstraint(name, title)
  ,_constraint(constraint)
  ,_particles(0)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Pt constraint will be calculated for those particles.

  if (particles) {
    _particles = (*particles);
  }
}

//--------------
// Destructor --
//--------------
TFitConstraintPt::~TFitConstraintPt() {

}

void TFitConstraintPt::addParticle( TAbsFitParticle* particle ) {
  // Add one particles to list of constrained particles

  _particles.push_back( particle );

}

void TFitConstraintPt::addParticles( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, TAbsFitParticle* p4,
				     TAbsFitParticle* p5, TAbsFitParticle* p6, TAbsFitParticle* p7, TAbsFitParticle* p8,
				     TAbsFitParticle* p9, TAbsFitParticle* p10) {
  // Add many particles to list of constrained particles

  if (p1) addParticle( p1 );
  if (p2) addParticle( p2 );
  if (p3) addParticle( p3 );
  if (p4) addParticle( p4 );
  if (p5) addParticle( p5 );
  if (p6) addParticle( p6 );
  if (p7) addParticle( p7 );
  if (p8) addParticle( p8 );
  if (p9) addParticle( p9 );
  if (p10) addParticle( p10 );

}

//--------------
// Operations --
//--------------
TMatrixD* TFitConstraintPt::getDerivative( TAbsFitParticle* particle ) {
  // returns derivative df/dP with P=(p,E) and f the constraint (f=0).
  // The matrix contains one row (df/dp, df/dE).

  UInt_t Npart = _particles.size();
  TLorentzVector CurrentVec(0., 0., 0., 0.);
  for (int i=0; i < Npart; i++) {
    CurrentVec += *(_particles[i]->getCurr4Vec());
  }

  Double_t Pt = CurrentVec.Pt();
  Double_t Px = CurrentVec.Px();
  Double_t Py = CurrentVec.Py();

  TMatrixD* DerivativeMatrix = new TMatrixD(1,4);
  (*DerivativeMatrix)(0,0) = Px / Pt;
  (*DerivativeMatrix)(0,1) = Py / Pt;
  (*DerivativeMatrix)(0,2) = 0.;
  (*DerivativeMatrix)(0,3) = 0.;
  return DerivativeMatrix;

}


Double_t TFitConstraintPt::getInitValue() {
  // Get initial value of constraint (before the fit)

  TLorentzVector InitVec(0., 0., 0., 0.);
  UInt_t Npart = _particles.size();
  for (int i=0;i<Npart;i++) {
    InitVec += *(_particles[i]->getIni4Vec());
  }
  Double_t InitValue = InitVec.Pt(); 
  InitValue -= _constraint;
  return InitValue;
}

Double_t TFitConstraintPt::getCurrentValue() {
  // Get value of constraint after the fit

  UInt_t Npart = _particles.size();
  TLorentzVector CurrentVec(0., 0., 0., 0.);
  for (int i=0; i < Npart; i++) {
    CurrentVec += *(_particles[i]->getCurr4Vec());
  }
  
  Double_t CurrentValue = CurrentVec.Pt();
  CurrentValue -= _constraint;
  return CurrentValue;
}

void TFitConstraintPt::print() {

  cout << "__________________________" << endl << endl;
  cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

  cout << "initial value: " << getInitValue() << endl;
  cout << "current value: " << getCurrentValue() << endl;
  cout << "constraint: " << _constraint << endl;

}
