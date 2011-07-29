// //Classname: TFitConstraintEp
// Author: Jan E. Sundermann, Verena Klose (TU Dresden)      


//________________________________________________________________
// 
// TFitConstraintEp::
// --------------------
//
// Fit constraint: energy and momentum conservation
//
//


using namespace std;

#include "KinFitter/TFitConstraintEp.h"

#include <iostream>
#include "TClass.h"

ClassImp(TFitConstraintEp)

//----------------
// Constructor --
//----------------



TFitConstraintEp::TFitConstraintEp()
  :TAbsFitConstraint()
  ,_particles1(0)
  ,_particles2(0)
  ,_constraint(0.)
  ,_component(TFitConstraintEp::pX)
{}

TFitConstraintEp::TFitConstraintEp(vector<TAbsFitParticle*>* particles1, 
				   vector<TAbsFitParticle*>* particles2,
				   TFitConstraintEp::component thecomponent, 
				   Double_t constraint)
  :TAbsFitConstraint()
  ,_particles1(0)
  ,_particles2(0)
  ,_constraint(constraint)
  ,_component(thecomponent)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEp
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  if (particles1) {
    _particles1 = (*particles1);
  }
  if (particles2) {
    _particles2 = (*particles2);
  }
}


TFitConstraintEp::TFitConstraintEp(vector<TAbsFitParticle*>* particles,
				   TFitConstraintEp::component thecomponent, 
				   Double_t constraint)
  :TAbsFitConstraint()
  ,_particles1(0)
  ,_particles2(0)
  ,_constraint(constraint)
  ,_component(thecomponent)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEp
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  if (particles) {
    _particles1 = (*particles);
  }
}


TFitConstraintEp::TFitConstraintEp(const TString &name, const TString &title,
				   vector<TAbsFitParticle*>* particles1, 
				   vector<TAbsFitParticle*>* particles2,
				   TFitConstraintEp::component thecomponent, 
				   Double_t constraint)
  :TAbsFitConstraint(name, title)
  ,_particles1(0)
  ,_particles2(0)
  ,_constraint(constraint)
  ,_component(thecomponent)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEp
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  if (particles1) {
    _particles1 = (*particles1);
  }
  if (particles2) {
    _particles2 = (*particles2);
  }
}

TFitConstraintEp::TFitConstraintEp(const TString &name, const TString &title,
				   vector<TAbsFitParticle*>* particles,
				   TFitConstraintEp::component thecomponent, 
				   Double_t constraint)
  :TAbsFitConstraint(name, title)
  ,_particles1(0)
  ,_particles2(0)
  ,_constraint(constraint)
  ,_component(thecomponent)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEp
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  if (particles) {
    _particles1 = (*particles);
  }
}


//--------------
// Destructor --
//--------------
TFitConstraintEp::~TFitConstraintEp() {

}

void TFitConstraintEp::addParticle1( TAbsFitParticle* particle ) {
  // Add one particles to list of constrained particles

  _particles1.push_back( particle );

}

void TFitConstraintEp::addParticle2( TAbsFitParticle* particle ) {
  // Add one particles to list of constrained particles

  _particles2.push_back( particle );

}

void TFitConstraintEp::addParticle( TAbsFitParticle* particle ) {
  // Add one particles to list of constrained particles

  _particles1.push_back( particle );

}



void TFitConstraintEp::addParticles1( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, TAbsFitParticle* p4,
				     TAbsFitParticle* p5, TAbsFitParticle* p6, TAbsFitParticle* p7, TAbsFitParticle* p8,
				     TAbsFitParticle* p9, TAbsFitParticle* p10) {
  // Add many particles to list of constrained particles

  if (p1) addParticle1( p1 );
  if (p2) addParticle1( p2 );
  if (p3) addParticle1( p3 );
  if (p4) addParticle1( p4 );
  if (p5) addParticle1( p5 );
  if (p6) addParticle1( p6 );
  if (p7) addParticle1( p7 );
  if (p8) addParticle1( p8 );
  if (p9) addParticle1( p9 );
  if (p10) addParticle1( p10 );

}

void TFitConstraintEp::addParticles2( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, TAbsFitParticle* p4,
				     TAbsFitParticle* p5, TAbsFitParticle* p6, TAbsFitParticle* p7, TAbsFitParticle* p8,
				     TAbsFitParticle* p9, TAbsFitParticle* p10) {
  // Add many particles to list of constrained particles

  if (p1) addParticle2( p1 );
  if (p2) addParticle2( p2 );
  if (p3) addParticle2( p3 );
  if (p4) addParticle2( p4 );
  if (p5) addParticle2( p5 );
  if (p6) addParticle2( p6 );
  if (p7) addParticle2( p7 );
  if (p8) addParticle2( p8 );
  if (p9) addParticle2( p9 );
  if (p10) addParticle2( p10 );

}


void TFitConstraintEp::addParticles( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, TAbsFitParticle* p4,
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
TMatrixD* TFitConstraintEp::getDerivative( TAbsFitParticle* particle ) {
  // returns derivative df/dP with P=(p,E) and f the constraint (f=0).
  // The matrix contains one row (df/dp, df/dE).

  TMatrixD* DerivativeMatrix = new TMatrixD(1,4);
  (*DerivativeMatrix) *= 0.;
  if( OnList( &_particles1, particle) ) {
    (*DerivativeMatrix)(0,(int) _component) = 1.;
  }
  if( OnList( &_particles2, particle) ) {
    (*DerivativeMatrix)(0,(int) _component) = -1.;
  }

  return DerivativeMatrix;
  
}

Bool_t TFitConstraintEp::OnList(vector<TAbsFitParticle*>* List,
			       TAbsFitParticle* particle) {
  // Checks whether list contains given particle

  Bool_t ok(false);
  UInt_t Npart = List->size();
  for (UInt_t i=0;i<Npart;i++) {
    ok = (particle == (*List)[i]);
    if (ok) break;
  }
  return ok;
}


Double_t TFitConstraintEp::getInitValue() {
  // Get initial value of constraint (before the fit)

  Double_t InitValue(0) ; 
  UInt_t Npart1 = _particles1.size();
  UInt_t Npart2 = _particles2.size();

  for (UInt_t i = 0; i < Npart1; i++) {
    const TLorentzVector* FourVec = _particles1[i]->getIni4Vec();
    InitValue += (*FourVec)[(int) _component];
  }  

  for (UInt_t i = 0; i < Npart2; i++) {
    const TLorentzVector* FourVec = _particles2[i]->getIni4Vec();
    InitValue -= (*FourVec)[(int) _component];
  }
  
  InitValue -= _constraint;
  return InitValue;

}

Double_t TFitConstraintEp::getCurrentValue() {
  // Get value of constraint after the fit

  Double_t CurrentValue(0);

  UInt_t Npart1 = _particles1.size();
  for (UInt_t i = 0; i < Npart1; i++) {
    const TLorentzVector* FourVec = _particles1[i]->getCurr4Vec();
    CurrentValue += (*FourVec)[(int) _component];
  } 

  UInt_t Npart2 = _particles2.size();
  for (UInt_t i = 0; i < Npart2; i++) {
    const TLorentzVector* FourVec = _particles2[i]->getCurr4Vec();
    CurrentValue -= (*FourVec)[(int) _component];
  } 

  CurrentValue -= _constraint;
  return CurrentValue;

}

void TFitConstraintEp::print() {

  cout << "__________________________" << endl << endl;
  cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

  cout << "initial value: " << getInitValue() << endl;
  cout << "current value: " << getCurrentValue() << endl;
  cout << "component: " << _component << endl;
  cout << "constraint: " << _constraint << endl;

}
