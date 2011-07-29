// //Classname: TFitConstraintEpGaus
// Author: Jan E. Sundermann, Verena Klose, Thomas Goepfert (TU Dresden)      


//________________________________________________________________
// 
// TFitConstraintEpGaus::
// --------------------
//
// Fit constraint: energy and momentum conservation
//                 with gaussian resolution
//


using namespace std;

#include "KinFitter/TFitConstraintEpGaus.h"
#include "KinFitter/TFitConstraintEp.h"

#include <iostream>
#include "TClass.h"

ClassImp(TFitConstraintEpGaus)

//----------------
// Constructor --
//----------------



TFitConstraintEpGaus::TFitConstraintEpGaus()
  :TFitConstraintEp()
{

  init();
  setConstraint(1., 1.);

}

TFitConstraintEpGaus::TFitConstraintEpGaus(vector<TAbsFitParticle*>* particles1, 
					   vector<TAbsFitParticle*>* particles2,
					   TFitConstraintEpGaus::component thecomponent, 
					   Double_t constraint,
					   Double_t width)
  :TFitConstraintEp( particles1, particles2, thecomponent, constraint )

{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEpGaus
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  init();
  setConstraint(constraint, width);

}


TFitConstraintEpGaus::TFitConstraintEpGaus(vector<TAbsFitParticle*>* particles,
					   TFitConstraintEpGaus::component thecomponent, 
					   Double_t constraint,
					   Double_t width)
  :TFitConstraintEp(particles, thecomponent, constraint)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEpGaus
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  init();
  setConstraint(constraint, width);

}


TFitConstraintEpGaus::TFitConstraintEpGaus(const TString &name, const TString &title,
					   vector<TAbsFitParticle*>* particles1, 
					   vector<TAbsFitParticle*>* particles2,
					   TFitConstraintEpGaus::component thecomponent, 
					   Double_t constraint,
					   Double_t width )
  :TFitConstraintEp(name, title, particles1, particles2, thecomponent, constraint)

{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEpGaus
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  init();
  setConstraint(constraint, width);

}

TFitConstraintEpGaus::TFitConstraintEpGaus(const TString &name, const TString &title,
					   vector<TAbsFitParticle*>* particles,
					   TFitConstraintEpGaus::component thecomponent, 
					   Double_t constraint,
					   Double_t width )
  :TFitConstraintEp(name, title, particles, thecomponent, constraint)
{
  // particles: vector containing pointer to TAbsFitParticle objects. 
  //            Energy or momentum conservation will be calculated for
  //            those particles.
  // thecomponent: conserved 4vector component ( pX, pY, pZ, E ). For
  //            full 4vector conservation four objects of type TFitConstraintEpGaus
  //            are needed (four constraints)
  // constraint:value of energy or momentum constraint ( e.g. sum[ pX_i ] = constraint )

  init();
  setConstraint(constraint, width);

}


void
TFitConstraintEpGaus::init() {

  _nPar = 1;
  _iniparameters.ResizeTo(1,1);
  _iniparameters(0,0) = 1.;
  _parameters.ResizeTo(1,1);
  _parameters = _iniparameters;

}



//--------------
// Destructor --
//--------------
TFitConstraintEpGaus::~TFitConstraintEpGaus() {

}



//--------------
// Operations --
//--------------
void TFitConstraintEpGaus::setConstraint( Double_t constraint, Double_t width ) {

  _constraint = constraint;
  _width = width;
  setCovMatrix( 0 );
  _covMatrix(0,0) = (_width*_width) / (_constraint*_constraint);

}


Double_t TFitConstraintEpGaus::getInitValue() {
  // Get initial value of constraint (before the fit)


  Double_t InitValue = 0.; 
  UInt_t Npart1 = _particles1.size();
  UInt_t Npart2 = _particles2.size();

  for (UInt_t i=0;i<Npart1;i++) {
    const TLorentzVector* FourVec = _particles1[i]->getIni4Vec();
    InitValue += (*FourVec)[(int) _component];
  }  

  for (UInt_t i=0;i<Npart2;i++) {
    const TLorentzVector* FourVec = _particles2[i]->getIni4Vec();
    InitValue -= (*FourVec)[(int) _component];
  }
  
  InitValue -= _iniparameters(0,0) * _constraint;
  return InitValue;

}


Double_t TFitConstraintEpGaus::getCurrentValue() {
  // Get value of constraint after the fit

  Double_t CurrentValue = 0.;

  UInt_t Npart1 = _particles1.size();
  for (UInt_t i=0;i<Npart1;i++) {
    const TLorentzVector* FourVec = _particles1[i]->getCurr4Vec();
    CurrentValue += (*FourVec)[(int) _component];
  } 

  UInt_t Npart2 = _particles2.size();
  for (UInt_t i=0;i<Npart2;i++) {
    const TLorentzVector* FourVec = _particles2[i]->getCurr4Vec();
    CurrentValue -= (*FourVec)[(int) _component];
  } 

  CurrentValue -= _parameters(0,0) * _constraint;
  return CurrentValue;

}



TMatrixD* TFitConstraintEpGaus::getDerivativeAlpha() { 
  // Calculate dF/dAlpha = -1 * C

  TMatrixD* DerivativeMatrix = new TMatrixD(1,1);
  DerivativeMatrix->Zero();

  (*DerivativeMatrix)(0,0) = -1. * _constraint;

  return DerivativeMatrix;

}


void TFitConstraintEpGaus::print() {

  cout << "__________________________" << endl << endl;
  cout <<"OBJ: " << IsA()->GetName() << "\t" << GetName() << "\t" << GetTitle() << endl;

  cout << "initial value: " << getInitValue() << endl;
  cout << "current value: " << getCurrentValue() << endl;
  cout << "component: " << _component << endl;
  cout << "constraint: " << _constraint << endl;
  cout << "width: " << _width << endl;

}
