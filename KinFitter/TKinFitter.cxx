// //Classname: TKinFitter
// Author: Jan E. Sundermann (U Freiburg), Verena Klose, Thomas Goepfert (TU Dresden)


//________________________________________________________________
// 
// TKinFitter::
// --------------------
//
// //Class to perm kinematic fit with non-linear constraints
//


using namespace std;

#include <iostream>
#include <iomanip>

#include "KinFitter/TKinFitter.h"
#include "KinFitter/TAbsFitParticle.h"
#include "KinFitter/TAbsFitConstraint.h"

#include "TArrayI.h"

ClassImp(TKinFitter)

TKinFitter::TKinFitter():
  TNamed("UnNamed", "UnTitled"),
  _A(1, 1),
  _AT(1, 1),
  _B(1, 1),
  _BT(1, 1),
  _V(1, 1),
  _Vinv(1, 1),
  _VB(1, 1),
  _VBinv(1, 1),
  _VA(1, 1),
  _VAinv(1, 1),
  _c(1, 1),
  _C11(1, 1),
  _C11T(1, 1),
  _C21(1, 1),
  _C21T(1, 1),
  _C22(1, 1),
  _C22T(1, 1),
  _C31(1, 1),
  _C31T(1, 1),
  _C32(1, 1),
  _C32T(1, 1),
  _C33(1, 1),
  _C33T(1, 1),
  _deltaA(1, 1),
  _deltaY(1, 1),
  _lambda(1, 1),
  _lambdaT(1, 1),
  _lambdaVFit(1, 1),
  _yaVFit(1, 1),
  _constraints(0),
  _particles(0),
  _paramMeasured(0)
{

  reset();

}

TKinFitter::TKinFitter(const TString &name, const TString &title):
  TNamed(name, title),
  _A(1, 1),
  _AT(1, 1),
  _B(1, 1),
  _BT(1, 1),
  _V(1, 1),
  _Vinv(1, 1),
  _VB(1, 1),
  _VBinv(1, 1),
  _VA(1, 1),
  _VAinv(1, 1),
  _c(1, 1),
  _C11(1, 1),
  _C11T(1, 1),
  _C21(1, 1),
  _C21T(1, 1),
  _C22(1, 1),
  _C22T(1, 1),
  _C31(1, 1),
  _C31T(1, 1),
  _C32(1, 1),
  _C32T(1, 1),
  _C33(1, 1),
  _C33T(1, 1),
  _deltaA(1, 1),
  _deltaY(1, 1),
  _lambda(1, 1),
  _lambdaT(1, 1),
  _lambdaVFit(1, 1),
  _yaVFit(1, 1),
  _constraints(0),
  _particles(0),
  _paramMeasured(0)
{

  reset();

}

void TKinFitter::reset() {
  // reset all internal parameters of the fitter

  _status = -1;
  _nbIter = 0;
  _nParA = 0;
  _nParB = 0;
  _verbosity = 1;
  _A.ResizeTo(1, 1);
  _AT.ResizeTo(1, 1);
  _B.ResizeTo(1, 1);
  _BT.ResizeTo(1, 1);
  _V.ResizeTo(1, 1);
  _Vinv.ResizeTo(1, 1);
  _VB.ResizeTo(1, 1);
  _VBinv.ResizeTo(1, 1);
  _VA.ResizeTo(1, 1);
  _VAinv.ResizeTo(1, 1);
  _c.ResizeTo(1, 1);
  _C11.ResizeTo(1, 1);
  _C11T.ResizeTo(1, 1);
  _C21.ResizeTo(1, 1);
  _C21T.ResizeTo(1, 1);
  _C22.ResizeTo(1, 1);
  _C22T.ResizeTo(1, 1);
  _C31.ResizeTo(1, 1);
  _C31T.ResizeTo(1, 1);
  _C32.ResizeTo(1, 1);
  _C32T.ResizeTo(1, 1);
  _C33.ResizeTo(1, 1);
  _C33T.ResizeTo(1, 1);
  _deltaA.ResizeTo(1, 1);
  _deltaY.ResizeTo(1, 1);
  _lambda.ResizeTo(1, 1);
  _lambdaT.ResizeTo(1, 1);
  _lambdaVFit.ResizeTo(1, 1);
  _yaVFit.ResizeTo(1, 1);

  _constraints.clear();
  _particles.clear();
  _paramMeasured.clear();

  // Set to default values
  _maxNbIter = 50;
  _maxDeltaS = 5e-3;
  _maxF =  1e-4;

}

void TKinFitter::resetStatus() {
  // reset status of the fit

  _status = -1;
  _nbIter = 0;

}

void TKinFitter::resetParams() {
  // reset all particles and contraints to their 
  // initial parameter values

  for (UInt_t iP = 0; iP < _particles.size(); iP++) {
    TAbsFitParticle* particle = _particles[iP];
    particle->reset();
  }
  for (UInt_t iC = 0; iC < _constraints.size(); iC++) {
    TAbsFitConstraint* constraint = _constraints[iC];
    constraint->reset();
  }


}

TKinFitter::~TKinFitter() {

}

void TKinFitter::countMeasParams() {
  // count number of measured parameters

  _nParB = 0;
  for (UInt_t indexParameter = 0; indexParameter < _paramMeasured.size(); indexParameter++) {
    if ( _paramMeasured[indexParameter] == true ) {
      _nParB++;
    }
  }
  for (UInt_t indexConstraint = 0; indexConstraint < _constraints.size(); indexConstraint++) {
    _nParB += _constraints[indexConstraint]->getNPar();
  }

}

void TKinFitter::countUnmeasParams() {
  // count number of unmeasured parameters

  _nParA = 0;
  for (UInt_t indexParameter = 0; indexParameter < _paramMeasured.size(); indexParameter++) {
    if ( _paramMeasured[indexParameter] == false ) {
      _nParA++;
    }
  }

}

void TKinFitter::addMeasParticle( TAbsFitParticle* particle ) {
  // add measured particle

  resetStatus();

  if ( particle != 0 ) {
    _particles.push_back( particle );
    for ( Int_t par = 1; par <= particle->getNPar(); par++ ) {
      _paramMeasured.push_back( true );
    }
  } else {
    cout << "Particle points to NULL." << endl;
  }

  countMeasParams();

}

void TKinFitter::addMeasParticles( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, 
				   TAbsFitParticle* p4, TAbsFitParticle* p5, TAbsFitParticle* p6,
				   TAbsFitParticle* p7, TAbsFitParticle* p8, TAbsFitParticle* p9) {
  // add many measured particles

  resetStatus();
  Int_t measParams = 0;

  // Add parameters and count added measured parameters
  if ( p1 != 0 ) {
    _particles.push_back( p1 );
    measParams += p1->getNPar();
  }
  if ( p2 != 0 ) {
    _particles.push_back( p2 );
    measParams += p2->getNPar();
  }
  if ( p3 != 0 ) {
    _particles.push_back( p3 );
    measParams += p3->getNPar();
  }
  if ( p4 != 0 ) {
    _particles.push_back( p4 );
    measParams += p4->getNPar();
  }
  if ( p5 != 0 ) {
    _particles.push_back( p5 );
    measParams += p5->getNPar();
  }
  if ( p6 != 0 ) {
    _particles.push_back( p6 );
    measParams += p6->getNPar();
  }
  if ( p7 != 0 ) {
    _particles.push_back( p7 );
    measParams += p7->getNPar();
  }
  if ( p8 != 0 ) {
    _particles.push_back( p8 );
    measParams += p8->getNPar();
  }
  if ( p9 != 0 ) {
    _particles.push_back( p9 );
    measParams += p9->getNPar();
  }

  // add measured parameters to list of parameters
  for ( Int_t par = 1; par <= measParams; par++ ) {
    _paramMeasured.push_back( true );
  }

  countMeasParams();

}

void TKinFitter::addUnmeasParticle( TAbsFitParticle* particle ) {
  // add one unmeasured particle

  resetStatus();

  if ( particle != 0 ) {
    _particles.push_back( particle );
    for ( Int_t par = 1; par <= particle->getNPar(); par++ ) {
      _paramMeasured.push_back( false );
    }
  } else {
    cout << "Particle points to NULL." << endl;
  }

  countUnmeasParams();

}

void TKinFitter::addUnmeasParticles( TAbsFitParticle* p1, TAbsFitParticle* p2, TAbsFitParticle* p3, 
				     TAbsFitParticle* p4, TAbsFitParticle* p5, TAbsFitParticle* p6,
				     TAbsFitParticle* p7, TAbsFitParticle* p8, TAbsFitParticle* p9) {
  // add many unmeasured particles

  resetStatus();
  Int_t unmeasParams = 0;

  // Add parameters and count added unmeasured parameters
  if ( p1 != 0 ) {
    _particles.push_back( p1 );
    unmeasParams += p1->getNPar();
  }
  if ( p2 != 0 ) {
    _particles.push_back( p2 );
    unmeasParams += p2->getNPar();
  }
  if ( p3 != 0 ) {
    _particles.push_back( p3 );
    unmeasParams += p3->getNPar();
  }
  if ( p4 != 0 ) {
    _particles.push_back( p4 );
    unmeasParams += p4->getNPar();
  }
  if ( p5 != 0 ) {
    _particles.push_back( p5 );
    unmeasParams += p5->getNPar();
  }
  if ( p6 != 0 ) {
    _particles.push_back( p6 );
    unmeasParams += p6->getNPar();
  }
  if ( p7 != 0 ) {
    _particles.push_back( p7 );
    unmeasParams += p7->getNPar();
  }
  if ( p8 != 0 ) {
    _particles.push_back( p8 );
    unmeasParams += p8->getNPar();
  }
  if ( p9 != 0 ) {
    _particles.push_back( p9 );
    unmeasParams += p9->getNPar();
  }

  // add unmeasured parameters to list of parameters
  for ( Int_t par = 1; par <= unmeasParams; par++ ) {
    _paramMeasured.push_back( false );
  }

  countUnmeasParams();

}

void TKinFitter::addConstraint( TAbsFitConstraint* constraint ) {
  // add constraint

  resetStatus();

  if ( constraint != 0 ) {
    _constraints.push_back( constraint );
  }

  countMeasParams();

}

void TKinFitter::addConstraints( TAbsFitConstraint* c1, TAbsFitConstraint* c2, TAbsFitConstraint* c3, 
				 TAbsFitConstraint* c4, TAbsFitConstraint* c5, TAbsFitConstraint* c6,
				 TAbsFitConstraint* c7, TAbsFitConstraint* c8, TAbsFitConstraint* c9) {
  // add many constraints

  resetStatus();

  if ( c1 != 0 ) _constraints.push_back( c1 );
  if ( c2 != 0 ) _constraints.push_back( c2 );
  if ( c3 != 0 ) _constraints.push_back( c3 );
  if ( c4 != 0 ) _constraints.push_back( c4 );
  if ( c5 != 0 ) _constraints.push_back( c5 );
  if ( c6 != 0 ) _constraints.push_back( c6 );
  if ( c7 != 0 ) _constraints.push_back( c7 );
  if ( c8 != 0 ) _constraints.push_back( c8 );
  if ( c9 != 0 ) _constraints.push_back( c9 );

  countMeasParams();

}

void TKinFitter::setParamMeas( Int_t index ) {
  // Declares a parameter as measured

  if ( index >= 0 && index < _paramMeasured.size() ) {
    resetStatus();
    _paramMeasured[index] = true;
    countMeasParams();
    countUnmeasParams();
  } 

}

void TKinFitter::setParamMeas( TAbsFitParticle* particle, Int_t index ) {
  // Declares parameter of particle as measured

  Int_t parIndex = getParameterIndex( particle, index );
  if ( parIndex >= 0 && parIndex < _paramMeasured.size() ) {
    resetStatus();
    _paramMeasured[parIndex] = true;
    countMeasParams();
    countUnmeasParams();
  }

}

void TKinFitter::setParamUnmeas( Int_t index ) {
  // Declares a parameter as unmeasured

  if ( index >= 0 && index < _paramMeasured.size() ) {
    resetStatus();
    _paramMeasured[index] = false;
    countMeasParams();
    countUnmeasParams();
  } 

}

void TKinFitter::setParamUnmeas( TAbsFitParticle* particle, Int_t index ) {
  // Declares parameter of particle as unmeasured

  Int_t parIndex = getParameterIndex( particle, index );
  if ( parIndex >= 0 && parIndex < _paramMeasured.size() ) {
    resetStatus();
    _paramMeasured[parIndex] = false;
    countMeasParams();
    countUnmeasParams();
  }

}

Bool_t TKinFitter::getParamMeas( Int_t index ) {
  // Returns whether parameter with given index is measured or unmeasured

  if ( index >= 0 && index < _paramMeasured.size() ) {
    return (Bool_t) _paramMeasured[index];
  } else {
    cout << "TKinFitter::getParamMeas:: index " << index << " out of bounds." << endl;
    return false;
  }

}

Bool_t TKinFitter::getParamMeas( TAbsFitParticle* particle, Int_t index ) {
  // Returns whether parameter with given particle is measured or unmeasured

  Int_t parIndex = getParameterIndex( particle, index );
  if ( parIndex >= 0 && parIndex < _paramMeasured.size() ) {
    return (Bool_t) _paramMeasured[parIndex];
  } else {
    cout << "TKinFitter::getParamMeas:: parameter with index " << index << " of particle " 
	 << particle->GetName() << " not found." << endl;
    return false;
  }

}

Int_t TKinFitter::getParticleIndex( TAbsFitParticle* particle ) {
  // Return index of given particle. If particle has not been added
  // to kinematic fitter -1 is returned

  for ( UInt_t index = 0; index < _particles.size(); index++ ) {
    TAbsFitParticle* aParticle = _particles[index];
    if ( aParticle == particle ) {
      return index;
    }
  }
  return -1;

}

Int_t TKinFitter::getParameterIndex( TAbsFitParticle* particle, Int_t aPar ) {
  // Returns index of parameter "aPar" of given particle

  // Get index of particle
  Int_t pIndex = getParticleIndex( particle );
  if ( pIndex == -1 ) return -1;

  // Check index
  Int_t nPar = particle->getNPar();
  if (aPar < 0 || aPar >= nPar) return -1;

  // Calculate index of particle
  Int_t parIndex = aPar;
  for ( UInt_t index = 0; index < pIndex; index++ ) {
    TAbsFitParticle* aParticle = _particles[index];
    parIndex += aParticle->getNPar();
  }

  return parIndex;

}

void TKinFitter::setVerbosity( Int_t verbosity ) { 
  // Set verbosity of the fitter:
  // 0: quiet
  // 1: print information before and after the fit
  // 2: print output for every iteration of the fit
  // 3: print debug information

  if ( verbosity < 0 ) verbosity = 0;
  if ( verbosity > 3 ) verbosity = 3;
  _verbosity = verbosity;

}


Int_t TKinFitter::fit() {
  // Perform the fit
  // Returns:
  // 0: converged
  // 1: not converged

  resetParams();
  resetStatus();

  Bool_t isConverged = false;
  Double_t prevF;
  Double_t currF = getF();
  Double_t prevS;
  Double_t currS = 0.;

  // Calculate covariance matrix V
  calcV();

  // print status
  if ( _verbosity >= 1 ) {
    print();
  }

  do {

    // Reset status to "RUNNING"
    _status = 10;

    // cout << "calcB();" << endl;
    calcB();
    // cout << "calcVB();" << endl;
    calcVB();
    if ( _nParA > 0 ) {
      // cout << "calcA();" << endl;
      calcA();
      // cout << "calcVA();" << endl;
      calcVA();
      // cout << "calcC32();" << endl;
      calcC32();
    }
    // cout << "calcC31();" << endl;
    calcC31();
    // cout << "calcC33();" << endl;
    calcC33();
    // cout << "calcC();" << endl;
    calcC();

    // Calculate corretion for a, y, and lambda
    // cout << "if ( _nParA > 0 ) calcDeltaA();" << endl;
    if ( _nParA > 0 ) calcDeltaA();
    // cout << "calcDeltaY();" << endl;
    calcDeltaY();
    // cout << "calcLambda();" << endl;
    calcLambda();
   
    if( _verbosity >= 3 ) printMatrices();

    if( _verbosity >= 2 ) {
      cout << endl << endl << endl << endl;
      print();   
      cout <<"---------" <<endl ;
      cout << endl << endl << endl << endl;
    }

    // Apply calculated corrections to measured and unmeasured particles
    applyDeltaAY();
    _nbIter++;
    
    //calculate F and S
    prevF = currF;
    currF = getF();
    prevS = currS;
    currS = getS();

    if( TMath::IsNaN(currF) ) {
      if(_verbosity>=2) cout << "The current value of F is NaN. Fit will be aborted." << endl;
      _status = -10;
    }
    if( TMath::IsNaN(currS) ) {
      if(_verbosity>=2) cout << "The current value of S is NaN. Fit will be aborted." << endl;
      _status = -10;
    }
    
    // If S or F are getting bigger reduce step width
//     Int_t nstep =0;
//     while ( currF >= prevF ) {
//       nstep++;
//       if (nstep < 6) {cout <<nstep <<" : currF: "<< currF << "\t , prevF: " << prevF << endl;}
//       //      cout << "Reducing step width ..." << endl;
//       _deltaA *= (1.-0.001);
//       _deltaY *= (1.-0.001);
// //       _lambda *= 0.00001;
// //       _lambdaT *= 0.00001;
//       applyDeltaAY();
//       currF = getF();
//       currS = getS();
//     }
    
    // Test convergence
    isConverged = converged(currF, prevS, currS);

  } while ( (! isConverged) && (_nbIter < _maxNbIter) && (_status!=-10) );

  // Calculate covariance matrices
  calcB();
  calcVB();
  if ( _nParA > 0 ) {
    calcA();
    calcVA();
    calcC21();
    calcC22();
    calcC32();
  }
  calcC11();
  calcC31();
  calcC33();
  calcVFit();
  applyVFit();

  if( _verbosity >= 3 ) printMatrices();

  // Set status information
  if (! isConverged ) {
    _status = 1;
  } else {
    _status = 0;
  }

  // print status
  if ( _verbosity >= 1 ) {
    print();
  }

  return _status;

}

void TKinFitter::setCovMatrix( TMatrixD &V ) {
  // Set the covariance matrix of the measured particles

  if ( (V.GetNrows() != _nParB) || (V.GetNcols() != _nParB) ) {
    cout << "Matrix needs to be a " << _nParB << "x" << _nParB << " matrix." << endl;
  } else {
    _V.ResizeTo( V );
    _V = V;
  }

}

Bool_t TKinFitter::calcV() {
  // Combines the covariance matrices of all measured particles

  _V.ResizeTo( _nParB, _nParB );
  _V.Zero();

  Int_t offsetMeasPar = -1;
  Int_t offsetPar = 0;

  for (UInt_t iP = 0; iP < _particles.size(); iP++) {

    // Get particle and its covariance matrix
    TAbsFitParticle* particle = _particles[iP];
    Int_t nParP = particle->getNPar();
    const TMatrixD* covMatrix =  particle->getCovMatrix();

    // Copy measured parameters
    Int_t vIndexX = offsetMeasPar;
    for (UInt_t iX = 0; iX < nParP; iX++) {

      Int_t vIndexY = offsetMeasPar;
      Bool_t measuredX =  (Bool_t) _paramMeasured[iX + offsetPar];
      if (measuredX) vIndexX++;

      for (UInt_t iY = 0; iY < nParP; iY++) {

	Bool_t measuredY =  (Bool_t) _paramMeasured[iY + offsetPar];
	if (measuredY) vIndexY++;

	// Only copy if both parameters are measured
	if ( measuredX && measuredY ) {
	  _V(vIndexX, vIndexY) = (*covMatrix)(iX, iY);
	}

      }
    }

    offsetPar += nParP;
    offsetMeasPar = vIndexX;

  }

  // Add parameters of constraints
  for (UInt_t iC = 0; iC < _constraints.size(); iC++) {

    // Get constraint and its covariance matrix
    TAbsFitConstraint* constraint = _constraints[iC];
    Int_t nParP = constraint->getNPar();
    
    if (nParP > 0) {
      const TMatrixD* covMatrix =  constraint->getCovMatrix();
      for (UInt_t iX = 0; iX < nParP; iX++) {
	for (UInt_t iY = 0; iY < nParP; iY++) {
	  _V(iX+offsetMeasPar+1, iY+offsetMeasPar+1) = (*covMatrix)(iX, iY);
	}
      }
      offsetMeasPar += nParP;
    }

  }

  _Vinv.ResizeTo( _V );
  _Vinv = _V;
  _Vinv.Invert();

  return true;

}

Bool_t TKinFitter::calcA() {
  // Calculate the Jacobi matrix of unmeasured parameters df_i/da_i
  // Row i contains the derivatives of constraint f_i. Column q contains
  // the derivative wrt. unmeasured parameter a_q.

  _A.ResizeTo( _constraints.size(), _nParA );

  for (UInt_t indexConstr = 0; indexConstr < _constraints.size(); indexConstr++) {

    Int_t offsetParam = 0;
    Int_t indexUnmeasParam = -1;
    for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {

      // Calculate matrix product  df/dP * dP/dy = (df/dr, df/dtheta, df/dphi, ...)
      TAbsFitParticle* particle = _particles[indexParticle];
      TMatrixD* derivParticle = particle->getDerivative();
      TMatrixD* derivConstr = _constraints[indexConstr]->getDerivative( particle );
      TMatrixD deriv( *derivConstr, TMatrixD::kMult, *derivParticle );

      // Copy columns with unmeasured parameters
      for (UInt_t indexParam = 0; indexParam < deriv.GetNcols(); indexParam++) {

	Bool_t measured =  (Bool_t) _paramMeasured[indexParam + offsetParam];
	if ( ! measured ) {
	  indexUnmeasParam++;
	  _A(indexConstr, indexUnmeasParam) = deriv(0, indexParam);
	}
      }
      offsetParam += deriv.GetNcols();

      delete derivParticle;
      delete derivConstr;

    }
  }

  // Calculate transposed matrix
  TMatrixD AT(TMatrixD::kTransposed, _A);
  _AT.ResizeTo( AT );
  _AT = AT;

  return true;

}

Bool_t TKinFitter::calcB() {
  // Calculate the Jacobi matrix of measured parameters df_i/da_i
  // Row i contains the derivatives of constraint f_i. Column q contains
  // the derivative wrt. a_q.

  _B.ResizeTo( _constraints.size(), _nParB );

  Int_t constrParamOffset = 0;

  for (UInt_t indexConstr = 0; indexConstr < _constraints.size(); indexConstr++) {

    Int_t offsetParam = 0;
    Int_t indexMeasParam = -1;

    // Copy particles Jacobi matrices
    for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {

      // Calculate matrix product  df/dP * dP/dy = (df/dr, df/dtheta, df/dphi, ...)
      TAbsFitParticle* particle = _particles[indexParticle];
      TMatrixD* derivParticle = particle->getDerivative();
      TMatrixD* derivConstr = _constraints[indexConstr]->getDerivative( particle );
      TMatrixD deriv( *derivConstr,  TMatrixD::kMult, *derivParticle );

//       if (indexConstr == 0) {
// 	cout << indexConstr << endl;
// 	cout << "derivConstr : " << endl;
// 	derivConstr->Print();
// 	cout << "derivParticle : " << endl;
// 	derivParticle->Print();
// 	cout << "derivConstr X derivParticle" << endl;
// 	deriv.Print();
//       }

      // Copy columns with measured parameters
      for (UInt_t indexParam = 0; indexParam < deriv.GetNcols(); indexParam++) {
	Bool_t measured =  (Bool_t) _paramMeasured[indexParam + offsetParam];
	if (measured) {
	  indexMeasParam++;
	  _B(indexConstr, indexMeasParam) = deriv(0, indexParam);
	}
      }
      offsetParam += deriv.GetNcols();

      delete derivParticle;
      delete derivConstr;

    }

    // Copy derivatives Jacobi matrices
    TAbsFitConstraint* constraint = _constraints[indexConstr];
    TMatrixD* deriv = constraint->getDerivativeAlpha();

    if (deriv != 0) {

      // Verbose output
      if (_verbosity >= 3) {
	cout << endl << "===  B deriv alpha: Constraint -> " 
	     << constraint->GetName() << "===" << endl;
	deriv->Print();
      }

      // Copy parameters
      for (Int_t indexParam = 0; indexParam < deriv->GetNcols(); indexParam++) {
	indexMeasParam++;
	_B( indexConstr, indexMeasParam+constrParamOffset) = (*deriv)(0, indexParam);
      }
      offsetParam += deriv->GetNcols();
      constrParamOffset += deriv->GetNcols();
      
      delete deriv;
    }

  }

  TMatrixD BT( TMatrixD::kTransposed,  _B );
  _BT.ResizeTo( BT );
  _BT = BT;

  return true;

}

Bool_t TKinFitter::calcVB() {
  // Calculate the matrix V_B = (B*V*B^T)^-1

  TMatrixD BV( _B, TMatrixD::kMult, _V );
  TMatrixD VBinv( BV, TMatrixD::kMult, _BT );
  _VBinv.ResizeTo( VBinv );
  _VBinv = VBinv;

  _VB.ResizeTo( _VBinv );
  _VB = _VBinv;
  _VB.Invert();

  return true;

}

Bool_t TKinFitter::calcVA() {
  // Calculate the matrix VA = (A^T*VB*A)

  TMatrixD ATVB( _AT, TMatrixD::kMult, _VB );
  TMatrixD VA(ATVB, TMatrixD::kMult, _A);
  _VA.ResizeTo( VA );
  _VA = VA;

  _VAinv.ResizeTo( _VA );
  _VAinv = _VA;
  _VAinv.Invert();

  return true;

}

Bool_t TKinFitter::calcC11() {
  // Calculate the matrix C11 = V^(-1) - V^(-1)*BT*VB*B*V^(-1) + V^(-1)*BT*VB*A*VA^(-1)*AT*VB*B*V^(-1)

  TMatrixD VBT( _V, TMatrixD::kMult, _BT );
  TMatrixD VBB( _VB, TMatrixD::kMult, _B );
  TMatrixD VBTVBB( VBT, TMatrixD::kMult, VBB );
  TMatrixD m2( VBTVBB,  TMatrixD::kMult, _V );

  _C11.ResizeTo( _V );
  _C11 = _V;
  _C11 -= m2;

  if ( _nParA > 0 ) {
    TMatrixD VBA( _VB, TMatrixD::kMult, _A );
    TMatrixD VBTVBA( VBT, TMatrixD::kMult, VBA );
    TMatrixD VAinvAT( _VAinv, TMatrixD::kMult, _AT );
    TMatrixD VBTVBAVAinvAT( VBTVBA, TMatrixD::kMult, VAinvAT );
    TMatrixD VBTVBAVAinvATVBB( VBTVBAVAinvAT, TMatrixD::kMult, VBB );
    TMatrixD m3( VBTVBAVAinvATVBB, TMatrixD::kMult, _V );
    _C11 += m3;
  }

  TMatrixD C11T( TMatrixD::kTransposed,  _C11 );
  _C11T.ResizeTo( C11T );
  _C11T = C11T;

  return true;

}

Bool_t TKinFitter::calcC21() {
  // Calculate the matrix  C21 = -VA^(-1)*AT*VB*B*V^(-1)

  TMatrixD VAinvAT( _VAinv, TMatrixD::kMult, _AT );
  TMatrixD VBB( _VB, TMatrixD::kMult, _B );
  TMatrixD VAinvATVBB( VAinvAT, TMatrixD::kMult, VBB );
  TMatrixD C21( VAinvATVBB, TMatrixD::kMult, _V );
  C21 *= -1.;
  _C21.ResizeTo( C21 );
  _C21 = C21;
  
  TMatrixD C21T( TMatrixD::kTransposed,  _C21 );
  _C21T.ResizeTo( C21T );
  _C21T = C21T;

  return true;

}

Bool_t TKinFitter::calcC22() {
  //  Calculate the matrix C22 = VA^(-1)

  _C22.ResizeTo( _VAinv );
  _C22 = _VAinv;

  TMatrixD C22T( TMatrixD::kTransposed,  _C22 );
  _C22T.ResizeTo( C22T );
  _C22T = C22T;

  return true;

}

Bool_t TKinFitter::calcC31() {
  // Calculate the matrix  C31 = VB*B*V - VB*A*VA^(-1)*AT*VB*B*V

  TMatrixD VbB(_VB, TMatrixD::kMult, _B);
  TMatrixD m1( VbB, TMatrixD::kMult, _V );
  _C31.ResizeTo( m1 );
  _C31 = m1;

  if ( _nParA > 0 ) {
    TMatrixD VbA(_VB, TMatrixD::kMult, _A);
    TMatrixD VAinvAT( _VAinv, TMatrixD::kMult, _AT );
    TMatrixD VbBV( VbB,  TMatrixD::kMult, _V );
    TMatrixD VbAVAinvAT(VbA, TMatrixD::kMult, VAinvAT); 
    TMatrixD m2(VbAVAinvAT, TMatrixD::kMult, VbBV);

    _C31 -= m2;
  }

  TMatrixD C31T( TMatrixD::kTransposed,  _C31 );
  _C31T.ResizeTo( C31T );
  _C31T = C31T;

  return true;

}

Bool_t TKinFitter::calcC32() {
  // Calculate the matrix  C32 = VB*A*VA^(-1)

  TMatrixD VbA( _VB, TMatrixD::kMult, _A );
  TMatrixD C32( VbA, TMatrixD::kMult, _VAinv );

  _C32.ResizeTo( C32 );
  _C32 = C32;

  TMatrixD C32T( TMatrixD::kTransposed,  _C32 );
  _C32T.ResizeTo( C32T );
  _C32T = C32T;

  return true;

}

Bool_t TKinFitter::calcC33() {
  // Calculate the matrix C33 = -VB + VB*A*VA^(-1)*AT*VB

  _C33.ResizeTo( _VB );
  _C33 = _VB;
  _C33 *= -1.;

  if ( _nParA > 0 ) {
    TMatrixD VbA(_VB, TMatrixD::kMult, _A );
    TMatrixD VAinvAT( _VAinv, TMatrixD::kMult, _AT );
    TMatrixD VbAVAinvAT( VbA, TMatrixD::kMult, VAinvAT );
    TMatrixD C33( VbAVAinvAT,  TMatrixD::kMult, _VB );
    _C33 += C33;
  }

  TMatrixD C33T( TMatrixD::kTransposed,  _C33 );
  _C33T.ResizeTo( C33T );
  _C33T = C33T;

  return true;
}

Bool_t TKinFitter::calcC() {
  // Calculate the matrix c = A*deltaAStar + B*deltaYStar - fStar

  int offsetParam = 0;

  // calculate delta(a*), = 0 in the first iteration
  TMatrixD deltaastar( 1, 1 );
  if ( _nParA > 0 ) {

    Int_t indexUnmeasParam = 0;
    deltaastar.ResizeTo( _nParA, 1 );
    for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {
    
      TAbsFitParticle* particle = _particles[indexParticle];
      const TMatrixD* astar = particle->getParCurr();
      const TMatrixD* a = particle->getParIni();
      TMatrixD deltaastarpart(*astar);
      deltaastarpart -= *a;
      
      // Copy unmeasured parameters
      for (int indexParam = 0; indexParam < deltaastarpart.GetNrows(); indexParam++) {
	Bool_t measured =  (Bool_t) _paramMeasured[indexParam + offsetParam];
	if (!measured) {
	  deltaastar(indexUnmeasParam, 0) = deltaastarpart(indexParam, 0);
	  indexUnmeasParam++;
	}
      }
      offsetParam += deltaastarpart.GetNrows();
      
    }

    if ( _verbosity >= 3 ) {
      cout << "  ==== deltaastar =====" << endl;
      deltaastar.Print();
      cout << endl;
    }

  }

  // calculate delta(y*), = 0 in the first iteration
  TMatrixD deltaystar( _nParB, 1 );
  offsetParam = 0;
  Int_t indexMeasParam = 0;

  for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {
    TAbsFitParticle* particle = _particles[indexParticle];
    const TMatrixD* ystar = particle->getParCurr();
    const TMatrixD* y = particle->getParIni();
    TMatrixD deltaystarpart(*ystar);
    deltaystarpart -= *y;

    // Copy measured parameters
    for (int indexParam = 0; indexParam < deltaystarpart.GetNrows(); indexParam++) {
      Bool_t measured =  (Bool_t) _paramMeasured[indexParam + offsetParam];
      if (measured) {
	deltaystar(indexMeasParam, 0) = deltaystarpart(indexParam, 0);
	indexMeasParam++;
      }
    }
    offsetParam += deltaystarpart.GetNrows();
  }

  for (UInt_t iC = 0; iC < _constraints.size(); iC++) {
    TAbsFitConstraint* constraint = _constraints[iC];
    if ( constraint->getNPar() > 0 ) {

      const TMatrixD* alphastar = constraint->getParCurr();
      const TMatrixD* alpha = constraint->getParIni();

      if ( _verbosity >= 3 ) {
	cout << "  ==== alphastar =====" << endl;
	alphastar->Print();
	cout << endl;
	cout << "  ==== alpha =====" << endl;
	alpha->Print();
	cout << endl;
      }

      TMatrixD deltaalphastarpart(*alphastar);
      deltaalphastarpart -= *alpha;

      for (int indexParam = 0; indexParam < deltaalphastarpart.GetNrows(); indexParam++) {
	deltaystar(indexMeasParam, 0) = deltaalphastarpart(indexParam, 0);
	indexMeasParam++;
      }
      offsetParam += deltaalphastarpart.GetNrows();
    }
  }

  if ( _verbosity >= 3 ) {
    cout << "  ==== deltaystar =====" << endl;
    deltaystar.Print();
    cout << endl;
  }

  // calculate f*
  TMatrixD fstar( _constraints.size(), 1 );
  for (UInt_t indexConstr = 0; indexConstr < _constraints.size(); indexConstr++) {
    fstar( indexConstr, 0 ) = _constraints[indexConstr]->getCurrentValue();
  }

  if ( _verbosity >= 3 ) {
    cout << "  ==== fstar =====" << endl;
    fstar.Print();
    cout << endl;
  }

  // calculate c
  _c.ResizeTo( fstar );
  _c = fstar;
  _c *= (-1.);
  TMatrixD Bdeltaystar( _B, TMatrixD::kMult, deltaystar );
  _c += Bdeltaystar;
  if ( _nParA ) {
    TMatrixD Adeltaastar( _A, TMatrixD::kMult, deltaastar );
    _c += Adeltaastar;
  }

  if ( _verbosity >= 3 ) {
    cout << "  ==== c =====" << endl;
    _c.Print();
  }

  return true;

}

Bool_t TKinFitter::calcDeltaA() {
  // Calculate the matrix deltaA = C32T * c
  // (corrections to unmeasured parameters)

  TMatrixD deltaA( _C32T, TMatrixD::kMult, _c );
  _deltaA.ResizeTo( deltaA );
  _deltaA = deltaA;

  return true;

}

Bool_t TKinFitter::calcDeltaY() {
  // Calculate the matrix deltaY = C31T * c 
  // (corrections to measured parameters)

  TMatrixD deltaY( _C31T, TMatrixD::kMult, _c );
  _deltaY.ResizeTo( deltaY );
  _deltaY = deltaY;

  return true;

}

Bool_t TKinFitter::calcLambda() {
  // Calculate the matrix Lambda = C33 * c 
  // (Lagrange Multipliers)

  TMatrixD lambda( _C33,  TMatrixD::kMult, _c );
  _lambda.ResizeTo( lambda );
  _lambda = lambda;

  TMatrixD lambdaT( TMatrixD::kTransposed,  _lambda );
  _lambdaT.ResizeTo( lambdaT );
  _lambdaT = lambdaT;

  return true;

}

Bool_t TKinFitter::calcVFit() {
  // Calculate the covariance matrix of fitted parameters
  //
  // Vfit(y) = ( C11  C21T )
  //     (a)   ( C21  C22  )
  //
  // Vfit(lambda) = (-C33)
  
  // Calculate covariance matrix of lambda
  _lambdaVFit.ResizeTo( _C33 );
  _lambdaVFit = _C33;
  _lambdaVFit *= -1.;


  // Calculate combined covariance matrix of y and a
  Int_t nbRows = _C11.GetNrows();
  Int_t nbCols = _C11.GetNcols();
  if ( _nParA > 0 ) {
    nbRows += _C21.GetNrows();
    nbCols += _C21T.GetNcols();
  }
  _yaVFit.ResizeTo( nbRows, nbCols );

  for (int iRow = 0; iRow < nbRows; iRow++) {
    for (int iCol = 0; iCol < nbCols; iCol++) {

      if (iRow >= _C11.GetNrows()) {
	if (iCol >= _C11.GetNcols()) {
	  _yaVFit(iRow, iCol) = _C22(iRow-_C11.GetNrows(), iCol-_C11.GetNcols());
	} else {
	  _yaVFit(iRow, iCol) = _C21(iRow-_C11.GetNrows(), iCol);
	}
      } else {
	if (iCol >= _C11.GetNcols()) {
	  _yaVFit(iRow, iCol) = _C21T(iRow, iCol-_C11.GetNcols());
	} else {
	  _yaVFit(iRow, iCol) = _C11(iRow, iCol);
	}
      }

    }
  }

  return true;

}

void TKinFitter::applyVFit() {
  // apply fit covariance matrix to measured and unmeasured  particles 
  // and constraints parameters

  Int_t offsetUnmeasured = _nParA;
  Int_t offsetMeasured = 0;
  Int_t offsetParam = 0;

  // Copy particles parameters
  for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {
    TAbsFitParticle* particle = _particles[indexParticle];
    Int_t nbParams = particle->getNPar();
    TArrayI params( nbParams );

    // Determine measured and unmeasured parameters
    for (Int_t indexParam = 0; indexParam < nbParams; indexParam++) {

      Bool_t measured =  (Bool_t) _paramMeasured[indexParam + offsetParam];
      if (measured) {
	params[indexParam] = offsetMeasured;
	offsetMeasured++;
      } else {
	params[indexParam] = offsetUnmeasured;
	offsetUnmeasured++;
      }

    }

    // Copy parameters
    TMatrixD vfit( nbParams, nbParams );
    for (Int_t c = 0; c < nbParams; c++) {
      for (Int_t r = 0; r < nbParams; r++) {
	vfit(c, r) = _yaVFit( params[c], params[r] );
      }
    }
    particle->setCovMatrixFit( &vfit );

    offsetParam += nbParams;

  }

  // Copy constraints parameters
  for (UInt_t indexConstraint = 0; indexConstraint < _constraints.size(); indexConstraint++) {
    TAbsFitConstraint* constraint = _constraints[indexConstraint];
    Int_t nbParams = constraint->getNPar();
    TMatrixD vfit( nbParams, nbParams );
    if (nbParams > 0) {
      TMatrixD vfit( nbParams, nbParams );
      for (Int_t c = 0; c < nbParams; c++) {
	for (Int_t r = 0; r < nbParams; r++) {
	  vfit(r, c) = _yaVFit( r + offsetMeasured, c + offsetMeasured );
	}
      }
      constraint->setCovMatrixFit( &vfit );
      offsetMeasured += nbParams;
    }
  }

}

Bool_t TKinFitter::applyDeltaAY() {
  // apply corrections to measured and unmeasured unmeasured 
  // particles and constraints

  Int_t offsetParam = 0;
  Int_t indexY = 0;
  Int_t indexA = 0;
  
  // apply corrections to particles
  for (UInt_t indexParticle = 0; indexParticle < _particles.size(); indexParticle++) {

    TAbsFitParticle* particle = _particles[indexParticle];
    Int_t nbParams = particle->getNPar();
    TMatrixD params( nbParams, 1 );
    for (Int_t index = 0; index < nbParams; index++) {

      Bool_t measured =  (Bool_t) _paramMeasured[index + offsetParam];
      if (measured) {
	params(index, 0) = _deltaY(indexY, 0);
	indexY++;
      } else {
	params(index, 0) = _deltaA(indexA, 0);
	indexA++;
      }

    }
    particle->applycorr( &params );
    offsetParam+=nbParams;

  }

  // apply corrections to constraints
  for (UInt_t indexConstraint = 0; indexConstraint < _constraints.size(); indexConstraint++) {

    TAbsFitConstraint* constraint = _constraints[indexConstraint];
    Int_t nbParams = constraint->getNPar();
    if ( nbParams > 0 ) {
      
      TMatrixD params( nbParams, 1 );
      for (Int_t index = 0; index < nbParams; index++) {
	params(index, 0) = _deltaY(indexY, 0);
	indexY++;
      }
      constraint->applyDeltaAlpha( &params );
      
    }
  }

  
  return true;
  
}

Double_t TKinFitter::getF() {
  // calculate current absolut value of constraints
  // F = Sum[ Abs(f_k( aStar, yStar)) ] 

  Double_t F = 0.;
  for (UInt_t indexConstr = 0; indexConstr < _constraints.size(); indexConstr++) {
    F += TMath::Abs( _constraints[indexConstr]->getCurrentValue() );
  }
  
  return F;

}

Double_t TKinFitter::getS() {
  // calculate current value of Chi2
  // S = deltaYT * V^-1 * deltaY

  Double_t S = 0.;
  if ( _nbIter > 0 ) {
    TMatrixD deltaYTVinv(_deltaY, TMatrixD::kTransposeMult, _Vinv);
    TMatrixD S2(deltaYTVinv, TMatrixD::kMult, _deltaY);
    S = S2(0,0);
  }

  return S;

}

Bool_t TKinFitter::converged( Double_t F, Double_t prevS, Double_t currS ) {
  // check whether convergence criteria are fulfilled

  Bool_t isConverged = false;
  
  // calculate F, delta(a) and delta(y) already applied
  isConverged = (F < _maxF);

  // Calculate current Chi^2 and delta(S)
  Double_t deltaS = currS - prevS;
  isConverged = isConverged && (TMath::Abs(deltaS) < _maxDeltaS);

  return isConverged;

}

TString TKinFitter::getStatusString() {

  TString statusstring = "";

  switch ( _status ) {
    
  case -1: {
    statusstring = "NO FIT PERFORMED";
    break;
  }
  case 10: {
    statusstring = "RUNNING";
    break;
  }
  case 0: {
    statusstring = "CONVERGED";
    break;
  }
  case 1: {
    statusstring = "NOT CONVERGED";
    break;
  }
  case -10: {
    statusstring = "ABORTED";
    break;
  }
  default:{
    statusstring = "NOT DEFINED";
    break;
  }
  }
    
  return statusstring;

}

void TKinFitter::print() {

  cout << endl << endl;
  cout << setprecision( 4 );
  cout << "Status: " << getStatusString();
  cout << "   F=" << getF() << "   S=" << getS() << "   N=" << _nbIter << "   NDF=" << getNDF() << endl;
  cout << "particles:" << endl;
  Int_t parIndex = 0;
  for (UInt_t iP = 0; iP < _particles.size(); iP++) {
    TAbsFitParticle* particle = _particles[iP];
    Int_t nParP = particle->getNPar();
    const TMatrixD* par = particle->getParCurr();
    const TMatrixD* parIni = particle->getParIni();
    const TMatrixD* covP = particle->getCovMatrix();
    cout << setw(3) << setiosflags(ios::right) << iP;
    cout << setw(15) << setiosflags(ios::right) << particle->GetName();
    cout << setw(3) << " ";
    for (int iPar = 0; iPar < nParP; iPar++) {
      if (iPar > 0) {
	cout << setiosflags(ios::right) << setw(21) << " ";
      }
      TString colstr = "";
      colstr += parIndex;
      colstr += ":";
      cout << setw(4) << colstr;
      cout << setw(2) << " ";   
      cout << setiosflags(ios::left) << setiosflags(ios::scientific) << setprecision(3);
      cout << setw(15) << (*par)(iPar, 0);
      if(_nbIter > 0 && _status < 10) {
	cout << setw(15) << TMath::Sqrt( _yaVFit(iPar, iPar) );
      } else {
	cout << setw(15) << " ";
      }
      Bool_t measured =  (Bool_t) _paramMeasured[parIndex];
      if (!measured) {
	cout << setw(15) << " ";
	cout << setw(13) << "   (unmeas.)";
      } else {
	cout << setw(15) << TMath::Sqrt( (*covP)(iPar, iPar) );
      }
      cout << endl;
      parIndex++;
    }
    if(_nbIter > 0 && _status < 10) {
      cout << "covariance matrix (fit): " << endl;
      const TMatrixD* covMatrix = particle->getCovMatrixFit();
      covMatrix->Print();
    }
    //    particle->print();
  }

  cout << endl;
  cout << "constraints: " << endl;
  for (UInt_t indexConstr = 0; indexConstr < _constraints.size(); indexConstr++) {
    _constraints[indexConstr]->print();
  }
  cout << endl;

}

void TKinFitter::printMatrices() {

  cout << endl << "================   A   ================" << endl;
  _A.Print();

  cout << endl << "================   B   ================" << endl;
  _B.Print();

  cout << endl << "================  V   ================" << endl;
  _V.Print();

  cout << endl << "================  Vinv   ================" << endl;
  _Vinv.Print();

  cout << endl << "================  VB   ================" << endl;
  _VB.Print();

  cout << endl << "================  VBinv   ================" << endl;
  _VBinv.Print();

  cout << endl << "================  VA   ================" << endl;
  _VA.Print();

  cout << endl << "================  VAinv   ================" << endl;
  _VAinv.Print();

  cout << endl << "================  C   ================" << endl;
  _c.Print();

  cout << endl << "================  C11   ================" << endl;
  _C11.Print();

  cout << endl << "================  C21   ================" << endl;
  _C21.Print();

  cout << endl << "================  C22   ================" << endl;
  _C22.Print();

  cout << endl << "================  C31   ================" << endl;
  _C31.Print();

  cout << endl << "================  C32   ================" << endl;
  _C32.Print();

  cout << endl << "================  C33   ================" << endl;
  _C33.Print();

  cout << endl << "================  deltaY   ================" << endl;
  _deltaY.Print();

  cout << endl << "================  deltaA   ================" << endl;
  _deltaA.Print();


  cout << endl << "================  lambda   ================" << endl;
  _lambda.Print();

  cout << endl << "================  lambdaVFit   ================" << endl;
  _lambdaVFit.Print();

  cout << endl << "================  yaVFit   ================" << endl;
  _yaVFit.Print();

}
