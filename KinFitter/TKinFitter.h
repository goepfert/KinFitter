
using namespace std;

#ifndef TKINFITTER_H
#define TKINFITTER_H


#include <vector>
#include "TMatrixD.h"
#include "TNamed.h"

class TAbsFitParticle;
class TAbsFitConstraint;
class TH1D;

class TKinFitter : public TNamed {

public :

  TKinFitter();
  TKinFitter(const TString &name, const TString &title);  
  ~TKinFitter();
  void reset();
  void resetStatus();   

  Int_t fit();

  void addMeasParticle( TAbsFitParticle* particle );
  void addMeasParticles( TAbsFitParticle* p1, TAbsFitParticle* p2 = 0, TAbsFitParticle* p3 = 0, 
			 TAbsFitParticle* p4 = 0, TAbsFitParticle* p5 = 0, TAbsFitParticle* p6 = 0,
			 TAbsFitParticle* p7 = 0, TAbsFitParticle* p8 = 0, TAbsFitParticle* p9 = 0);
  void addUnmeasParticle( TAbsFitParticle* particle );
  void addUnmeasParticles( TAbsFitParticle* p1, TAbsFitParticle* p2 = 0, TAbsFitParticle* p3 = 0, 
			   TAbsFitParticle* p4 = 0, TAbsFitParticle* p5 = 0, TAbsFitParticle* p6 = 0,
			   TAbsFitParticle* p7 = 0, TAbsFitParticle* p8 = 0, TAbsFitParticle* p9 = 0);
  void addConstraint( TAbsFitConstraint* constraint );
  void addConstraints( TAbsFitConstraint* c1, TAbsFitConstraint* c2 = 0, TAbsFitConstraint* c3 = 0, 
		       TAbsFitConstraint* c4 = 0, TAbsFitConstraint* c5 = 0, TAbsFitConstraint* c6 = 0,
		       TAbsFitConstraint* c7 = 0, TAbsFitConstraint* c8 = 0, TAbsFitConstraint* c9 = 0);

  void setParamMeas( Int_t index );
  void setParamMeas( TAbsFitParticle* particle, Int_t index );
  void setParamUnmeas( Int_t index );
  void setParamUnmeas( TAbsFitParticle* particle, Int_t index );
  Bool_t getParamMeas( Int_t index );
  Bool_t getParamMeas( TAbsFitParticle* particle, Int_t index );

  Int_t getNDF() { return  (_constraints.size() - _nParA); }
  Int_t getNParA() { return _nParA; }                              // Number of unmeasured parameters
  Int_t getNParB() { return _nParB; }                              // Number of measured parameters
  Int_t getNConstraints() { return _constraints.size(); }          // Number of measured parameters
  void setMaxNbIter( Int_t maxNbIter ) { _maxNbIter = maxNbIter; } // Set maximum number of iterations possible
  Int_t getMaxNumberIter() { return _maxNbIter; }                  // Maximum number of iterations possible
  Int_t getNbIter() { return _nbIter; }                            // Number of iterations performed in the previous fit
  Int_t getStatus() { return _status; }                            // Status of the kinematic fitter

  void setMaxDeltaS( Double_t maxDeltaS ) { _maxDeltaS = TMath::Abs( maxDeltaS ); }
  Double_t getMaxDeltaS() { return _maxDeltaS; }
  void setMaxF( Double_t maxF ) { _maxF = TMath::Abs( maxF ); }
  Double_t getMaxF() { return _maxF; }
  const TMatrixD* getCovMatrix() { return &_V; }            // Covariance matrix of measured parameters
  void setCovMatrix( TMatrixD &V );                         // Set covariance matrix of measured parameters
  const TMatrixD* getCovMatrixFit() { return &_yaVFit; }    // Covariance matrix of the fit
  Double_t getS();
  Double_t getF();
  void setVerbosity( Int_t verbosity = 1 ); 
  Int_t getVerbosity( ) { return _verbosity; }

  Int_t getNParticles() { return _particles.size(); }                                // Number of particles
  const TAbsFitParticle* getParticle( Int_t index ) { return _particles[index]; }    // Particle with given index
  Int_t getParticleIndex( TAbsFitParticle* particle );                               // Index of given particle
  Int_t getParameterIndex( TAbsFitParticle* particle, Int_t aPar );                  // Index of given particles parameter
 
  void print();
  void printMatrices();

protected:

  Bool_t calcA();
  Bool_t calcB();
  Bool_t calcVA();
  Bool_t calcVB();
  Bool_t calcC();

  Bool_t calcC11();
  Bool_t calcC21();
  Bool_t calcC22();
  Bool_t calcC31();
  Bool_t calcC32();
  Bool_t calcC33();

  Bool_t calcDeltaA();
  Bool_t calcDeltaY();
  Bool_t calcLambda();

  Bool_t calcV();
  Bool_t calcVFit();

  Bool_t applyDeltaAY();
  void applyVFit();

  Bool_t converged(Double_t F, Double_t prevS, Double_t currS);

  TString getStatusString();
  void countMeasParams();
  void countUnmeasParams();
  void resetParams();

private :

  Int_t _maxNbIter;    // Maximum number of iterations
  Double_t _maxDeltaS; // Convergence criterium for deltaS
  Double_t _maxF;      // Convergence criterium for F
  Int_t _verbosity;    // Verbosty of the fitter 0: quiet, 1: print result, 2: print iterations, 3: print also matrices

  TMatrixD _A;      // Jacobi Matrix of unmeasured parameters
  TMatrixD _AT;     // Transposed Jacobi Matrix of unmeasured parameters
  TMatrixD _B;      // Jacobi Matrix of measured parameters
  TMatrixD _BT;     // Transposed Jacobi Matrix of measured parameters
  TMatrixD _V;      // Covariance matrix
  TMatrixD _Vinv;   // Inverse covariance matrix
  TMatrixD _VB;     // VB = ( B*V^(-1)*BT )^(-1)
  TMatrixD _VBinv;  // VB = ( B*V^(-1)*BT )
  TMatrixD _VA;     // VA = ( AT*VB*A )
  TMatrixD _VAinv;  // VA = ( AT*VB*A )^(-1)
  TMatrixD _c;      // Vector c = A*delta(a*) + B*delta(y*) - f*

  TMatrixD _C11;     // Matrix C11
  TMatrixD _C11T;    // Matrix C11T
  TMatrixD _C21;     // Matrix C21
  TMatrixD _C21T;    // Matrix C21T
  TMatrixD _C22;     // Matrix C22
  TMatrixD _C22T;    // Matrix C22T
  TMatrixD _C31;     // Matrix C31
  TMatrixD _C31T;    // Matrix C31T
  TMatrixD _C32;     // Matrix C32
  TMatrixD _C32T;    // Matrix C32T
  TMatrixD _C33;     // Matrix C33
  TMatrixD _C33T;    // Matrix C33T

  TMatrixD _deltaA;  // The correction vector deltaA for unmeasured particles of the last iteration
  TMatrixD _deltaY;  // The correction vector deltaY for measured particles of the last iteration
  TMatrixD _lambda;  // The column vector of Lagrange multiplicators (likelihood L = S + 2 sum_i lambda_i * f_i)
  TMatrixD _lambdaT; // The row vector of Lagrange multiplicators (likelihood L = S + 2 sum_i lambda_i * f_i)

  TMatrixD _lambdaVFit;   // Covariance matrix of lambda after the fit
  TMatrixD _yaVFit;       // Combined covariance matrix of y and a after the fit

  Int_t _nParA;     // Number of unmeasured parameters
  Int_t _nParB;     // Number of measured parameters

  vector<TAbsFitConstraint*> _constraints;    // vector with constraints
  vector<TAbsFitParticle*>   _particles;      // vector with particles
  vector<Bool_t>             _paramMeasured;  // holds a list of measured an unmeasured particle parameters

  Int_t _status;   // Status of the last fit;_
  Int_t _nbIter;   // number of iteration performed in the fit

  ClassDef(TKinFitter, 1) // Class to perform kinematic fit with non-linear constraints
};

#endif
