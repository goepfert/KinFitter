// -*- C++ -*-
// This linkdef file contains all "pragma link" for CLHEP inclusion
// into root including non-member operators and functions
// of Vector, Matrix, DiagMatrix and SymMatrix:
#ifdef __CINT__
// ##################################################
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// ################## Functions #####################

#pragma link C++ class TAbsFitConstraint;
#pragma link C++ class TFitConstraintEp;
#pragma link C++ class TFitConstraintEpGaus;
#pragma link C++ class TFitConstraintPt;
#pragma link C++ class TFitConstraintM;
#pragma link C++ class TFitConstraintMGaus;
#pragma link C++ class TFitConstraintMBW;
#pragma link C++ class TAbsFitParticle;
#pragma link C++ class TFitParticlePxPyPz;
#pragma link C++ class TFitParticlePxPyPzE;
#pragma link C++ class TFitParticlePxPyPzM;
#pragma link C++ class TFitParticleRelPxPyPz;
#pragma link C++ class TFitParticleRelPxPyPzE;
#pragma link C++ class TFitParticleRelPxPyPzM;
#pragma link C++ class TFitParticlePThetaPhi;
#pragma link C++ class TFitParticlePThetaPhiE;
#pragma link C++ class TFitParticlePThetaPhiM;
#pragma link C++ class TFitParticlePInvThetaPhi;
#pragma link C++ class TFitParticlePInvThetaPhiE;
#pragma link C++ class TFitParticlePInvThetaPhiM;
#pragma link C++ class TFitParticlePtEtaPhi;
#pragma link C++ class TFitParticlePtEtaPhiE;
#pragma link C++ class TFitParticlePtEtaPhiM;
#pragma link C++ class TFitParticleRelE;
#pragma link C++ class TFitParticleRelPtEtaPhi;
#pragma link C++ class TFitParticleRelPtEtaPhiE;
#pragma link C++ class TFitParticleRelPtEtaPhiM;
#pragma link C++ class TKinFitter;
#pragma link C++ class TToyGentT;
#pragma link C++ class TGetErrorMatrix;

#endif
