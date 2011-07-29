#ifndef TToyGentT_h
#define TToyGentT_h

/** 
    @class TToyGentT
    class for generating Toy MC ttbar events
    @author Thomas Goepfert
*/

#include <vector>

#include "TMatrixD.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

#include "KinFitter/TKinFitter.h"

using namespace std;

class TAbsFitParticle;
class TAbsFitConstraint;

class TToyGentT : public TNamed {

public :

  TToyGentT( const TAbsFitParticle* Whad_d1, 
	     const TAbsFitParticle* Whad_d2, 
	     const TAbsFitParticle* Tophad_d2, 
	     const TAbsFitParticle* Wlep_d1, 
	     const TAbsFitParticle* Wlep_d2, 
	     const TAbsFitParticle* Toplep_d2);
  TToyGentT(){};
  virtual ~TToyGentT();
  Bool_t doToyExperiments( Int_t nbExperiments = 1000 );

  TH1D* _histStatus;
  TH1D* _histNIter;
  TH1D* _histPChi2;
  TH1D* _histChi2;

  TH1D* _histMWhadTrue;
  TH1D* _histMWhadSmear;
  TH1D* _histMWhadFit;
  TH1D* _histMTophadTrue;
  TH1D* _histMTophadSmear;
  TH1D* _histMTophadFit;

  TH1D* _histMWlepTrue;
  TH1D* _histMWlepSmear;
  TH1D* _histMWlepFit;
  TH1D* _histMToplepTrue;
  TH1D* _histMToplepSmear;
  TH1D* _histMToplepFit;


  TObjArray _histsParTrue;
  TObjArray _histsParSmear;
  TObjArray _histsParFit;

  TObjArray _histsPull1;
  TObjArray _histsError1;
  TObjArray _histsDiff1;
  TObjArray _histsPull2;
  TObjArray _histsError2;
  TObjArray _histsDiff2;

  TObjArray _histsPull3;
  TObjArray _histsError3;
  TObjArray _histsDiff3;
   
  int _nPar;
  
  void setprintPartIni(Bool_t value) { _printPartIni = value; }
  void setprintConsIni(Bool_t value) { _printConsIni = value; }
  void setprintSmearedPartBefore(Bool_t value) { _printSmearedPartBefore = value; }
  void setprintPartAfter(Bool_t value) { _printPartAfter = value; } 
  void setprintConsBefore(Bool_t value) { _printConsBefore = value; }
  void setprintConsAfter(Bool_t value) { _printConsAfter = value; }

  void setCheckConstraintsTruth(Bool_t value) { _doCheckConstraintsTruth = value; }

  void setVerbosity(int verbosity){_verbosity=verbosity;};

protected:

  void smearParticles();

  void createHists(TKinFitter* fitter);

  void fillPull1();
  void fillPull2();
  void fillPar();
  void fillM();
  void fillConstraints();

private :
  
  vector<TAbsFitParticle*> _inimeasParticles;    //!< vector that contains all true measured particles
  vector<TAbsFitParticle*> _iniunmeasParticles;  //!< vector that contains all true unmeasured particles

  TAbsFitParticle* _iniWhad_d1;
  TAbsFitParticle* _iniWhad_d2;
  TAbsFitParticle* _iniTophad_d2;
  TAbsFitParticle* _iniWlep_d1;
  TAbsFitParticle* _iniWlep_d2;
  TAbsFitParticle* _iniToplep_d2;

  vector<TAbsFitParticle*> _measParticles;     //!< vector that contains all smeared measured particles
  vector<TAbsFitParticle*> _unmeasParticles;   //!< vector that contains all smeared unmeasured particles
  vector<TAbsFitConstraint*> _measConstraints; //!< constraints which contains measured parameters

  TAbsFitParticle* _Whad_d1;
  TAbsFitParticle* _Whad_d2;
  TAbsFitParticle* _Tophad_d2;
  TAbsFitParticle* _Wlep_d1;
  TAbsFitParticle* _Wlep_d2;
  TAbsFitParticle* _Toplep_d2;

  Bool_t _printPartIni;
  Bool_t _printConsIni;
  Bool_t _printSmearedPartBefore ;
  Bool_t _printConsBefore;
  Bool_t _printConsAfter;
  Bool_t _printPartAfter;
  Bool_t _doCheckConstraintsTruth;

  int _verbosity;

  ClassDef(TToyGentT, 1)
};

#endif
