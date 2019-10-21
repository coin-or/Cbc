#ifndef MIP_SECCUTGEN_H
#define MIP_SECCUTGEN_H

#include <vector>
#include <string>
#include <iomanip>

#include "MIP_wrap.hh"

/// SEC cut generator for circuit
class MIP_SECCutGen : public MIP_wrapper::CutGen {
  MIP_SECCutGen() { }
  MIP_wrapper* pMIP=0;
public:
  MIP_SECCutGen( MIP_wrapper* pw ) : pMIP(pw) { }
  /// Say what type of cuts
  virtual int getMask() { return MIP_wrapper::MaskConsType_Lazy | MIP_wrapper::MaskConsType_Usercut; }
  /// The flow variables
  std::vector<MIP_wrapper::VarId> varXij;
  int nN=0;                        // N nodes
  /// returns error message if fails
  std::string validate() const;
  void generate(const MIP_wrapper::Output&, MIP_wrapper::CutInput&);
  void print( std::ostream& );
};

#endif // MIP_SECCUTGEN_H
