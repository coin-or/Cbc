// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/22/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------

#ifndef CglLandPValidator_H
#define CglLandPValidator_H
#include "OsiSolverInterface.hpp"
#include "CglParam.hpp"
#include <vector>

/** constants describing rejection codes*/
//[5] = {"Accepted", "violation too small", "small coefficient too small", "big dynamic","too dense"}


namespace LAP
{

/** Class to validate or reject a cut */
class CGLLIB_EXPORT Validator
{
public:
    /** Reasons for rejecting a cut */
    enum RejectionsReasons
    {
        NoneAccepted=0 /**Cut was accepted*/,
        SmallViolation /** Violation of the cut is too small */,
        SmallCoefficient /** There is a small coefficient we can not get rid off.*/,
        BigDynamic /** Dynamic of coefficinet is too important. */,
        DenseCut/**cut is too dense */,
        EmptyCut/**After cleaning cut has become empty*/,
        DummyEnd/** dummy*/
    };

    /** Constructor with default values */
    Validator(double maxFillIn = 1.,
              double maxRatio = 1e8,
              double minViolation = 0,
              bool scale = false,
              double rhsScale = 1);

    /** Clean an OsiCut */
    int cleanCut(OsiRowCut & aCut, const double * solCut,const OsiSolverInterface &si, const CglParam & par,
                 const double * colLower, const double * colUpper);
    /** Clean an OsiCut by another method */
    int cleanCut2(OsiRowCut & aCut, const double * solCut, const OsiSolverInterface &si, const CglParam & par,
                  const double * colLower, const double * colUpper);
    /** Call the cut cleaner */
    int operator()(OsiRowCut & aCut, const double * solCut,const OsiSolverInterface &si, const CglParam & par,
                   const double * colLower, const double * colUpper)
    {
        return cleanCut(aCut, solCut, si, par, colLower, colUpper);
    }
    /** @name set functions */
    /** @{ */
    void setMaxFillIn(double value)
    {
        maxFillIn_ = value;
    }
    void setMaxRatio(double value)
    {
        maxRatio_ = value;
    }
    void setMinViolation(double value)
    {
        minViolation_ = value;
    }

    void setRhsScale(double v)
    {
        rhsScale_ = v;
    }
    /** @} */
    /** @name get functions */
    /** @{ */
    double getMaxFillIn()
    {
        return maxFillIn_;
    }
    double getMaxRatio()
    {
        return maxRatio_;
    }
    double getMinViolation()
    {
        return minViolation_;
    }
    /** @} */

    const char* failureString(RejectionsReasons code) const
    {
        return rejections_[static_cast<int> (code)];
    }
    const char* failureString(int code) const
    {
        return rejections_[ code];
    }
    int numRejected(RejectionsReasons code)const
    {
        return numRejected_[static_cast<int> (code)];
    }
    int numRejected(int code)const
    {
        return numRejected_[ code];
    }
private:
    /** max percentage of given formulation fillIn should be accepted for cut fillin.*/
    double maxFillIn_;
    /** max ratio between smallest and biggest coefficient */
    double maxRatio_;
    /** minimum violation for accepting a cut */
    double minViolation_;
    /** Do we do scaling? */
    //bool scale_;
    /** Scale of right-hand-side.*/
    double rhsScale_;
    /** Strings explaining reason for rejections */
    static const char* rejections_[DummyEnd];
    /** Number of cut rejected for each of the reasons.*/
    std::vector<int> numRejected_;
};

}/* Ends namespace LAP.*/
#endif
