// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#include "CglLandPMessages.hpp"
#include <cstring>
#define REMOVE_LOG 0
namespace LAP
{

#ifdef LAP_ADD_MSG
#error "Macro ADD_MSG already defined"    
#endif
#define LAP_ADD_MSG(Id,Type,Level,MSG) addMessage(Id, CoinOneMessage( Type(Id), Level, MSG))
inline int std_m(int n)
{
    return 1 + n;
}
inline int warn_m(int n)
{
    return  3000 + n;
}
inline int err_m(int n)
{
    return n + 6000;
}


LandPMessages::LandPMessages()
        :
        CoinMessages(DUMMY_END)
{
    strcpy(source_,"Lap");
    LAP_ADD_MSG(Separating,std_m,3+REMOVE_LOG,"Starting separation on variable %d, initial depth of cut %f");
    LAP_ADD_MSG(FoundImprovingRow, std_m,4,
                "Found improving row (leaving variable). Row %d (basic var %d), "
                "leaving status %d, sign of gamma %d, reduced cost %f");
    LAP_ADD_MSG(FoundBestImprovingCol, std_m, 4,
                " Found best improvement (entering variable). Var %d, "
                "value of gamma %f, expected depth of next cut %f");
    LAP_ADD_MSG(WarnFailedBestImprovingCol, err_m, 3,
                "Failed to find an improving entering variable while reduced cost was %f, "
                "depth of current cut %f, best cut depth with pivot %f");
    //Log line is cut number time pivot number,  cut depth, leaving, incoming gamma degenerate
    LAP_ADD_MSG(LogHead, std_m, 3+REMOVE_LOG,
                "Pivot no \t cut depth \t leaving var \t incoming var \t direction \t gamma \t degenerate");
    LAP_ADD_MSG(PivotLog, std_m, 3+REMOVE_LOG,
                "%8d\t %9f\t %11d \t %11d \t %11d \t %8f \t %12d \t %.5g \t %11d");
    LAP_ADD_MSG(FinishedOptimal, std_m, 2,
                "Found optimal lift-and-project cut, depth %f number of pivots performed %d");
    LAP_ADD_MSG(HitLimit, std_m, 2, "Stopping lift-and-project optimization hit %s limit. Number of pivots %d");
    LAP_ADD_MSG(WarnBadSigmaComputation, err_m, 1,
                "Cut depth after pivot is not what was expected by computations before, difference %.15f");
    LAP_ADD_MSG(WarnBadRowComputation, err_m, 1,
                "Row obtained after pivot is not what was expected (distance between the two %f in norm inf).");
    LAP_ADD_MSG(WarnGiveUpRow, err_m,1,"Limit of %d negative reduced costs with no strict improvement");
    LAP_ADD_MSG(PivotFailedSigmaUnchanged, err_m, 1,
                "A pivot failed to be performed (probably refactorization was performed) but sigma is unchanged continue...");
    LAP_ADD_MSG(PivotFailedSigmaIncreased, err_m,
                1,"A pivot failed to be performed, and sigma has changed exit without generating cut");
    LAP_ADD_MSG(FailedSigmaIncreased, err_m,1,"Cut violation has increased in last pivot");
    LAP_ADD_MSG
    (WarnBadRhsComputation,
     err_m,1,
     "rhs obtained  after pivot is not what was expected (distance between the two %f).");
    LAP_ADD_MSG(WarnFailedPivotTol,
                err_m, 2,"All pivots are below tolerance");
    LAP_ADD_MSG(WarnFailedPivotIIf,
                err_m, 2,"There is no possible pivot within tolerance (every pivot make rhs for current row %f too close to integer feasibility");
    LAP_ADD_MSG(NumberNegRc, std_m, 4,
                "Number of rows with negative reduced cost %i");
    LAP_ADD_MSG(NumberZeroRc, std_m, 4,
                "Number of rows with zero reduced cost %i");
    LAP_ADD_MSG(NumberPosRc, std_m, 4,
                "Number of rows with positive reduced cost %i");
    LAP_ADD_MSG(WeightsStats, std_m, 2,
                "Maximal weight %g minimal weight %g");
    LAP_ADD_MSG(RoundStats, std_m, 1,
                "Separated %i cuts with %i pivots, source entered %i times, %i sigma increases.");
    LAP_ADD_MSG(CutStat, std_m, 1,
                "Separated cut %i with %i pivots, source entered %i times, %i sigma increases, %i potential cycles.%g");
}
}

