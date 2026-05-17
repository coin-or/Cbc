// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           LIF
//           CNRS, Aix-Marseille Universites
// Date:     02/23/08
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#ifndef CglLandPMessages_H
#define CglLandPMessages_H

#include "CoinMessage.hpp"
#include "CoinMessageHandler.hpp"

namespace LAP
{
/** Forward declaration of class to store extra debug data.*/
class DebugData;
/** Types of messages for lift-and-project simplex.*/
enum LAP_messages
{
    Separating,
    FoundImprovingRow,
    FoundBestImprovingCol,
    WarnFailedBestImprovingCol,
    LogHead,
    PivotLog,
    FinishedOptimal,
    HitLimit,
    NumberNegRc,
    NumberZeroRc,
    NumberPosRc,
    WeightsStats,
    WarnBadSigmaComputation,
    WarnBadRowComputation,
    WarnGiveUpRow,
    PivotFailedSigmaUnchanged,
    PivotFailedSigmaIncreased,
    FailedSigmaIncreased,
    WarnBadRhsComputation,
    WarnFailedPivotTol,
    WarnFailedPivotIIf,
    RoundStats,
    CutStat,
    DUMMY_END
};
/** Message handler for lift-and-project simplex. */
class LandPMessages : public CoinMessages
{
public:

    /** Constructor */
    LandPMessages();
};
}
#endif
