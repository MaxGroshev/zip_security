/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#include "zstd_compress_internal.h"
#include "zstd_lazy.h"
#include "../common/bits.h" /* ZSTD_countTrailingZeros64 */

#if !defined(ZSTD_EXCLUDE_GREEDY_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_LAZY_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_LAZY2_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_BTLAZY2_BLOCK_COMPRESSOR)

#define kLazySkippingStep 8


/*-*************************************
*  Binary Tree search
***************************************/

static
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
void ZSTD_updateDUBT(ZSTD_matchState_t* ms,
                const BYTE* ip, const BYTE* iend,
                U32 mls)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const hashTable = ms->hashTable;
    U32  const hashLog = cParams->hashLog;

    U32* const bt = ms->chainTable;
    U32  const btLog  = cParams->chainLog - 1;
    U32  const btMask = (1 << btLog) - 1;

    const BYTE* const base = ms->window.base;
    U32 const target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;

    if (idx != target)
        DEBUGLOG(7, "ZSTD_updateDUBT, from %u to %u (dictLimit:%u)",
                    idx, target, ms->window.dictLimit);
    assert(ip + 8 <= iend);   /* condition for ZSTD_hashPtr */
    (void)iend;

    assert(idx >= ms->window.dictLimit);   /* condition for valid base+idx */
    for ( ; idx < target ; idx++) {
        size_t const h  = ZSTD_hashPtr(base + idx, hashLog, mls);   /* assumption : ip + 8 <= iend */
        U32    const matchIndex = hashTable[h];

        U32*   const nextCandidatePtr = bt + 2*(idx&btMask);
        U32*   const sortMarkPtr  = nextCandidatePtr + 1;

        DEBUGLOG(8, "ZSTD_updateDUBT: insert %u", idx);
        hashTable[h] = idx;   /* Update Hash Table */
        *nextCandidatePtr = matchIndex;   /* update BT like a chain */
        *sortMarkPtr = ZSTD_DUBT_UNSORTED_MARK;
    }
    ms->nextToUpdate = target;
}


/** ZSTD_insertDUBT1() :
 *  sort one already inserted but unsorted position
 *  assumption : curr >= btlow == (curr - btmask)
 *  doesn't fail */
static
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
void ZSTD_insertDUBT1(const ZSTD_matchState_t* ms,
                 U32 curr, const BYTE* inputEnd,
                 U32 nbCompares, U32 btLow,
                 const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const bt = ms->chainTable;
    U32  const btLog  = cParams->chainLog - 1;
    U32  const btMask = (1 << btLog) - 1;
    size_t commonLengthSmaller=0, commonLengthLarger=0;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const ip = (curr>=dictLimit) ? base + curr : dictBase + curr;
    const BYTE* const iend = (curr>=dictLimit) ? inputEnd : dictBase + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* match;
    U32* smallerPtr = bt + 2*(curr&btMask);
    U32* largerPtr  = smallerPtr + 1;
    U32 matchIndex = *smallerPtr;   /* this candidate is unsorted : next sorted candidate is reached through *smallerPtr, while *largerPtr contains previous unsorted candidate (which is already saved and can be overwritten) */
    U32 dummy32;   /* to be nullified at the end */
    U32 const windowValid = ms->window.lowLimit;
    U32 const maxDistance = 1U << cParams->windowLog;
    U32 const windowLow = (curr - windowValid > maxDistance) ? curr - maxDistance : windowValid;


    DEBUGLOG(8, "ZSTD_insertDUBT1(%u) (dictLimit=%u, lowLimit=%u)",
                curr, dictLimit, windowLow);
    assert(curr >= btLow);
    assert(ip < iend);   /* condition for ZSTD_count */

    for (; nbCompares && (matchIndex > windowLow); --nbCompares) {
        U32* const nextPtr = bt + 2*(matchIndex & btMask);
        size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
        assert(matchIndex < curr);
        /* note : all candidates are now supposed sorted,
         * but it's still possible to have nextPtr[1] == ZSTD_DUBT_UNSORTED_MARK
         * when a real index has the same value as ZSTD_DUBT_UNSORTED_MARK */

        if ( (dictMode != ZSTD_extDict)
          || (matchIndex+matchLength >= dictLimit)  /* both in current segment*/
          || (curr < dictLimit) /* both in extDict */) {
            const BYTE* const mBase = ( (dictMode != ZSTD_extDict)
                                     || (matchIndex+matchLength >= dictLimit)) ?
                                        base : dictBase;
            assert( (matchIndex+matchLength >= dictLimit)   /* might be wrong if extDict is incorrectly set to 0 */
                 || (curr < dictLimit) );
            match = mBase + matchIndex;
            matchLength += ZSTD_count(ip+matchLength, match+matchLength, iend);
        } else {
            match = dictBase + matchIndex;
            matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
            if (matchIndex+matchLength >= dictLimit)
                match = base + matchIndex;   /* preparation for next read of match[matchLength] */
        }

        DEBUGLOG(8, "ZSTD_insertDUBT1: comparing %u with %u : found %u common bytes ",
                    curr, matchIndex, (U32)matchLength);

        if (ip+matchLength == iend) {   /* equal : no way to know if inf or sup */
            break;   /* drop , to guarantee consistency ; miss a bit of compression, but other solutions can corrupt tree */
        }

        if (match[matchLength] < ip[matchLength]) {  /* necessarily within buffer */
            /* match is smaller than current */
            *smallerPtr = matchIndex;             /* update smaller idx */
            commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
            if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            DEBUGLOG(8, "ZSTD_insertDUBT1: %u (>btLow=%u) is smaller : next => %u",
                        matchIndex, btLow, nextPtr[1]);
            smallerPtr = nextPtr+1;               /* new "candidate" => larger than match, which was smaller than target */
            matchIndex = nextPtr[1];              /* new matchIndex, larger than previous and closer to current */
        } else {
            /* match is larger than current */
            *largerPtr = matchIndex;
            commonLengthLarger = matchLength;
            if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop searching */
            DEBUGLOG(8, "ZSTD_insertDUBT1: %u (>btLow=%u) is larger => %u",
                        matchIndex, btLow, nextPtr[0]);
            largerPtr = nextPtr;
            matchIndex = nextPtr[0];
    }   }

    *smallerPtr = *largerPtr = 0;
}


static
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_DUBT_findBetterDictMatch (
        const ZSTD_matchState_t* ms,
        const BYTE* const ip, const BYTE* const iend,
        size_t* offsetPtr,
        size_t bestLength,
        U32 nbCompares,
        U32 const mls,
        const ZSTD_dictMode_e dictMode)
{
    const ZSTD_matchState_t * const dms = ms->dictMatchState;
    const ZSTD_compressionParameters* const dmsCParams = &dms->cParams;
    const U32 * const dictHashTable = dms->hashTable;
    U32         const hashLog = dmsCParams->hashLog;
    size_t      const h  = ZSTD_hashPtr(ip, hashLog, mls);
    U32               dictMatchIndex = dictHashTable[h];

    const BYTE* const base = ms->window.base;
    const BYTE* const prefixStart = base + ms->window.dictLimit;
    U32         const curr = (U32)(ip-base);
    const BYTE* const dictBase = dms->window.base;
    const BYTE* const dictEnd = dms->window.nextSrc;
    U32         const dictHighLimit = (U32)(dms->window.nextSrc - dms->window.base);
    U32         const dictLowLimit = dms->window.lowLimit;
    U32         const dictIndexDelta = ms->window.lowLimit - dictHighLimit;

    U32*        const dictBt = dms->chainTable;
    U32         const btLog  = dmsCParams->chainLog - 1;
    U32         const btMask = (1 << btLog) - 1;
    U32         const btLow = (btMask >= dictHighLimit - dictLowLimit) ? dictLowLimit : dictHighLimit - btMask;

    size_t commonLengthSmaller=0, commonLengthLarger=0;

    (void)dictMode;
    assert(dictMode == ZSTD_dictMatchState);

    for (; nbCompares && (dictMatchIndex > dictLowLimit); --nbCompares) {
        U32* const nextPtr = dictBt + 2*(dictMatchIndex & btMask);
        size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
        const BYTE* match = dictBase + dictMatchIndex;
        matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
        if (dictMatchIndex+matchLength >= dictHighLimit)
            match = base + dictMatchIndex + dictIndexDelta;   /* to prepare for next usage of match[matchLength] */

        if (matchLength > bestLength) {
            U32 matchIndex = dictMatchIndex + dictIndexDelta;
            if ( (4*(int)(matchLength-bestLength)) > (int)(ZSTD_highbit32(curr-matchIndex+1) - ZSTD_highbit32((U32)offsetPtr[0]+1)) ) {
                DEBUGLOG(9, "ZSTD_DUBT_findBetterDictMatch(%u) : found better match length %u -> %u and offsetCode %u -> %u (dictMatchIndex %u, matchIndex %u)",
                    curr, (U32)bestLength, (U32)matchLength, (U32)*offsetPtr, OFFSET_TO_OFFBASE(curr - matchIndex), dictMatchIndex, matchIndex);
                bestLength = matchLength, *offsetPtr = OFFSET_TO_OFFBASE(curr - matchIndex);
            }
            if (ip+matchLength == iend) {   /* reached end of input : ip[matchLength] is not valid, no way to know if it's larger or smaller than match */
                break;   /* drop, to guarantee consistency (miss a little bit of compression) */
            }
        }

        if (match[matchLength] < ip[matchLength]) {
            if (dictMatchIndex <= btLow) { break; }   /* beyond tree size, stop the search */
            commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
            dictMatchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
        } else {
            /* match is larger than current */
            if (dictMatchIndex <= btLow) { break; }   /* beyond tree size, stop the search */
            commonLengthLarger = matchLength;
            dictMatchIndex = nextPtr[0];
        }
    }

    if (bestLength >= MINMATCH) {
        U32 const mIndex = curr - (U32)OFFBASE_TO_OFFSET(*offsetPtr); (void)mIndex;
        DEBUGLOG(8, "ZSTD_DUBT_findBetterDictMatch(%u) : found match of length %u and offsetCode %u (pos %u)",
                    curr, (U32)bestLength, (U32)*offsetPtr, mIndex);
    }
    return bestLength;

}


static
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_DUBT_findBestMatch(ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iend,
                        size_t* offBasePtr,
                        U32 const mls,
                        const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32*   const hashTable = ms->hashTable;
    U32    const hashLog = cParams->hashLog;
    size_t const h  = ZSTD_hashPtr(ip, hashLog, mls);
    U32          matchIndex  = hashTable[h];

    const BYTE* const base = ms->window.base;
    U32    const curr = (U32)(ip-base);
    U32    const windowLow = ZSTD_getLowestMatchIndex(ms, curr, cParams->windowLog);

    U32*   const bt = ms->chainTable;
    U32    const btLog  = cParams->chainLog - 1;
    U32    const btMask = (1 << btLog) - 1;
    U32    const btLow = (btMask >= curr) ? 0 : curr - btMask;
    U32    const unsortLimit = MAX(btLow, windowLow);

    U32*         nextCandidate = bt + 2*(matchIndex&btMask);
    U32*         unsortedMark = bt + 2*(matchIndex&btMask) + 1;
    U32          nbCompares = 1U << cParams->searchLog;
    U32          nbCandidates = nbCompares;
    U32          previousCandidate = 0;

    DEBUGLOG(7, "ZSTD_DUBT_findBestMatch (%u) ", curr);
    assert(ip <= iend-8);   /* required for h calculation */
    assert(dictMode != ZSTD_dedicatedDictSearch);

    /* reach end of unsorted candidates list */
    while ( (matchIndex > unsortLimit)
         && (*unsortedMark == ZSTD_DUBT_UNSORTED_MARK)
         && (nbCandidates > 1) ) {
        DEBUGLOG(8, "ZSTD_DUBT_findBestMatch: candidate %u is unsorted",
                    matchIndex);
        *unsortedMark = previousCandidate;  /* the unsortedMark becomes a reversed chain, to move up back to original position */
        previousCandidate = matchIndex;
        matchIndex = *nextCandidate;
        nextCandidate = bt + 2*(matchIndex&btMask);
        unsortedMark = bt + 2*(matchIndex&btMask) + 1;
        nbCandidates --;
    }

    /* nullify last candidate if it's still unsorted
     * simplification, detrimental to compression ratio, beneficial for speed */
    if ( (matchIndex > unsortLimit)
      && (*unsortedMark==ZSTD_DUBT_UNSORTED_MARK) ) {
        DEBUGLOG(7, "ZSTD_DUBT_findBestMatch: nullify last unsorted candidate %u",
                    matchIndex);
        *nextCandidate = *unsortedMark = 0;
    }

    /* batch sort stacked candidates */
    matchIndex = previousCandidate;
    while (matchIndex) {  /* will end on matchIndex == 0 */
        U32* const nextCandidateIdxPtr = bt + 2*(matchIndex&btMask) + 1;
        U32 const nextCandidateIdx = *nextCandidateIdxPtr;
        ZSTD_insertDUBT1(ms, matchIndex, iend,
                         nbCandidates, unsortLimit, dictMode);
        matchIndex = nextCandidateIdx;
        nbCandidates++;
    }

    /* find longest match */
    {   size_t commonLengthSmaller = 0, commonLengthLarger = 0;
        const BYTE* const dictBase = ms->window.dictBase;
        const U32 dictLimit = ms->window.dictLimit;
        const BYTE* const dictEnd = dictBase + dictLimit;
        const BYTE* const prefixStart = base + dictLimit;
        U32* smallerPtr = bt + 2*(curr&btMask);
        U32* largerPtr  = bt + 2*(curr&btMask) + 1;
        U32 matchEndIdx = curr + 8 + 1;
        U32 dummy32;   /* to be nullified at the end */
        size_t bestLength = 0;

        matchIndex  = hashTable[h];
        hashTable[h] = curr;   /* Update Hash Table */

        for (; nbCompares && (matchIndex > windowLow); --nbCompares) {
            U32* const nextPtr = bt + 2*(matchIndex & btMask);
            size_t matchLength = MIN(commonLengthSmaller, commonLengthLarger);   /* guaranteed minimum nb of common bytes */
            const BYTE* match;

            if ((dictMode != ZSTD_extDict) || (matchIndex+matchLength >= dictLimit)) {
                match = base + matchIndex;
                matchLength += ZSTD_count(ip+matchLength, match+matchLength, iend);
            } else {
                match = dictBase + matchIndex;
                matchLength += ZSTD_count_2segments(ip+matchLength, match+matchLength, iend, dictEnd, prefixStart);
                if (matchIndex+matchLength >= dictLimit)
                    match = base + matchIndex;   /* to prepare for next usage of match[matchLength] */
            }

            if (matchLength > bestLength) {
                if (matchLength > matchEndIdx - matchIndex)
                    matchEndIdx = matchIndex + (U32)matchLength;
                if ( (4*(int)(matchLength-bestLength)) > (int)(ZSTD_highbit32(curr - matchIndex + 1) - ZSTD_highbit32((U32)*offBasePtr)) )
                    bestLength = matchLength, *offBasePtr = OFFSET_TO_OFFBASE(curr - matchIndex);
                if (ip+matchLength == iend) {   /* equal : no way to know if inf or sup */
                    if (dictMode == ZSTD_dictMatchState) {
                        nbCompares = 0; /* in addition to avoiding checking any
                                         * further in this loop, make sure we
                                         * skip checking in the dictionary. */
                    }
                    break;   /* drop, to guarantee consistency (miss a little bit of compression) */
                }
            }

            if (match[matchLength] < ip[matchLength]) {
                /* match is smaller than current */
                *smallerPtr = matchIndex;             /* update smaller idx */
                commonLengthSmaller = matchLength;    /* all smaller will now have at least this guaranteed common length */
                if (matchIndex <= btLow) { smallerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
                smallerPtr = nextPtr+1;               /* new "smaller" => larger of match */
                matchIndex = nextPtr[1];              /* new matchIndex larger than previous (closer to current) */
            } else {
                /* match is larger than current */
                *largerPtr = matchIndex;
                commonLengthLarger = matchLength;
                if (matchIndex <= btLow) { largerPtr=&dummy32; break; }   /* beyond tree size, stop the search */
                largerPtr = nextPtr;
                matchIndex = nextPtr[0];
        }   }

        *smallerPtr = *largerPtr = 0;

        assert(nbCompares <= (1U << ZSTD_SEARCHLOG_MAX)); /* Check we haven't underflowed. */
        if (dictMode == ZSTD_dictMatchState && nbCompares) {
            bestLength = ZSTD_DUBT_findBetterDictMatch(
                    ms, ip, iend,
                    offBasePtr, bestLength, nbCompares,
                    mls, dictMode);
        }

        assert(matchEndIdx > curr+8); /* ensure nextToUpdate is increased */
        ms->nextToUpdate = matchEndIdx - 8;   /* skip repetitive patterns */
        if (bestLength >= MINMATCH) {
            U32 const mIndex = curr - (U32)OFFBASE_TO_OFFSET(*offBasePtr); (void)mIndex;
            DEBUGLOG(8, "ZSTD_DUBT_findBestMatch(%u) : found match of length %u and offsetCode %u (pos %u)",
                        curr, (U32)bestLength, (U32)*offBasePtr, mIndex);
        }
        return bestLength;
    }
}


/** ZSTD_BtFindBestMatch() : Tree updater, providing best match */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_BtFindBestMatch( ZSTD_matchState_t* ms,
                const BYTE* const ip, const BYTE* const iLimit,
                      size_t* offBasePtr,
                const U32 mls /* template */,
                const ZSTD_dictMode_e dictMode)
{
    DEBUGLOG(7, "ZSTD_BtFindBestMatch");
    if (ip < ms->window.base + ms->nextToUpdate) return 0;   /* skipped area */
    ZSTD_updateDUBT(ms, ip, iLimit, mls);
    return ZSTD_DUBT_findBestMatch(ms, ip, iLimit, offBasePtr, mls, dictMode);
}

/***********************************
* Dedicated dict search
***********************************/

void ZSTD_dedicatedDictSearch_lazy_loadDictionary(ZSTD_matchState_t* ms, const BYTE* const ip)
{
    const BYTE* const base = ms->window.base;
    U32 const target = (U32)(ip - base);
    U32* const hashTable = ms->hashTable;
    U32* const chainTable = ms->chainTable;
    U32 const chainSize = 1 << ms->cParams.chainLog;
    U32 idx = ms->nextToUpdate;
    U32 const minChain = chainSize < target - idx ? target - chainSize : idx;
    U32 const bucketSize = 1 << ZSTD_LAZY_DDSS_BUCKET_LOG;
    U32 const cacheSize = bucketSize - 1;
    U32 const chainAttempts = (1 << ms->cParams.searchLog) - cacheSize;
    U32 const chainLimit = chainAttempts > 255 ? 255 : chainAttempts;

    /* We know the hashtable is oversized by a factor of `bucketSize`.
     * We are going to temporarily pretend `bucketSize == 1`, keeping only a
     * single entry. We will use the rest of the space to construct a temporary
     * chaintable.
     */
    U32 const hashLog = ms->cParams.hashLog - ZSTD_LAZY_DDSS_BUCKET_LOG;
    U32* const tmpHashTable = hashTable;
    U32* const tmpChainTable = hashTable + ((size_t)1 << hashLog);
    U32 const tmpChainSize = (U32)((1 << ZSTD_LAZY_DDSS_BUCKET_LOG) - 1) << hashLog;
    U32 const tmpMinChain = tmpChainSize < target ? target - tmpChainSize : idx;
    U32 hashIdx;

    assert(ms->cParams.chainLog <= 24);
    assert(ms->cParams.hashLog > ms->cParams.chainLog);
    assert(idx != 0);
    assert(tmpMinChain <= minChain);

    /* fill conventional hash table and conventional chain table */
    for ( ; idx < target; idx++) {
        U32 const h = (U32)ZSTD_hashPtr(base + idx, hashLog, ms->cParams.minMatch);
        if (idx >= tmpMinChain) {
            tmpChainTable[idx - tmpMinChain] = hashTable[h];
        }
        tmpHashTable[h] = idx;
    }

    /* sort chains into ddss chain table */
    {
        U32 chainPos = 0;
        for (hashIdx = 0; hashIdx < (1U << hashLog); hashIdx++) {
            U32 count;
            U32 countBeyondMinChain = 0;
            U32 i = tmpHashTable[hashIdx];
            for (count = 0; i >= tmpMinChain && count < cacheSize; count++) {
                /* skip through the chain to the first position that won't be
                 * in the hash cache bucket */
                if (i < minChain) {
                    countBeyondMinChain++;
                }
                i = tmpChainTable[i - tmpMinChain];
            }
            if (count == cacheSize) {
                for (count = 0; count < chainLimit;) {
                    if (i < minChain) {
                        if (!i || ++countBeyondMinChain > cacheSize) {
                            /* only allow pulling `cacheSize` number of entries
                             * into the cache or chainTable beyond `minChain`,
                             * to replace the entries pulled out of the
                             * chainTable into the cache. This lets us reach
                             * back further without increasing the total number
                             * of entries in the chainTable, guaranteeing the
                             * DDSS chain table will fit into the space
                             * allocated for the regular one. */
                            break;
                        }
                    }
                    chainTable[chainPos++] = i;
                    count++;
                    if (i < tmpMinChain) {
                        break;
                    }
                    i = tmpChainTable[i - tmpMinChain];
                }
            } else {
                count = 0;
            }
            if (count) {
                tmpHashTable[hashIdx] = ((chainPos - count) << 8) + count;
            } else {
                tmpHashTable[hashIdx] = 0;
            }
        }
        assert(chainPos <= chainSize); /* I believe this is guaranteed... */
    }

    /* move chain pointers into the last entry of each hash bucket */
    for (hashIdx = (1 << hashLog); hashIdx; ) {
        U32 const bucketIdx = --hashIdx << ZSTD_LAZY_DDSS_BUCKET_LOG;
        U32 const chainPackedPointer = tmpHashTable[hashIdx];
        U32 i;
        for (i = 0; i < cacheSize; i++) {
            hashTable[bucketIdx + i] = 0;
        }
        hashTable[bucketIdx + bucketSize - 1] = chainPackedPointer;
    }

    /* fill the buckets of the hash table */
    for (idx = ms->nextToUpdate; idx < target; idx++) {
        U32 const h = (U32)ZSTD_hashPtr(base + idx, hashLog, ms->cParams.minMatch)
                   << ZSTD_LAZY_DDSS_BUCKET_LOG;
        U32 i;
        /* Shift hash cache down 1. */
        for (i = cacheSize - 1; i; i--)
            hashTable[h + i] = hashTable[h + i - 1];
        hashTable[h] = idx;
    }

    ms->nextToUpdate = target;
}

/* Returns the longest match length found in the dedicated dict search structure.
 * If none are longer than the argument ml, then ml will be returned.
 */
FORCE_INLINE_TEMPLATE
size_t ZSTD_dedicatedDictSearch_lazy_search(size_t* offsetPtr, size_t ml, U32 nbAttempts,
                                            const ZSTD_matchState_t* const dms,
                                            const BYTE* const ip, const BYTE* const iLimit,
                                            const BYTE* const prefixStart, const U32 curr,
                                            const U32 dictLimit, const size_t ddsIdx) {
    const U32 ddsLowestIndex  = dms->window.dictLimit;
    const BYTE* const ddsBase = dms->window.base;
    const BYTE* const ddsEnd  = dms->window.nextSrc;
    const U32 ddsSize         = (U32)(ddsEnd - ddsBase);
    const U32 ddsIndexDelta   = dictLimit - ddsSize;
    const U32 bucketSize      = (1 << ZSTD_LAZY_DDSS_BUCKET_LOG);
    const U32 bucketLimit     = nbAttempts < bucketSize - 1 ? nbAttempts : bucketSize - 1;
    U32 ddsAttempt;
    U32 matchIndex;

    for (ddsAttempt = 0; ddsAttempt < bucketSize - 1; ddsAttempt++) {
        PREFETCH_L1(ddsBase + dms->hashTable[ddsIdx + ddsAttempt]);
    }

    {
        U32 const chainPackedPointer = dms->hashTable[ddsIdx + bucketSize - 1];
        U32 const chainIndex = chainPackedPointer >> 8;

        PREFETCH_L1(&dms->chainTable[chainIndex]);
    }

    for (ddsAttempt = 0; ddsAttempt < bucketLimit; ddsAttempt++) {
        size_t currentMl=0;
        const BYTE* match;
        matchIndex = dms->hashTable[ddsIdx + ddsAttempt];
        match = ddsBase + matchIndex;

        if (!matchIndex) {
            return ml;
        }

        /* guaranteed by table construction */
        (void)ddsLowestIndex;
        assert(matchIndex >= ddsLowestIndex);
        assert(match+4 <= ddsEnd);
        if (MEM_read32(match) == MEM_read32(ip)) {
            /* assumption : matchIndex <= dictLimit-4 (by table construction) */
            currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, ddsEnd, prefixStart) + 4;
        }

        /* save best solution */
        if (currentMl > ml) {
            ml = currentMl;
            *offsetPtr = OFFSET_TO_OFFBASE(curr - (matchIndex + ddsIndexDelta));
            if (ip+currentMl == iLimit) {
                /* best possible, avoids read overflow on next attempt */
                return ml;
            }
        }
    }

    {
        U32 const chainPackedPointer = dms->hashTable[ddsIdx + bucketSize - 1];
        U32 chainIndex = chainPackedPointer >> 8;
        U32 const chainLength = chainPackedPointer & 0xFF;
        U32 const chainAttempts = nbAttempts - ddsAttempt;
        U32 const chainLimit = chainAttempts > chainLength ? chainLength : chainAttempts;
        U32 chainAttempt;

        for (chainAttempt = 0 ; chainAttempt < chainLimit; chainAttempt++) {
            PREFETCH_L1(ddsBase + dms->chainTable[chainIndex + chainAttempt]);
        }

        for (chainAttempt = 0 ; chainAttempt < chainLimit; chainAttempt++, chainIndex++) {
            size_t currentMl=0;
            const BYTE* match;
            matchIndex = dms->chainTable[chainIndex];
            match = ddsBase + matchIndex;

            /* guaranteed by table construction */
            assert(matchIndex >= ddsLowestIndex);
            assert(match+4 <= ddsEnd);
            if (MEM_read32(match) == MEM_read32(ip)) {
                /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, ddsEnd, prefixStart) + 4;
            }

            /* save best solution */
            if (currentMl > ml) {
                ml = currentMl;
                *offsetPtr = OFFSET_TO_OFFBASE(curr - (matchIndex + ddsIndexDelta));
                if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
            }
        }
    }
    return ml;
}


/* *********************************
*  Hash Chain
***********************************/
#define NEXT_IN_CHAIN(d, mask)   chainTable[(d) & (mask)]

/* Update chains up to ip (excluded)
   Assumption : always within prefix (i.e. not within extDict) */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
U32 ZSTD_insertAndFindFirstIndex_internal(
                        ZSTD_matchState_t* ms,
                        const ZSTD_compressionParameters* const cParams,
                        const BYTE* ip, U32 const mls, U32 const lazySkipping)
{
    U32* const hashTable  = ms->hashTable;
    const U32 hashLog = cParams->hashLog;
    U32* const chainTable = ms->chainTable;
    const U32 chainMask = (1 << cParams->chainLog) - 1;
    const BYTE* const base = ms->window.base;
    const U32 target = (U32)(ip - base);
    U32 idx = ms->nextToUpdate;

    while(idx < target) { /* catch up */
        size_t const h = ZSTD_hashPtr(base+idx, hashLog, mls);
        NEXT_IN_CHAIN(idx, chainMask) = hashTable[h];
        hashTable[h] = idx;
        idx++;
        /* Stop inserting every position when in the lazy skipping mode. */
        if (lazySkipping)
            break;
    }

    ms->nextToUpdate = target;
    return hashTable[ZSTD_hashPtr(ip, hashLog, mls)];
}

U32 ZSTD_insertAndFindFirstIndex(ZSTD_matchState_t* ms, const BYTE* ip) {
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    return ZSTD_insertAndFindFirstIndex_internal(ms, cParams, ip, ms->cParams.minMatch, /* lazySkipping*/ 0);
}

/* inlining is important to hardwire a hot branch (template emulation) */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_HcFindBestMatch(
                        ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iLimit,
                        size_t* offsetPtr,
                        const U32 mls, const ZSTD_dictMode_e dictMode)
{
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    U32* const chainTable = ms->chainTable;
    const U32 chainSize = (1 << cParams->chainLog);
    const U32 chainMask = chainSize-1;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const U32 curr = (U32)(ip-base);
    const U32 maxDistance = 1U << cParams->windowLog;
    const U32 lowestValid = ms->window.lowLimit;
    const U32 withinMaxDistance = (curr - lowestValid > maxDistance) ? curr - maxDistance : lowestValid;
    const U32 isDictionary = (ms->loadedDictEnd != 0);
    const U32 lowLimit = isDictionary ? lowestValid : withinMaxDistance;
    const U32 minChain = curr > chainSize ? curr - chainSize : 0;
    U32 nbAttempts = 1U << cParams->searchLog;
    size_t ml=4-1;

    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const U32 ddsHashLog = dictMode == ZSTD_dedicatedDictSearch
                         ? dms->cParams.hashLog - ZSTD_LAZY_DDSS_BUCKET_LOG : 0;
    const size_t ddsIdx = dictMode == ZSTD_dedicatedDictSearch
                        ? ZSTD_hashPtr(ip, ddsHashLog, mls) << ZSTD_LAZY_DDSS_BUCKET_LOG : 0;

    U32 matchIndex;

    if (dictMode == ZSTD_dedicatedDictSearch) {
        const U32* entry = &dms->hashTable[ddsIdx];
        PREFETCH_L1(entry);
    }

    /* HC4 match finder */
    matchIndex = ZSTD_insertAndFindFirstIndex_internal(ms, cParams, ip, mls, ms->lazySkipping);

    for ( ; (matchIndex>=lowLimit) & (nbAttempts>0) ; nbAttempts--) {
        size_t currentMl=0;
        if ((dictMode != ZSTD_extDict) || matchIndex >= dictLimit) {
            const BYTE* const match = base + matchIndex;
            assert(matchIndex >= dictLimit);   /* ensures this is true if dictMode != ZSTD_extDict */
            /* read 4B starting from (match + ml + 1 - sizeof(U32)) */
            if (MEM_read32(match + ml - 3) == MEM_read32(ip + ml - 3))   /* potentially better */
                currentMl = ZSTD_count(ip, match, iLimit);
        } else {
            const BYTE* const match = dictBase + matchIndex;
            assert(match+4 <= dictEnd);
            if (MEM_read32(match) == MEM_read32(ip))   /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dictEnd, prefixStart) + 4;
        }

        /* save best solution */
        if (currentMl > ml) {
            ml = currentMl;
            *offsetPtr = OFFSET_TO_OFFBASE(curr - matchIndex);
            if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
        }

        if (matchIndex <= minChain) break;
        matchIndex = NEXT_IN_CHAIN(matchIndex, chainMask);
    }

    assert(nbAttempts <= (1U << ZSTD_SEARCHLOG_MAX)); /* Check we haven't underflowed. */
    if (dictMode == ZSTD_dedicatedDictSearch) {
        ml = ZSTD_dedicatedDictSearch_lazy_search(offsetPtr, ml, nbAttempts, dms,
                                                  ip, iLimit, prefixStart, curr, dictLimit, ddsIdx);
    } else if (dictMode == ZSTD_dictMatchState) {
        const U32* const dmsChainTable = dms->chainTable;
        const U32 dmsChainSize         = (1 << dms->cParams.chainLog);
        const U32 dmsChainMask         = dmsChainSize - 1;
        const U32 dmsLowestIndex       = dms->window.dictLimit;
        const BYTE* const dmsBase      = dms->window.base;
        const BYTE* const dmsEnd       = dms->window.nextSrc;
        const U32 dmsSize              = (U32)(dmsEnd - dmsBase);
        const U32 dmsIndexDelta        = dictLimit - dmsSize;
        const U32 dmsMinChain = dmsSize > dmsChainSize ? dmsSize - dmsChainSize : 0;

        matchIndex = dms->hashTable[ZSTD_hashPtr(ip, dms->cParams.hashLog, mls)];

        for ( ; (matchIndex>=dmsLowestIndex) & (nbAttempts>0) ; nbAttempts--) {
            size_t currentMl=0;
            const BYTE* const match = dmsBase + matchIndex;
            assert(match+4 <= dmsEnd);
            if (MEM_read32(match) == MEM_read32(ip))   /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dmsEnd, prefixStart) + 4;

            /* save best solution */
            if (currentMl > ml) {
                ml = currentMl;
                assert(curr > matchIndex + dmsIndexDelta);
                *offsetPtr = OFFSET_TO_OFFBASE(curr - (matchIndex + dmsIndexDelta));
                if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
            }

            if (matchIndex <= dmsMinChain) break;

            matchIndex = dmsChainTable[matchIndex & dmsChainMask];
        }
    }

    return ml;
}

/* *********************************
* (SIMD) Row-based matchfinder
***********************************/
/* Constants for row-based hash */
#define ZSTD_ROW_HASH_TAG_MASK ((1u << ZSTD_ROW_HASH_TAG_BITS) - 1)
#define ZSTD_ROW_HASH_MAX_ENTRIES 64    /* absolute maximum number of entries per row, for all configurations */

#define ZSTD_ROW_HASH_CACHE_MASK (ZSTD_ROW_HASH_CACHE_SIZE - 1)

typedef U64 ZSTD_VecMask;   /* Clarifies when we are interacting with a U64 representing a mask of matches */

/* ZSTD_VecMask_next():
 * Starting from the LSB, returns the idx of the next non-zero bit.
 * Basically counting the nb of trailing zeroes.
 */
MEM_STATIC U32 ZSTD_VecMask_next(ZSTD_VecMask val) {
    return ZSTD_countTrailingZeros64(val);
}

/* ZSTD_row_nextIndex():
 * Returns the next index to insert at within a tagTable row, and updates the "head"
 * value to reflect the update. Essentially cycles backwards from [1, {entries per row})
 */
FORCE_INLINE_TEMPLATE U32 ZSTD_row_nextIndex(BYTE* const tagRow, U32 const rowMask) {
    U32 next = (*tagRow-1) & rowMask;
    next += (next == 0) ? rowMask : 0; /* skip first position */
    *tagRow = (BYTE)next;
    return next;
}

/* ZSTD_isAligned():
 * Checks that a pointer is aligned to "align" bytes which must be a power of 2.
 */
MEM_STATIC int ZSTD_isAligned(void const* ptr, size_t align) {
    assert((align & (align - 1)) == 0);
    return (((size_t)ptr) & (align - 1)) == 0;
}

/* ZSTD_row_prefetch():
 * Performs prefetching for the hashTable and tagTable at a given row.
 */
FORCE_INLINE_TEMPLATE void ZSTD_row_prefetch(U32 const* hashTable, BYTE const* tagTable, U32 const relRow, U32 const rowLog) {
    PREFETCH_L1(hashTable + relRow);
    if (rowLog >= 5) {
        PREFETCH_L1(hashTable + relRow + 16);
        /* Note: prefetching more of the hash table does not appear to be beneficial for 128-entry rows */
    }
    PREFETCH_L1(tagTable + relRow);
    if (rowLog == 6) {
        PREFETCH_L1(tagTable + relRow + 32);
    }
    assert(rowLog == 4 || rowLog == 5 || rowLog == 6);
    assert(ZSTD_isAligned(hashTable + relRow, 64));                 /* prefetched hash row always 64-byte aligned */
    assert(ZSTD_isAligned(tagTable + relRow, (size_t)1 << rowLog)); /* prefetched tagRow sits on correct multiple of bytes (32,64,128) */
}

/* ZSTD_row_fillHashCache():
 * Fill up the hash cache starting at idx, prefetching up to ZSTD_ROW_HASH_CACHE_SIZE entries,
 * but not beyond iLimit.
 */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
void ZSTD_row_fillHashCache(ZSTD_matchState_t* ms, const BYTE* base,
                                   U32 const rowLog, U32 const mls,
                                   U32 idx, const BYTE* const iLimit)
{
    U32 const* const hashTable = ms->hashTable;
    BYTE const* const tagTable = ms->tagTable;
    U32 const hashLog = ms->rowHashLog;
    U32 const maxElemsToPrefetch = (base + idx) > iLimit ? 0 : (U32)(iLimit - (base + idx) + 1);
    U32 const lim = idx + MIN(ZSTD_ROW_HASH_CACHE_SIZE, maxElemsToPrefetch);

    for (; idx < lim; ++idx) {
        U32 const hash = (U32)ZSTD_hashPtrSalted(base + idx, hashLog + ZSTD_ROW_HASH_TAG_BITS, mls, ms->hashSalt);
        U32 const row = (hash >> ZSTD_ROW_HASH_TAG_BITS) << rowLog;
        ZSTD_row_prefetch(hashTable, tagTable, row, rowLog);
        ms->hashCache[idx & ZSTD_ROW_HASH_CACHE_MASK] = hash;
    }

    DEBUGLOG(6, "ZSTD_row_fillHashCache(): [%u %u %u %u %u %u %u %u]", ms->hashCache[0], ms->hashCache[1],
                                                     ms->hashCache[2], ms->hashCache[3], ms->hashCache[4],
                                                     ms->hashCache[5], ms->hashCache[6], ms->hashCache[7]);
}

/* ZSTD_row_nextCachedHash():
 * Returns the hash of base + idx, and replaces the hash in the hash cache with the byte at
 * base + idx + ZSTD_ROW_HASH_CACHE_SIZE. Also prefetches the appropriate rows from hashTable and tagTable.
 */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
U32 ZSTD_row_nextCachedHash(U32* cache, U32 const* hashTable,
                                                  BYTE const* tagTable, BYTE const* base,
                                                  U32 idx, U32 const hashLog,
                                                  U32 const rowLog, U32 const mls,
                                                  U64 const hashSalt)
{
    U32 const newHash = (U32)ZSTD_hashPtrSalted(base+idx+ZSTD_ROW_HASH_CACHE_SIZE, hashLog + ZSTD_ROW_HASH_TAG_BITS, mls, hashSalt);
    U32 const row = (newHash >> ZSTD_ROW_HASH_TAG_BITS) << rowLog;
    ZSTD_row_prefetch(hashTable, tagTable, row, rowLog);
    {   U32 const hash = cache[idx & ZSTD_ROW_HASH_CACHE_MASK];
        cache[idx & ZSTD_ROW_HASH_CACHE_MASK] = newHash;
        return hash;
    }
}

/* ZSTD_row_update_internalImpl():
 * Updates the hash table with positions starting from updateStartIdx until updateEndIdx.
 */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
void ZSTD_row_update_internalImpl(ZSTD_matchState_t* ms,
                                  U32 updateStartIdx, U32 const updateEndIdx,
                                  U32 const mls, U32 const rowLog,
                                  U32 const rowMask, U32 const useCache)
{
    U32* const hashTable = ms->hashTable;
    BYTE* const tagTable = ms->tagTable;
    U32 const hashLog = ms->rowHashLog;
    const BYTE* const base = ms->window.base;

    DEBUGLOG(6, "ZSTD_row_update_internalImpl(): updateStartIdx=%u, updateEndIdx=%u", updateStartIdx, updateEndIdx);
    for (; updateStartIdx < updateEndIdx; ++updateStartIdx) {
        U32 const hash = useCache ? ZSTD_row_nextCachedHash(ms->hashCache, hashTable, tagTable, base, updateStartIdx, hashLog, rowLog, mls, ms->hashSalt)
                                  : (U32)ZSTD_hashPtrSalted(base + updateStartIdx, hashLog + ZSTD_ROW_HASH_TAG_BITS, mls, ms->hashSalt);
        U32 const relRow = (hash >> ZSTD_ROW_HASH_TAG_BITS) << rowLog;
        U32* const row = hashTable + relRow;
        BYTE* tagRow = tagTable + relRow;
        U32 const pos = ZSTD_row_nextIndex(tagRow, rowMask);

        assert(hash == ZSTD_hashPtrSalted(base + updateStartIdx, hashLog + ZSTD_ROW_HASH_TAG_BITS, mls, ms->hashSalt));
        tagRow[pos] = hash & ZSTD_ROW_HASH_TAG_MASK;
        row[pos] = updateStartIdx;
    }
}

/* ZSTD_row_update_internal():
 * Inserts the byte at ip into the appropriate position in the hash table, and updates ms->nextToUpdate.
 * Skips sections of long matches as is necessary.
 */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
void ZSTD_row_update_internal(ZSTD_matchState_t* ms, const BYTE* ip,
                              U32 const mls, U32 const rowLog,
                              U32 const rowMask, U32 const useCache)
{
    U32 idx = ms->nextToUpdate;
    const BYTE* const base = ms->window.base;
    const U32 target = (U32)(ip - base);
    const U32 kSkipThreshold = 384;
    const U32 kMaxMatchStartPositionsToUpdate = 96;
    const U32 kMaxMatchEndPositionsToUpdate = 32;

    if (useCache) {
        /* Only skip positions when using hash cache, i.e.
         * if we are loading a dict, don't skip anything.
         * If we decide to skip, then we only update a set number
         * of positions at the beginning and end of the match.
         */
        if (UNLIKELY(target - idx > kSkipThreshold)) {
            U32 const bound = idx + kMaxMatchStartPositionsToUpdate;
            ZSTD_row_update_internalImpl(ms, idx, bound, mls, rowLog, rowMask, useCache);
            idx = target - kMaxMatchEndPositionsToUpdate;
            ZSTD_row_fillHashCache(ms, base, rowLog, mls, idx, ip+1);
        }
    }
    assert(target >= idx);
    ZSTD_row_update_internalImpl(ms, idx, target, mls, rowLog, rowMask, useCache);
    ms->nextToUpdate = target;
}

/* ZSTD_row_update():
 * External wrapper for ZSTD_row_update_internal(). Used for filling the hashtable during dictionary
 * processing.
 */
void ZSTD_row_update(ZSTD_matchState_t* const ms, const BYTE* ip) {
    const U32 rowLog = BOUNDED(4, ms->cParams.searchLog, 6);
    const U32 rowMask = (1u << rowLog) - 1;
    const U32 mls = MIN(ms->cParams.minMatch, 6 /* mls caps out at 6 */);

    DEBUGLOG(5, "ZSTD_row_update(), rowLog=%u", rowLog);
    ZSTD_row_update_internal(ms, ip, mls, rowLog, rowMask, 0 /* don't use cache */);
}

/* Returns the mask width of bits group of which will be set to 1. Given not all
 * architectures have easy movemask instruction, this helps to iterate over
 * groups of bits easier and faster.
 */
FORCE_INLINE_TEMPLATE U32
ZSTD_row_matchMaskGroupWidth(const U32 rowEntries)
{
    assert((rowEntries == 16) || (rowEntries == 32) || rowEntries == 64);
    assert(rowEntries <= ZSTD_ROW_HASH_MAX_ENTRIES);
    (void)rowEntries;
#if defined(ZSTD_ARCH_ARM_NEON)
    /* NEON path only works for little endian */
    if (!MEM_isLittleEndian()) {
        return 1;
    }
    if (rowEntries == 16) {
        return 4;
    }
    if (rowEntries == 32) {
        return 2;
    }
    if (rowEntries == 64) {
        return 1;
    }
#endif
    return 1;
}

#if defined(ZSTD_ARCH_X86_SSE2)
FORCE_INLINE_TEMPLATE ZSTD_VecMask
ZSTD_row_getSSEMask(int nbChunks, const BYTE* const src, const BYTE tag, const U32 head)
{
    const __m128i comparisonMask = _mm_set1_epi8((char)tag);
    int matches[4] = {0};
    int i;
    assert(nbChunks == 1 || nbChunks == 2 || nbChunks == 4);
    for (i=0; i<nbChunks; i++) {
        const __m128i chunk = _mm_loadu_si128((const __m128i*)(const void*)(src + 16*i));
        const __m128i equalMask = _mm_cmpeq_epi8(chunk, comparisonMask);
        matches[i] = _mm_movemask_epi8(equalMask);
    }
    if (nbChunks == 1) return ZSTD_rotateRight_U16((U16)matches[0], head);
    if (nbChunks == 2) return ZSTD_rotateRight_U32((U32)matches[1] << 16 | (U32)matches[0], head);
    assert(nbChunks == 4);
    return ZSTD_rotateRight_U64((U64)matches[3] << 48 | (U64)matches[2] << 32 | (U64)matches[1] << 16 | (U64)matches[0], head);
}
#endif

#if defined(ZSTD_ARCH_ARM_NEON)
FORCE_INLINE_TEMPLATE ZSTD_VecMask
ZSTD_row_getNEONMask(const U32 rowEntries, const BYTE* const src, const BYTE tag, const U32 headGrouped)
{
    assert((rowEntries == 16) || (rowEntries == 32) || rowEntries == 64);
    if (rowEntries == 16) {
        /* vshrn_n_u16 shifts by 4 every u16 and narrows to 8 lower bits.
         * After that groups of 4 bits represent the equalMask. We lower
         * all bits except the highest in these groups by doing AND with
         * 0x88 = 0b10001000.
         */
        const uint8x16_t chunk = vld1q_u8(src);
        const uint16x8_t equalMask = vreinterpretq_u16_u8(vceqq_u8(chunk, vdupq_n_u8(tag)));
        const uint8x8_t res = vshrn_n_u16(equalMask, 4);
        const U64 matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
        return ZSTD_rotateRight_U64(matches, headGrouped) & 0x8888888888888888ull;
    } else if (rowEntries == 32) {
        /* Same idea as with rowEntries == 16 but doing AND with
         * 0x55 = 0b01010101.
         */
        const uint16x8x2_t chunk = vld2q_u16((const uint16_t*)(const void*)src);
        const uint8x16_t chunk0 = vreinterpretq_u8_u16(chunk.val[0]);
        const uint8x16_t chunk1 = vreinterpretq_u8_u16(chunk.val[1]);
        const uint8x16_t dup = vdupq_n_u8(tag);
        const uint8x8_t t0 = vshrn_n_u16(vreinterpretq_u16_u8(vceqq_u8(chunk0, dup)), 6);
        const uint8x8_t t1 = vshrn_n_u16(vreinterpretq_u16_u8(vceqq_u8(chunk1, dup)), 6);
        const uint8x8_t res = vsli_n_u8(t0, t1, 4);
        const U64 matches = vget_lane_u64(vreinterpret_u64_u8(res), 0) ;
        return ZSTD_rotateRight_U64(matches, headGrouped) & 0x5555555555555555ull;
    } else { /* rowEntries == 64 */
        const uint8x16x4_t chunk = vld4q_u8(src);
        const uint8x16_t dup = vdupq_n_u8(tag);
        const uint8x16_t cmp0 = vceqq_u8(chunk.val[0], dup);
        const uint8x16_t cmp1 = vceqq_u8(chunk.val[1], dup);
        const uint8x16_t cmp2 = vceqq_u8(chunk.val[2], dup);
        const uint8x16_t cmp3 = vceqq_u8(chunk.val[3], dup);

        const uint8x16_t t0 = vsriq_n_u8(cmp1, cmp0, 1);
        const uint8x16_t t1 = vsriq_n_u8(cmp3, cmp2, 1);
        const uint8x16_t t2 = vsriq_n_u8(t1, t0, 2);
        const uint8x16_t t3 = vsriq_n_u8(t2, t2, 4);
        const uint8x8_t t4 = vshrn_n_u16(vreinterpretq_u16_u8(t3), 4);
        const U64 matches = vget_lane_u64(vreinterpret_u64_u8(t4), 0);
        return ZSTD_rotateRight_U64(matches, headGrouped);
    }
}
#endif

/* Returns a ZSTD_VecMask (U64) that has the nth group (determined by
 * ZSTD_row_matchMaskGroupWidth) of bits set to 1 if the newly-computed "tag"
 * matches the hash at the nth position in a row of the tagTable.
 * Each row is a circular buffer beginning at the value of "headGrouped". So we
 * must rotate the "matches" bitfield to match up with the actual layout of the
 * entries within the hashTable */
FORCE_INLINE_TEMPLATE ZSTD_VecMask
ZSTD_row_getMatchMask(const BYTE* const tagRow, const BYTE tag, const U32 headGrouped, const U32 rowEntries)
{
    const BYTE* const src = tagRow;
    assert((rowEntries == 16) || (rowEntries == 32) || rowEntries == 64);
    assert(rowEntries <= ZSTD_ROW_HASH_MAX_ENTRIES);
    assert(ZSTD_row_matchMaskGroupWidth(rowEntries) * rowEntries <= sizeof(ZSTD_VecMask) * 8);

#if defined(ZSTD_ARCH_X86_SSE2)

    return ZSTD_row_getSSEMask(rowEntries / 16, src, tag, headGrouped);

#else /* SW or NEON-LE */

# if defined(ZSTD_ARCH_ARM_NEON)
  /* This NEON path only works for little endian - otherwise use SWAR below */
    if (MEM_isLittleEndian()) {
        return ZSTD_row_getNEONMask(rowEntries, src, tag, headGrouped);
    }
# endif /* ZSTD_ARCH_ARM_NEON */
    /* SWAR */
    {   const int chunkSize = sizeof(size_t);
        const size_t shiftAmount = ((chunkSize * 8) - chunkSize);
        const size_t xFF = ~((size_t)0);
        const size_t x01 = xFF / 0xFF;
        const size_t x80 = x01 << 7;
        const size_t splatChar = tag * x01;
        ZSTD_VecMask matches = 0;
        int i = rowEntries - chunkSize;
        assert((sizeof(size_t) == 4) || (sizeof(size_t) == 8));
        if (MEM_isLittleEndian()) { /* runtime check so have two loops */
            const size_t extractMagic = (xFF / 0x7F) >> chunkSize;
            do {
                size_t chunk = MEM_readST(&src[i]);
                chunk ^= splatChar;
                chunk = (((chunk | x80) - x01) | chunk) & x80;
                matches <<= chunkSize;
                matches |= (chunk * extractMagic) >> shiftAmount;
                i -= chunkSize;
            } while (i >= 0);
        } else { /* big endian: reverse bits during extraction */
            const size_t msb = xFF ^ (xFF >> 1);
            const size_t extractMagic = (msb / 0x1FF) | msb;
            do {
                size_t chunk = MEM_readST(&src[i]);
                chunk ^= splatChar;
                chunk = (((chunk | x80) - x01) | chunk) & x80;
                matches <<= chunkSize;
                matches |= ((chunk >> 7) * extractMagic) >> shiftAmount;
                i -= chunkSize;
            } while (i >= 0);
        }
        matches = ~matches;
        if (rowEntries == 16) {
            return ZSTD_rotateRight_U16((U16)matches, headGrouped);
        } else if (rowEntries == 32) {
            return ZSTD_rotateRight_U32((U32)matches, headGrouped);
        } else {
            return ZSTD_rotateRight_U64((U64)matches, headGrouped);
        }
    }
#endif
}

/* The high-level approach of the SIMD row based match finder is as follows:
 * - Figure out where to insert the new entry:
 *      - Generate a hash from a byte along with an additional 1-byte "short hash". The additional byte is our "tag"
 *      - The hashTable is effectively split into groups or "rows" of 16 or 32 entries of U32, and the hash determines
 *        which row to insert into.
 *      - Determine the correct position within the row to insert the entry into. Each row of 16 or 32 can
 *        be considered as a circular buffer with a "head" index that resides in the tagTable.
 *      - Also insert the "tag" into the equivalent row and position in the tagTable.
 *          - Note: The tagTable has 17 or 33 1-byte entries per row, due to 16 or 32 tags, and 1 "head" entry.
 *                  The 17 or 33 entry rows are spaced out to occur every 32 or 64 bytes, respectively,
 *                  for alignment/performance reasons, leaving some bytes unused.
 * - Use SIMD to efficiently compare the tags in the tagTable to the 1-byte "short hash" and
 *   generate a bitfield that we can cycle through to check the collisions in the hash table.
 * - Pick the longest match.
 */
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_RowFindBestMatch(
                        ZSTD_matchState_t* ms,
                        const BYTE* const ip, const BYTE* const iLimit,
                        size_t* offsetPtr,
                        const U32 mls, const ZSTD_dictMode_e dictMode,
                        const U32 rowLog)
{
    U32* const hashTable = ms->hashTable;
    BYTE* const tagTable = ms->tagTable;
    U32* const hashCache = ms->hashCache;
    const U32 hashLog = ms->rowHashLog;
    const ZSTD_compressionParameters* const cParams = &ms->cParams;
    const BYTE* const base = ms->window.base;
    const BYTE* const dictBase = ms->window.dictBase;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictEnd = dictBase + dictLimit;
    const U32 curr = (U32)(ip-base);
    const U32 maxDistance = 1U << cParams->windowLog;
    const U32 lowestValid = ms->window.lowLimit;
    const U32 withinMaxDistance = (curr - lowestValid > maxDistance) ? curr - maxDistance : lowestValid;
    const U32 isDictionary = (ms->loadedDictEnd != 0);
    const U32 lowLimit = isDictionary ? lowestValid : withinMaxDistance;
    const U32 rowEntries = (1U << rowLog);
    const U32 rowMask = rowEntries - 1;
    const U32 cappedSearchLog = MIN(cParams->searchLog, rowLog); /* nb of searches is capped at nb entries per row */
    const U32 groupWidth = ZSTD_row_matchMaskGroupWidth(rowEntries);
    const U64 hashSalt = ms->hashSalt;
    U32 nbAttempts = 1U << cappedSearchLog;
    size_t ml=4-1;
    U32 hash;

    /* DMS/DDS variables that may be referenced laster */
    const ZSTD_matchState_t* const dms = ms->dictMatchState;

    /* Initialize the following variables to satisfy static analyzer */
    size_t ddsIdx = 0;
    U32 ddsExtraAttempts = 0; /* cctx hash tables are limited in searches, but allow extra searches into DDS */
    U32 dmsTag = 0;
    U32* dmsRow = NULL;
    BYTE* dmsTagRow = NULL;

    if (dictMode == ZSTD_dedicatedDictSearch) {
        const U32 ddsHashLog = dms->cParams.hashLog - ZSTD_LAZY_DDSS_BUCKET_LOG;
        {   /* Prefetch DDS hashtable entry */
            ddsIdx = ZSTD_hashPtr(ip, ddsHashLog, mls) << ZSTD_LAZY_DDSS_BUCKET_LOG;
            PREFETCH_L1(&dms->hashTable[ddsIdx]);
        }
        ddsExtraAttempts = cParams->searchLog > rowLog ? 1U << (cParams->searchLog - rowLog) : 0;
    }

    if (dictMode == ZSTD_dictMatchState) {
        /* Prefetch DMS rows */
        U32* const dmsHashTable = dms->hashTable;
        BYTE* const dmsTagTable = dms->tagTable;
        U32 const dmsHash = (U32)ZSTD_hashPtr(ip, dms->rowHashLog + ZSTD_ROW_HASH_TAG_BITS, mls);
        U32 const dmsRelRow = (dmsHash >> ZSTD_ROW_HASH_TAG_BITS) << rowLog;
        dmsTag = dmsHash & ZSTD_ROW_HASH_TAG_MASK;
        dmsTagRow = (BYTE*)(dmsTagTable + dmsRelRow);
        dmsRow = dmsHashTable + dmsRelRow;
        ZSTD_row_prefetch(dmsHashTable, dmsTagTable, dmsRelRow, rowLog);
    }

    /* Update the hashTable and tagTable up to (but not including) ip */
    if (!ms->lazySkipping) {
        ZSTD_row_update_internal(ms, ip, mls, rowLog, rowMask, 1 /* useCache */);
        hash = ZSTD_row_nextCachedHash(hashCache, hashTable, tagTable, base, curr, hashLog, rowLog, mls, hashSalt);
    } else {
        /* Stop inserting every position when in the lazy skipping mode.
         * The hash cache is also not kept up to date in this mode.
         */
        hash = (U32)ZSTD_hashPtrSalted(ip, hashLog + ZSTD_ROW_HASH_TAG_BITS, mls, hashSalt);
        ms->nextToUpdate = curr;
    }
    ms->hashSaltEntropy += hash; /* collect salt entropy */

    {   /* Get the hash for ip, compute the appropriate row */
        U32 const relRow = (hash >> ZSTD_ROW_HASH_TAG_BITS) << rowLog;
        U32 const tag = hash & ZSTD_ROW_HASH_TAG_MASK;
        U32* const row = hashTable + relRow;
        BYTE* tagRow = (BYTE*)(tagTable + relRow);
        U32 const headGrouped = (*tagRow & rowMask) * groupWidth;
        U32 matchBuffer[ZSTD_ROW_HASH_MAX_ENTRIES];
        size_t numMatches = 0;
        size_t currMatch = 0;
        ZSTD_VecMask matches = ZSTD_row_getMatchMask(tagRow, (BYTE)tag, headGrouped, rowEntries);

        /* Cycle through the matches and prefetch */
        for (; (matches > 0) && (nbAttempts > 0); matches &= (matches - 1)) {
            U32 const matchPos = ((headGrouped + ZSTD_VecMask_next(matches)) / groupWidth) & rowMask;
            U32 const matchIndex = row[matchPos];
            if(matchPos == 0) continue;
            assert(numMatches < rowEntries);
            if (matchIndex < lowLimit)
                break;
            if ((dictMode != ZSTD_extDict) || matchIndex >= dictLimit) {
                PREFETCH_L1(base + matchIndex);
            } else {
                PREFETCH_L1(dictBase + matchIndex);
            }
            matchBuffer[numMatches++] = matchIndex;
            --nbAttempts;
        }

        /* Speed opt: insert current byte into hashtable too. This allows us to avoid one iteration of the loop
           in ZSTD_row_update_internal() at the next search. */
        {
            U32 const pos = ZSTD_row_nextIndex(tagRow, rowMask);
            tagRow[pos] = (BYTE)tag;
            row[pos] = ms->nextToUpdate++;
        }

        /* Return the longest match */
        for (; currMatch < numMatches; ++currMatch) {
            U32 const matchIndex = matchBuffer[currMatch];
            size_t currentMl=0;
            assert(matchIndex < curr);
            assert(matchIndex >= lowLimit);

            if ((dictMode != ZSTD_extDict) || matchIndex >= dictLimit) {
                const BYTE* const match = base + matchIndex;
                assert(matchIndex >= dictLimit);   /* ensures this is true if dictMode != ZSTD_extDict */
                /* read 4B starting from (match + ml + 1 - sizeof(U32)) */
                if (MEM_read32(match + ml - 3) == MEM_read32(ip + ml - 3))   /* potentially better */
                    currentMl = ZSTD_count(ip, match, iLimit);
            } else {
                const BYTE* const match = dictBase + matchIndex;
                assert(match+4 <= dictEnd);
                if (MEM_read32(match) == MEM_read32(ip))   /* assumption : matchIndex <= dictLimit-4 (by table construction) */
                    currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dictEnd, prefixStart) + 4;
            }

            /* Save best solution */
            if (currentMl > ml) {
                ml = currentMl;
                *offsetPtr = OFFSET_TO_OFFBASE(curr - matchIndex);
                if (ip+currentMl == iLimit) break; /* best possible, avoids read overflow on next attempt */
            }
        }
    }

    assert(nbAttempts <= (1U << ZSTD_SEARCHLOG_MAX)); /* Check we haven't underflowed. */
    if (dictMode == ZSTD_dedicatedDictSearch) {
        ml = ZSTD_dedicatedDictSearch_lazy_search(offsetPtr, ml, nbAttempts + ddsExtraAttempts, dms,
                                                  ip, iLimit, prefixStart, curr, dictLimit, ddsIdx);
    } else if (dictMode == ZSTD_dictMatchState) {
        /* TODO: Measure and potentially add prefetching to DMS */
        const U32 dmsLowestIndex       = dms->window.dictLimit;
        const BYTE* const dmsBase      = dms->window.base;
        const BYTE* const dmsEnd       = dms->window.nextSrc;
        const U32 dmsSize              = (U32)(dmsEnd - dmsBase);
        const U32 dmsIndexDelta        = dictLimit - dmsSize;

        {   U32 const headGrouped = (*dmsTagRow & rowMask) * groupWidth;
            U32 matchBuffer[ZSTD_ROW_HASH_MAX_ENTRIES];
            size_t numMatches = 0;
            size_t currMatch = 0;
            ZSTD_VecMask matches = ZSTD_row_getMatchMask(dmsTagRow, (BYTE)dmsTag, headGrouped, rowEntries);

            for (; (matches > 0) && (nbAttempts > 0); matches &= (matches - 1)) {
                U32 const matchPos = ((headGrouped + ZSTD_VecMask_next(matches)) / groupWidth) & rowMask;
                U32 const matchIndex = dmsRow[matchPos];
                if(matchPos == 0) continue;
                if (matchIndex < dmsLowestIndex)
                    break;
                PREFETCH_L1(dmsBase + matchIndex);
                matchBuffer[numMatches++] = matchIndex;
                --nbAttempts;
            }

            /* Return the longest match */
            for (; currMatch < numMatches; ++currMatch) {
                U32 const matchIndex = matchBuffer[currMatch];
                size_t currentMl=0;
                assert(matchIndex >= dmsLowestIndex);
                assert(matchIndex < curr);

                {   const BYTE* const match = dmsBase + matchIndex;
                    assert(match+4 <= dmsEnd);
                    if (MEM_read32(match) == MEM_read32(ip))
                        currentMl = ZSTD_count_2segments(ip+4, match+4, iLimit, dmsEnd, prefixStart) + 4;
                }

                if (currentMl > ml) {
                    ml = currentMl;
                    assert(curr > matchIndex + dmsIndexDelta);
                    *offsetPtr = OFFSET_TO_OFFBASE(curr - (matchIndex + dmsIndexDelta));
                    if (ip+currentMl == iLimit) break;
                }
            }
        }
    }
    return ml;
}


/**
 * Generate search functions templated on (dictMode, mls, rowLog).
 * These functions are outlined for code size & compilation time.
 * ZSTD_searchMax() dispatches to the correct implementation function.
 *
 * TODO: The start of the search function involves loading and calculating a
 * bunch of constants from the ZSTD_matchState_t. These computations could be
 * done in an initialization function, and saved somewhere in the match state.
 * Then we could pass a pointer to the saved state instead of the match state,
 * and avoid duplicate computations.
 *
 * TODO: Move the match re-winding into searchMax. This improves compression
 * ratio, and unlocks further simplifications with the next TODO.
 *
 * TODO: Try moving the repcode search into searchMax. After the re-winding
 * and repcode search are in searchMax, there is no more logic in the match
 * finder loop that requires knowledge about the dictMode. So we should be
 * able to avoid force inlining it, and we can join the extDict loop with
 * the single segment loop. It should go in searchMax instead of its own
 * function to avoid having multiple virtual function calls per search.
 */

#define ZSTD_BT_SEARCH_FN(dictMode, mls) ZSTD_BtFindBestMatch_##dictMode##_##mls
#define ZSTD_HC_SEARCH_FN(dictMode, mls) ZSTD_HcFindBestMatch_##dictMode##_##mls
#define ZSTD_ROW_SEARCH_FN(dictMode, mls, rowLog) ZSTD_RowFindBestMatch_##dictMode##_##mls##_##rowLog

#define ZSTD_SEARCH_FN_ATTRS FORCE_NOINLINE

#define GEN_ZSTD_BT_SEARCH_FN(dictMode, mls)                                           \
    ZSTD_SEARCH_FN_ATTRS size_t ZSTD_BT_SEARCH_FN(dictMode, mls)(                      \
            ZSTD_matchState_t* ms,                                                     \
            const BYTE* ip, const BYTE* const iLimit,                                  \
            size_t* offBasePtr)                                                        \
    {                                                                                  \
        assert(MAX(4, MIN(6, ms->cParams.minMatch)) == mls);                           \
        return ZSTD_BtFindBestMatch(ms, ip, iLimit, offBasePtr, mls, ZSTD_##dictMode); \
    }                                                                                  \

#define GEN_ZSTD_HC_SEARCH_FN(dictMode, mls)                                          \
    ZSTD_SEARCH_FN_ATTRS size_t ZSTD_HC_SEARCH_FN(dictMode, mls)(                     \
            ZSTD_matchState_t* ms,                                                    \
            const BYTE* ip, const BYTE* const iLimit,                                 \
            size_t* offsetPtr)                                                        \
    {                                                                                 \
        assert(MAX(4, MIN(6, ms->cParams.minMatch)) == mls);                          \
        return ZSTD_HcFindBestMatch(ms, ip, iLimit, offsetPtr, mls, ZSTD_##dictMode); \
    }                                                                                 \

#define GEN_ZSTD_ROW_SEARCH_FN(dictMode, mls, rowLog)                                          \
    ZSTD_SEARCH_FN_ATTRS size_t ZSTD_ROW_SEARCH_FN(dictMode, mls, rowLog)(                     \
            ZSTD_matchState_t* ms,                                                             \
            const BYTE* ip, const BYTE* const iLimit,                                          \
            size_t* offsetPtr)                                                                 \
    {                                                                                          \
        assert(MAX(4, MIN(6, ms->cParams.minMatch)) == mls);                                   \
        assert(MAX(4, MIN(6, ms->cParams.searchLog)) == rowLog);                               \
        return ZSTD_RowFindBestMatch(ms, ip, iLimit, offsetPtr, mls, ZSTD_##dictMode, rowLog); \
    }                                                                                          \

#define ZSTD_FOR_EACH_ROWLOG(X, dictMode, mls) \
    X(dictMode, mls, 4)                        \
    X(dictMode, mls, 5)                        \
    X(dictMode, mls, 6)

#define ZSTD_FOR_EACH_MLS_ROWLOG(X, dictMode) \
    ZSTD_FOR_EACH_ROWLOG(X, dictMode, 4)      \
    ZSTD_FOR_EACH_ROWLOG(X, dictMode, 5)      \
    ZSTD_FOR_EACH_ROWLOG(X, dictMode, 6)

#define ZSTD_FOR_EACH_MLS(X, dictMode) \
    X(dictMode, 4)                     \
    X(dictMode, 5)                     \
    X(dictMode, 6)

#define ZSTD_FOR_EACH_DICT_MODE(X, ...) \
    X(__VA_ARGS__, noDict)              \
    X(__VA_ARGS__, extDict)             \
    X(__VA_ARGS__, dictMatchState)      \
    X(__VA_ARGS__, dedicatedDictSearch)

/* Generate row search fns for each combination of (dictMode, mls, rowLog) */
ZSTD_FOR_EACH_DICT_MODE(ZSTD_FOR_EACH_MLS_ROWLOG, GEN_ZSTD_ROW_SEARCH_FN)
/* Generate binary Tree search fns for each combination of (dictMode, mls) */
ZSTD_FOR_EACH_DICT_MODE(ZSTD_FOR_EACH_MLS, GEN_ZSTD_BT_SEARCH_FN)
/* Generate hash chain search fns for each combination of (dictMode, mls) */
ZSTD_FOR_EACH_DICT_MODE(ZSTD_FOR_EACH_MLS, GEN_ZSTD_HC_SEARCH_FN)

typedef enum { search_hashChain=0, search_binaryTree=1, search_rowHash=2 } searchMethod_e;

#define GEN_ZSTD_CALL_BT_SEARCH_FN(dictMode, mls)                         \
    case mls:                                                             \
        return ZSTD_BT_SEARCH_FN(dictMode, mls)(ms, ip, iend, offsetPtr);
#define GEN_ZSTD_CALL_HC_SEARCH_FN(dictMode, mls)                         \
    case mls:                                                             \
        return ZSTD_HC_SEARCH_FN(dictMode, mls)(ms, ip, iend, offsetPtr);
#define GEN_ZSTD_CALL_ROW_SEARCH_FN(dictMode, mls, rowLog)                         \
    case rowLog:                                                                   \
        return ZSTD_ROW_SEARCH_FN(dictMode, mls, rowLog)(ms, ip, iend, offsetPtr);

#define ZSTD_SWITCH_MLS(X, dictMode)   \
    switch (mls) {                     \
        ZSTD_FOR_EACH_MLS(X, dictMode) \
    }

#define ZSTD_SWITCH_ROWLOG(dictMode, mls)                                    \
    case mls:                                                                \
        switch (rowLog) {                                                    \
            ZSTD_FOR_EACH_ROWLOG(GEN_ZSTD_CALL_ROW_SEARCH_FN, dictMode, mls) \
        }                                                                    \
        ZSTD_UNREACHABLE;                                                    \
        break;

#define ZSTD_SWITCH_SEARCH_METHOD(dictMode)                       \
    switch (searchMethod) {                                       \
        case search_hashChain:                                    \
            ZSTD_SWITCH_MLS(GEN_ZSTD_CALL_HC_SEARCH_FN, dictMode) \
            break;                                                \
        case search_binaryTree:                                   \
            ZSTD_SWITCH_MLS(GEN_ZSTD_CALL_BT_SEARCH_FN, dictMode) \
            break;                                                \
        case search_rowHash:                                      \
            ZSTD_SWITCH_MLS(ZSTD_SWITCH_ROWLOG, dictMode)         \
            break;                                                \
    }                                                             \
    ZSTD_UNREACHABLE;

/**
 * Searches for the longest match at @p ip.
 * Dispatches to the correct implementation function based on the
 * (searchMethod, dictMode, mls, rowLog). We use switch statements
 * here instead of using an indirect function call through a function
 * pointer because after Spectre and Meltdown mitigations, indirect
 * function calls can be very costly, especially in the kernel.
 *
 * NOTE: dictMode and searchMethod should be templated, so those switch
 * statements should be optimized out. Only the mls & rowLog switches
 * should be left.
 *
 * @param ms The match state.
 * @param ip The position to search at.
 * @param iend The end of the input data.
 * @param[out] offsetPtr Stores the match offset into this pointer.
 * @param mls The minimum search length, in the range [4, 6].
 * @param rowLog The row log (if applicable), in the range [4, 6].
 * @param searchMethod The search method to use (templated).
 * @param dictMode The dictMode (templated).
 *
 * @returns The length of the longest match found, or < mls if no match is found.
 * If a match is found its offset is stored in @p offsetPtr.
 */
FORCE_INLINE_TEMPLATE size_t ZSTD_searchMax(
    ZSTD_matchState_t* ms,
    const BYTE* ip,
    const BYTE* iend,
    size_t* offsetPtr,
    U32 const mls,
    U32 const rowLog,
    searchMethod_e const searchMethod,
    ZSTD_dictMode_e const dictMode)
{
    if (dictMode == ZSTD_noDict) {
        ZSTD_SWITCH_SEARCH_METHOD(noDict)
    } else if (dictMode == ZSTD_extDict) {
        ZSTD_SWITCH_SEARCH_METHOD(extDict)
    } else if (dictMode == ZSTD_dictMatchState) {
        ZSTD_SWITCH_SEARCH_METHOD(dictMatchState)
    } else if (dictMode == ZSTD_dedicatedDictSearch) {
        ZSTD_SWITCH_SEARCH_METHOD(dedicatedDictSearch)
    }
    ZSTD_UNREACHABLE;
    return 0;
}

/* *******************************
*  Common parser - lazy strategy
*********************************/

FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_compressBlock_lazy_generic(
                        ZSTD_matchState_t* ms, seqStore_t* seqStore,
                        U32 rep[ZSTD_REP_NUM],
                        const void* src, size_t srcSize,
                        const searchMethod_e searchMethod, const U32 depth,
                        ZSTD_dictMode_e const dictMode)
{
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = (searchMethod == search_rowHash) ? iend - 8 - ZSTD_ROW_HASH_CACHE_SIZE : iend - 8;
    const BYTE* const base = ms->window.base;
    const U32 prefixLowestIndex = ms->window.dictLimit;
    const BYTE* const prefixLowest = base + prefixLowestIndex;
    const U32 mls = BOUNDED(4, ms->cParams.minMatch, 6);
    const U32 rowLog = BOUNDED(4, ms->cParams.searchLog, 6);

    U32 offset_1 = rep[0], offset_2 = rep[1];
    U32 offsetSaved1 = 0, offsetSaved2 = 0;

    const int isDMS = dictMode == ZSTD_dictMatchState;
    const int isDDS = dictMode == ZSTD_dedicatedDictSearch;
    const int isDxS = isDMS || isDDS;
    const ZSTD_matchState_t* const dms = ms->dictMatchState;
    const U32 dictLowestIndex      = isDxS ? dms->window.dictLimit : 0;
    const BYTE* const dictBase     = isDxS ? dms->window.base : NULL;
    const BYTE* const dictLowest   = isDxS ? dictBase + dictLowestIndex : NULL;
    const BYTE* const dictEnd      = isDxS ? dms->window.nextSrc : NULL;
    const U32 dictIndexDelta       = isDxS ?
                                     prefixLowestIndex - (U32)(dictEnd - dictBase) :
                                     0;
    const U32 dictAndPrefixLength = (U32)((ip - prefixLowest) + (dictEnd - dictLowest));

    DEBUGLOG(5, "ZSTD_compressBlock_lazy_generic (dictMode=%u) (searchFunc=%u)", (U32)dictMode, (U32)searchMethod);
    ip += (dictAndPrefixLength == 0);
    if (dictMode == ZSTD_noDict) {
        U32 const curr = (U32)(ip - base);
        U32 const windowLow = ZSTD_getLowestPrefixIndex(ms, curr, ms->cParams.windowLog);
        U32 const maxRep = curr - windowLow;
        if (offset_2 > maxRep) offsetSaved2 = offset_2, offset_2 = 0;
        if (offset_1 > maxRep) offsetSaved1 = offset_1, offset_1 = 0;
    }
    if (isDxS) {
        /* dictMatchState repCode checks don't currently handle repCode == 0
         * disabling. */
        assert(offset_1 <= dictAndPrefixLength);
        assert(offset_2 <= dictAndPrefixLength);
    }

    /* Reset the lazy skipping state */
    ms->lazySkipping = 0;

    if (searchMethod == search_rowHash) {
        ZSTD_row_fillHashCache(ms, base, rowLog, mls, ms->nextToUpdate, ilimit);
    }

    /* Match Loop */
#if defined(__GNUC__) && defined(__x86_64__)
    /* I've measured random a 5% speed loss on levels 5 & 6 (greedy) when the
     * code alignment is perturbed. To fix the instability align the loop on 32-bytes.
     */
    __asm__(".p2align 5");
#endif
    while (ip < ilimit) {
        size_t matchLength=0;
        size_t offBase = REPCODE1_TO_OFFBASE;
        const BYTE* start=ip+1;
        DEBUGLOG(7, "search baseline (depth 0)");

        /* check repCode */
        if (isDxS) {
            const U32 repIndex = (U32)(ip - base) + 1 - offset_1;
            const BYTE* repMatch = ((dictMode == ZSTD_dictMatchState || dictMode == ZSTD_dedicatedDictSearch)
                                && repIndex < prefixLowestIndex) ?
                                   dictBase + (repIndex - dictIndexDelta) :
                                   base + repIndex;
            if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                && (MEM_read32(repMatch) == MEM_read32(ip+1)) ) {
                const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                if (depth==0) goto _storeSequence;
            }
        }
        if ( dictMode == ZSTD_noDict
          && ((offset_1 > 0) & (MEM_read32(ip+1-offset_1) == MEM_read32(ip+1)))) {
            matchLength = ZSTD_count(ip+1+4, ip+1+4-offset_1, iend) + 4;
            if (depth==0) goto _storeSequence;
        }

        /* first search (depth 0) */
        {   size_t offbaseFound = 999999999;
            size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &offbaseFound, mls, rowLog, searchMethod, dictMode);
            if (ml2 > matchLength)
                matchLength = ml2, start = ip, offBase = offbaseFound;
        }

        if (matchLength < 4) {
            size_t const step = ((size_t)(ip-anchor) >> kSearchStrength) + 1;   /* jump faster over incompressible sections */;
            ip += step;
            /* Enter the lazy skipping mode once we are skipping more than 8 bytes at a time.
             * In this mode we stop inserting every position into our tables, and only insert
             * positions that we search, which is one in step positions.
             * The exact cutoff is flexible, I've just chosen a number that is reasonably high,
             * so we minimize the compression ratio loss in "normal" scenarios. This mode gets
             * triggered once we've gone 2KB without finding any matches.
             */
            ms->lazySkipping = step > kLazySkippingStep;
            continue;
        }

        /* let's try to find a better solution */
        if (depth>=1)
        while (ip<ilimit) {
            DEBUGLOG(7, "search depth 1");
            ip ++;
            if ( (dictMode == ZSTD_noDict)
              && (offBase) && ((offset_1>0) & (MEM_read32(ip) == MEM_read32(ip - offset_1)))) {
                size_t const mlRep = ZSTD_count(ip+4, ip+4-offset_1, iend) + 4;
                int const gain2 = (int)(mlRep * 3);
                int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offBase) + 1);
                if ((mlRep >= 4) && (gain2 > gain1))
                    matchLength = mlRep, offBase = REPCODE1_TO_OFFBASE, start = ip;
            }
            if (isDxS) {
                const U32 repIndex = (U32)(ip - base) - offset_1;
                const BYTE* repMatch = repIndex < prefixLowestIndex ?
                               dictBase + (repIndex - dictIndexDelta) :
                               base + repIndex;
                if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                    && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                    const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                    size_t const mlRep = ZSTD_count_2segments(ip+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                    int const gain2 = (int)(mlRep * 3);
                    int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offBase) + 1);
                    if ((mlRep >= 4) && (gain2 > gain1))
                        matchLength = mlRep, offBase = REPCODE1_TO_OFFBASE, start = ip;
                }
            }
            {   size_t ofbCandidate=999999999;
                size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &ofbCandidate, mls, rowLog, searchMethod, dictMode);
                int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)ofbCandidate));   /* raw approx */
                int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 4);
                if ((ml2 >= 4) && (gain2 > gain1)) {
                    matchLength = ml2, offBase = ofbCandidate, start = ip;
                    continue;   /* search a better one */
            }   }

            /* let's find an even better one */
            if ((depth==2) && (ip<ilimit)) {
                DEBUGLOG(7, "search depth 2");
                ip ++;
                if ( (dictMode == ZSTD_noDict)
                  && (offBase) && ((offset_1>0) & (MEM_read32(ip) == MEM_read32(ip - offset_1)))) {
                    size_t const mlRep = ZSTD_count(ip+4, ip+4-offset_1, iend) + 4;
                    int const gain2 = (int)(mlRep * 4);
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 1);
                    if ((mlRep >= 4) && (gain2 > gain1))
                        matchLength = mlRep, offBase = REPCODE1_TO_OFFBASE, start = ip;
                }
                if (isDxS) {
                    const U32 repIndex = (U32)(ip - base) - offset_1;
                    const BYTE* repMatch = repIndex < prefixLowestIndex ?
                                   dictBase + (repIndex - dictIndexDelta) :
                                   base + repIndex;
                    if (((U32)((prefixLowestIndex-1) - repIndex) >= 3 /* intentional underflow */)
                        && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                        const BYTE* repMatchEnd = repIndex < prefixLowestIndex ? dictEnd : iend;
                        size_t const mlRep = ZSTD_count_2segments(ip+4, repMatch+4, iend, repMatchEnd, prefixLowest) + 4;
                        int const gain2 = (int)(mlRep * 4);
                        int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 1);
                        if ((mlRep >= 4) && (gain2 > gain1))
                            matchLength = mlRep, offBase = REPCODE1_TO_OFFBASE, start = ip;
                    }
                }
                {   size_t ofbCandidate=999999999;
                    size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &ofbCandidate, mls, rowLog, searchMethod, dictMode);
                    int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)ofbCandidate));   /* raw approx */
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 7);
                    if ((ml2 >= 4) && (gain2 > gain1)) {
                        matchLength = ml2, offBase = ofbCandidate, start = ip;
                        continue;
            }   }   }
            break;  /* nothing found : store previous solution */
        }

        /* NOTE:
         * Pay attention that `start[-value]` can lead to strange undefined behavior
         * notably if `value` is unsigned, resulting in a large positive `-value`.
         */
        /* catch up */
        if (OFFBASE_IS_OFFSET(offBase)) {
            if (dictMode == ZSTD_noDict) {
                while ( ((start > anchor) & (start - OFFBASE_TO_OFFSET(offBase) > prefixLowest))
                     && (start[-1] == (start-OFFBASE_TO_OFFSET(offBase))[-1]) )  /* only search for offset within prefix */
                    { start--; matchLength++; }
            }
            if (isDxS) {
                U32 const matchIndex = (U32)((size_t)(start-base) - OFFBASE_TO_OFFSET(offBase));
                const BYTE* match = (matchIndex < prefixLowestIndex) ? dictBase + matchIndex - dictIndexDelta : base + matchIndex;
                const BYTE* const mStart = (matchIndex < prefixLowestIndex) ? dictLowest : prefixLowest;
                while ((start>anchor) && (match>mStart) && (start[-1] == match[-1])) { start--; match--; matchLength++; }  /* catch up */
            }
            offset_2 = offset_1; offset_1 = (U32)OFFBASE_TO_OFFSET(offBase);
        }
        /* store sequence */
_storeSequence:
        {   size_t const litLength = (size_t)(start - anchor);
            ZSTD_storeSeq(seqStore, litLength, anchor, iend, (U32)offBase, matchLength);
            anchor = ip = start + matchLength;
        }
        if (ms->lazySkipping) {
            /* We've found a match, disable lazy skipping mode, and refill the hash cache. */
            if (searchMethod == search_rowHash) {
                ZSTD_row_fillHashCache(ms, base, rowLog, mls, ms->nextToUpdate, ilimit);
            }
            ms->lazySkipping = 0;
        }

        /* check immediate repcode */
        if (isDxS) {
            while (ip <= ilimit) {
                U32 const current2 = (U32)(ip-base);
                U32 const repIndex = current2 - offset_2;
                const BYTE* repMatch = repIndex < prefixLowestIndex ?
                        dictBase - dictIndexDelta + repIndex :
                        base + repIndex;
                if ( ((U32)((prefixLowestIndex-1) - (U32)repIndex) >= 3 /* intentional overflow */)
                   && (MEM_read32(repMatch) == MEM_read32(ip)) ) {
                    const BYTE* const repEnd2 = repIndex < prefixLowestIndex ? dictEnd : iend;
                    matchLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd2, prefixLowest) + 4;
                    offBase = offset_2; offset_2 = offset_1; offset_1 = (U32)offBase;   /* swap offset_2 <=> offset_1 */
                    ZSTD_storeSeq(seqStore, 0, anchor, iend, REPCODE1_TO_OFFBASE, matchLength);
                    ip += matchLength;
                    anchor = ip;
                    continue;
                }
                break;
            }
        }

        if (dictMode == ZSTD_noDict) {
            while ( ((ip <= ilimit) & (offset_2>0))
                 && (MEM_read32(ip) == MEM_read32(ip - offset_2)) ) {
                /* store sequence */
                matchLength = ZSTD_count(ip+4, ip+4-offset_2, iend) + 4;
                offBase = offset_2; offset_2 = offset_1; offset_1 = (U32)offBase; /* swap repcodes */
                ZSTD_storeSeq(seqStore, 0, anchor, iend, REPCODE1_TO_OFFBASE, matchLength);
                ip += matchLength;
                anchor = ip;
                continue;   /* faster when present ... (?) */
    }   }   }

    /* If offset_1 started invalid (offsetSaved1 != 0) and became valid (offset_1 != 0),
     * rotate saved offsets. See comment in ZSTD_compressBlock_fast_noDict for more context. */
    offsetSaved2 = ((offsetSaved1 != 0) && (offset_1 != 0)) ? offsetSaved1 : offsetSaved2;

    /* save reps for next block */
    rep[0] = offset_1 ? offset_1 : offsetSaved1;
    rep[1] = offset_2 ? offset_2 : offsetSaved2;

    /* Return the last literals size */
    return (size_t)(iend - anchor);
}
#endif /* build exclusions */


#ifndef ZSTD_EXCLUDE_GREEDY_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_greedy(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 0, ZSTD_noDict);
}

size_t ZSTD_compressBlock_greedy_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 0, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_greedy_dedicatedDictSearch(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 0, ZSTD_dedicatedDictSearch);
}

size_t ZSTD_compressBlock_greedy_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 0, ZSTD_noDict);
}

size_t ZSTD_compressBlock_greedy_dictMatchState_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 0, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_greedy_dedicatedDictSearch_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 0, ZSTD_dedicatedDictSearch);
}
#endif

#ifndef ZSTD_EXCLUDE_LAZY_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_lazy(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 1, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 1, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy_dedicatedDictSearch(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 1, ZSTD_dedicatedDictSearch);
}

size_t ZSTD_compressBlock_lazy_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 1, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy_dictMatchState_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 1, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy_dedicatedDictSearch_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 1, ZSTD_dedicatedDictSearch);
}
#endif

#ifndef ZSTD_EXCLUDE_LAZY2_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_lazy2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 2, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 2, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy2_dedicatedDictSearch(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 2, ZSTD_dedicatedDictSearch);
}

size_t ZSTD_compressBlock_lazy2_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 2, ZSTD_noDict);
}

size_t ZSTD_compressBlock_lazy2_dictMatchState_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 2, ZSTD_dictMatchState);
}

size_t ZSTD_compressBlock_lazy2_dedicatedDictSearch_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 2, ZSTD_dedicatedDictSearch);
}
#endif

#ifndef ZSTD_EXCLUDE_BTLAZY2_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_btlazy2(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_binaryTree, 2, ZSTD_noDict);
}

size_t ZSTD_compressBlock_btlazy2_dictMatchState(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_generic(ms, seqStore, rep, src, srcSize, search_binaryTree, 2, ZSTD_dictMatchState);
}
#endif

#if !defined(ZSTD_EXCLUDE_GREEDY_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_LAZY_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_LAZY2_BLOCK_COMPRESSOR) \
 || !defined(ZSTD_EXCLUDE_BTLAZY2_BLOCK_COMPRESSOR)
FORCE_INLINE_TEMPLATE
ZSTD_ALLOW_POINTER_OVERFLOW_ATTR
size_t ZSTD_compressBlock_lazy_extDict_generic(
                        ZSTD_matchState_t* ms, seqStore_t* seqStore,
                        U32 rep[ZSTD_REP_NUM],
                        const void* src, size_t srcSize,
                        const searchMethod_e searchMethod, const U32 depth)
{
    const BYTE* const istart = (const BYTE*)src;
    const BYTE* ip = istart;
    const BYTE* anchor = istart;
    const BYTE* const iend = istart + srcSize;
    const BYTE* const ilimit = searchMethod == search_rowHash ? iend - 8 - ZSTD_ROW_HASH_CACHE_SIZE : iend - 8;
    const BYTE* const base = ms->window.base;
    const U32 dictLimit = ms->window.dictLimit;
    const BYTE* const prefixStart = base + dictLimit;
    const BYTE* const dictBase = ms->window.dictBase;
    const BYTE* const dictEnd  = dictBase + dictLimit;
    const BYTE* const dictStart  = dictBase + ms->window.lowLimit;
    const U32 windowLog = ms->cParams.windowLog;
    const U32 mls = BOUNDED(4, ms->cParams.minMatch, 6);
    const U32 rowLog = BOUNDED(4, ms->cParams.searchLog, 6);

    U32 offset_1 = rep[0], offset_2 = rep[1];

    DEBUGLOG(5, "ZSTD_compressBlock_lazy_extDict_generic (searchFunc=%u)", (U32)searchMethod);

    /* Reset the lazy skipping state */
    ms->lazySkipping = 0;

    /* init */
    ip += (ip == prefixStart);
    if (searchMethod == search_rowHash) {
        ZSTD_row_fillHashCache(ms, base, rowLog, mls, ms->nextToUpdate, ilimit);
    }

    /* Match Loop */
#if defined(__GNUC__) && defined(__x86_64__)
    /* I've measured random a 5% speed loss on levels 5 & 6 (greedy) when the
     * code alignment is perturbed. To fix the instability align the loop on 32-bytes.
     */
    __asm__(".p2align 5");
#endif
    while (ip < ilimit) {
        size_t matchLength=0;
        size_t offBase = REPCODE1_TO_OFFBASE;
        const BYTE* start=ip+1;
        U32 curr = (U32)(ip-base);

        /* check repCode */
        {   const U32 windowLow = ZSTD_getLowestMatchIndex(ms, curr+1, windowLog);
            const U32 repIndex = (U32)(curr+1 - offset_1);
            const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
            const BYTE* const repMatch = repBase + repIndex;
            if ( ((U32)((dictLimit-1) - repIndex) >= 3) /* intentional overflow */
               & (offset_1 <= curr+1 - windowLow) ) /* note: we are searching at curr+1 */
            if (MEM_read32(ip+1) == MEM_read32(repMatch)) {
                /* repcode detected we should take it */
                const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+1+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                if (depth==0) goto _storeSequence;
        }   }

        /* first search (depth 0) */
        {   size_t ofbCandidate = 999999999;
            size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &ofbCandidate, mls, rowLog, searchMethod, ZSTD_extDict);
            if (ml2 > matchLength)
                matchLength = ml2, start = ip, offBase = ofbCandidate;
        }

        if (matchLength < 4) {
            size_t const step = ((size_t)(ip-anchor) >> kSearchStrength);
            ip += step + 1;   /* jump faster over incompressible sections */
            /* Enter the lazy skipping mode once we are skipping more than 8 bytes at a time.
             * In this mode we stop inserting every position into our tables, and only insert
             * positions that we search, which is one in step positions.
             * The exact cutoff is flexible, I've just chosen a number that is reasonably high,
             * so we minimize the compression ratio loss in "normal" scenarios. This mode gets
             * triggered once we've gone 2KB without finding any matches.
             */
            ms->lazySkipping = step > kLazySkippingStep;
            continue;
        }

        /* let's try to find a better solution */
        if (depth>=1)
        while (ip<ilimit) {
            ip ++;
            curr++;
            /* check repCode */
            if (offBase) {
                const U32 windowLow = ZSTD_getLowestMatchIndex(ms, curr, windowLog);
                const U32 repIndex = (U32)(curr - offset_1);
                const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
                const BYTE* const repMatch = repBase + repIndex;
                if ( ((U32)((dictLimit-1) - repIndex) >= 3) /* intentional overflow : do not test positions overlapping 2 memory segments  */
                   & (offset_1 <= curr - windowLow) ) /* equivalent to `curr > repIndex >= windowLow` */
                if (MEM_read32(ip) == MEM_read32(repMatch)) {
                    /* repcode detected */
                    const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                    size_t const repLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                    int const gain2 = (int)(repLength * 3);
                    int const gain1 = (int)(matchLength*3 - ZSTD_highbit32((U32)offBase) + 1);
                    if ((repLength >= 4) && (gain2 > gain1))
                        matchLength = repLength, offBase = REPCODE1_TO_OFFBASE, start = ip;
            }   }

            /* search match, depth 1 */
            {   size_t ofbCandidate = 999999999;
                size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &ofbCandidate, mls, rowLog, searchMethod, ZSTD_extDict);
                int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)ofbCandidate));   /* raw approx */
                int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 4);
                if ((ml2 >= 4) && (gain2 > gain1)) {
                    matchLength = ml2, offBase = ofbCandidate, start = ip;
                    continue;   /* search a better one */
            }   }

            /* let's find an even better one */
            if ((depth==2) && (ip<ilimit)) {
                ip ++;
                curr++;
                /* check repCode */
                if (offBase) {
                    const U32 windowLow = ZSTD_getLowestMatchIndex(ms, curr, windowLog);
                    const U32 repIndex = (U32)(curr - offset_1);
                    const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
                    const BYTE* const repMatch = repBase + repIndex;
                    if ( ((U32)((dictLimit-1) - repIndex) >= 3) /* intentional overflow : do not test positions overlapping 2 memory segments  */
                       & (offset_1 <= curr - windowLow) ) /* equivalent to `curr > repIndex >= windowLow` */
                    if (MEM_read32(ip) == MEM_read32(repMatch)) {
                        /* repcode detected */
                        const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                        size_t const repLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                        int const gain2 = (int)(repLength * 4);
                        int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 1);
                        if ((repLength >= 4) && (gain2 > gain1))
                            matchLength = repLength, offBase = REPCODE1_TO_OFFBASE, start = ip;
                }   }

                /* search match, depth 2 */
                {   size_t ofbCandidate = 999999999;
                    size_t const ml2 = ZSTD_searchMax(ms, ip, iend, &ofbCandidate, mls, rowLog, searchMethod, ZSTD_extDict);
                    int const gain2 = (int)(ml2*4 - ZSTD_highbit32((U32)ofbCandidate));   /* raw approx */
                    int const gain1 = (int)(matchLength*4 - ZSTD_highbit32((U32)offBase) + 7);
                    if ((ml2 >= 4) && (gain2 > gain1)) {
                        matchLength = ml2, offBase = ofbCandidate, start = ip;
                        continue;
            }   }   }
            break;  /* nothing found : store previous solution */
        }

        /* catch up */
        if (OFFBASE_IS_OFFSET(offBase)) {
            U32 const matchIndex = (U32)((size_t)(start-base) - OFFBASE_TO_OFFSET(offBase));
            const BYTE* match = (matchIndex < dictLimit) ? dictBase + matchIndex : base + matchIndex;
            const BYTE* const mStart = (matchIndex < dictLimit) ? dictStart : prefixStart;
            while ((start>anchor) && (match>mStart) && (start[-1] == match[-1])) { start--; match--; matchLength++; }  /* catch up */
            offset_2 = offset_1; offset_1 = (U32)OFFBASE_TO_OFFSET(offBase);
        }

        /* store sequence */
_storeSequence:
        {   size_t const litLength = (size_t)(start - anchor);
            ZSTD_storeSeq(seqStore, litLength, anchor, iend, (U32)offBase, matchLength);
            anchor = ip = start + matchLength;
        }
        if (ms->lazySkipping) {
            /* We've found a match, disable lazy skipping mode, and refill the hash cache. */
            if (searchMethod == search_rowHash) {
                ZSTD_row_fillHashCache(ms, base, rowLog, mls, ms->nextToUpdate, ilimit);
            }
            ms->lazySkipping = 0;
        }

        /* check immediate repcode */
        while (ip <= ilimit) {
            const U32 repCurrent = (U32)(ip-base);
            const U32 windowLow = ZSTD_getLowestMatchIndex(ms, repCurrent, windowLog);
            const U32 repIndex = repCurrent - offset_2;
            const BYTE* const repBase = repIndex < dictLimit ? dictBase : base;
            const BYTE* const repMatch = repBase + repIndex;
            if ( ((U32)((dictLimit-1) - repIndex) >= 3) /* intentional overflow : do not test positions overlapping 2 memory segments  */
               & (offset_2 <= repCurrent - windowLow) ) /* equivalent to `curr > repIndex >= windowLow` */
            if (MEM_read32(ip) == MEM_read32(repMatch)) {
                /* repcode detected we should take it */
                const BYTE* const repEnd = repIndex < dictLimit ? dictEnd : iend;
                matchLength = ZSTD_count_2segments(ip+4, repMatch+4, iend, repEnd, prefixStart) + 4;
                offBase = offset_2; offset_2 = offset_1; offset_1 = (U32)offBase;   /* swap offset history */
                ZSTD_storeSeq(seqStore, 0, anchor, iend, REPCODE1_TO_OFFBASE, matchLength);
                ip += matchLength;
                anchor = ip;
                continue;   /* faster when present ... (?) */
            }
            break;
    }   }

    /* Save reps for next block */
    rep[0] = offset_1;
    rep[1] = offset_2;

    /* Return the last literals size */
    return (size_t)(iend - anchor);
}
#endif /* build exclusions */

#ifndef ZSTD_EXCLUDE_GREEDY_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_greedy_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 0);
}

size_t ZSTD_compressBlock_greedy_extDict_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 0);
}
#endif

#ifndef ZSTD_EXCLUDE_LAZY_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_lazy_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 1);
}

size_t ZSTD_compressBlock_lazy_extDict_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 1);
}
#endif

#ifndef ZSTD_EXCLUDE_LAZY2_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_lazy2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_hashChain, 2);
}

size_t ZSTD_compressBlock_lazy2_extDict_row(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)
{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_rowHash, 2);
}
#endif

#ifndef ZSTD_EXCLUDE_BTLAZY2_BLOCK_COMPRESSOR
size_t ZSTD_compressBlock_btlazy2_extDict(
        ZSTD_matchState_t* ms, seqStore_t* seqStore, U32 rep[ZSTD_REP_NUM],
        void const* src, size_t srcSize)

{
    return ZSTD_compressBlock_lazy_extDict_generic(ms, seqStore, rep, src, srcSize, search_binaryTree, 2);
}
#endif

/* Dictionary format :
 * See :
 * https://github.com/facebook/zstd/blob/release/doc/zstd_compression_format.md#dictionary-format
 */
/*! ZSTD_loadZstdDictionary() :
 * @return : dictID, or an error code
 *  assumptions : magic number supposed already checked
 *                dictSize supposed >= 8
 */
static size_t ZSTD_loadZstdDictionary(ZSTD_compressedBlockState_t* bs,
                                      ZSTD_matchState_t* ms,
                                      ZSTD_cwksp* ws,
                                      ZSTD_CCtx_params const* params,
                                      const void* dict, size_t dictSize,
                                      ZSTD_dictTableLoadMethod_e dtlm,
                                      ZSTD_tableFillPurpose_e tfp,
                                      void* workspace)
{
    const BYTE* dictPtr = (const BYTE*)dict;
    const BYTE* const dictEnd = dictPtr + dictSize;
    size_t dictID;
    size_t eSize;
    ZSTD_STATIC_ASSERT(HUF_WORKSPACE_SIZE >= (1<<MAX(MLFSELog,LLFSELog)));
    assert(dictSize >= 8);
    assert(MEM_readLE32(dictPtr) == ZSTD_MAGIC_DICTIONARY);

    dictID = params->fParams.noDictIDFlag ? 0 :  MEM_readLE32(dictPtr + 4 /* skip magic number */ );
    eSize = ZSTD_loadCEntropy(bs, workspace, dict, dictSize);
    FORWARD_IF_ERROR(eSize, "ZSTD_loadCEntropy failed");
    dictPtr += eSize;

    {
        size_t const dictContentSize = (size_t)(dictEnd - dictPtr);
        FORWARD_IF_ERROR(ZSTD_loadDictionaryContent(
            ms, NULL, ws, params, dictPtr, dictContentSize, dtlm, tfp), "");
    }
    return dictID;
}

/** ZSTD_compress_insertDictionary() :
*   @return : dictID, or an error code */
static size_t
ZSTD_compress_insertDictionary(ZSTD_compressedBlockState_t* bs,
                               ZSTD_matchState_t* ms,
                               ldmState_t* ls,
                               ZSTD_cwksp* ws,
                         const ZSTD_CCtx_params* params,
                         const void* dict, size_t dictSize,
                               ZSTD_dictContentType_e dictContentType,
                               ZSTD_dictTableLoadMethod_e dtlm,
                               ZSTD_tableFillPurpose_e tfp,
                               void* workspace)
{
    DEBUGLOG(4, "ZSTD_compress_insertDictionary (dictSize=%u)", (U32)dictSize);
    if ((dict==NULL) || (dictSize<8)) {
        RETURN_ERROR_IF(dictContentType == ZSTD_dct_fullDict, dictionary_wrong, "");
        return 0;
    }

    ZSTD_reset_compressedBlockState(bs);

    /* dict restricted modes */
    if (dictContentType == ZSTD_dct_rawContent)
        return ZSTD_loadDictionaryContent(ms, ls, ws, params, dict, dictSize, dtlm, tfp);

    if (MEM_readLE32(dict) != ZSTD_MAGIC_DICTIONARY) {
        if (dictContentType == ZSTD_dct_auto) {
            DEBUGLOG(4, "raw content dictionary detected");
            return ZSTD_loadDictionaryContent(
                ms, ls, ws, params, dict, dictSize, dtlm, tfp);
        }
        RETURN_ERROR_IF(dictContentType == ZSTD_dct_fullDict, dictionary_wrong, "");
        assert(0);   /* impossible */
    }

    /* dict as full zstd dictionary */
    return ZSTD_loadZstdDictionary(
        bs, ms, ws, params, dict, dictSize, dtlm, tfp, workspace);
}

#define ZSTD_USE_CDICT_PARAMS_SRCSIZE_CUTOFF (128 KB)
#define ZSTD_USE_CDICT_PARAMS_DICTSIZE_MULTIPLIER (6ULL)

/*! ZSTD_compressBegin_internal() :
 * Assumption : either @dict OR @cdict (or none) is non-NULL, never both
 * @return : 0, or an error code */
static size_t ZSTD_compressBegin_internal(ZSTD_CCtx* cctx,
                                    const void* dict, size_t dictSize,
                                    ZSTD_dictContentType_e dictContentType,
                                    ZSTD_dictTableLoadMethod_e dtlm,
                                    const ZSTD_CDict* cdict,
                                    const ZSTD_CCtx_params* params, U64 pledgedSrcSize,
                                    ZSTD_buffered_policy_e zbuff)
{
    size_t const dictContentSize = cdict ? cdict->dictContentSize : dictSize;
#if ZSTD_TRACE
    cctx->traceCtx = (ZSTD_trace_compress_begin != NULL) ? ZSTD_trace_compress_begin(cctx) : 0;
#endif
    DEBUGLOG(4, "ZSTD_compressBegin_internal: wlog=%u", params->cParams.windowLog);
    /* params are supposed to be fully validated at this point */
    assert(!ZSTD_isError(ZSTD_checkCParams(params->cParams)));
    assert(!((dict) && (cdict)));  /* either dict or cdict, not both */
    if ( (cdict)
      && (cdict->dictContentSize > 0)
      && ( pledgedSrcSize < ZSTD_USE_CDICT_PARAMS_SRCSIZE_CUTOFF
        || pledgedSrcSize < cdict->dictContentSize * ZSTD_USE_CDICT_PARAMS_DICTSIZE_MULTIPLIER
        || pledgedSrcSize == ZSTD_CONTENTSIZE_UNKNOWN
        || cdict->compressionLevel == 0)
      && (params->attachDictPref != ZSTD_dictForceLoad) ) {
        return ZSTD_resetCCtx_usingCDict(cctx, cdict, params, pledgedSrcSize, zbuff);
    }

    FORWARD_IF_ERROR( ZSTD_resetCCtx_internal(cctx, params, pledgedSrcSize,
                                     dictContentSize,
                                     ZSTDcrp_makeClean, zbuff) , "");
    {   size_t const dictID = cdict ?
                ZSTD_compress_insertDictionary(
                        cctx->blockState.prevCBlock, &cctx->blockState.matchState,
                        &cctx->ldmState, &cctx->workspace, &cctx->appliedParams, cdict->dictContent,
                        cdict->dictContentSize, cdict->dictContentType, dtlm,
                        ZSTD_tfp_forCCtx, cctx->entropyWorkspace)
              : ZSTD_compress_insertDictionary(
                        cctx->blockState.prevCBlock, &cctx->blockState.matchState,
                        &cctx->ldmState, &cctx->workspace, &cctx->appliedParams, dict, dictSize,
                        dictContentType, dtlm, ZSTD_tfp_forCCtx, cctx->entropyWorkspace);
        FORWARD_IF_ERROR(dictID, "ZSTD_compress_insertDictionary failed");
        assert(dictID <= UINT_MAX);
        cctx->dictID = (U32)dictID;
        cctx->dictContentSize = dictContentSize;
    }
    return 0;
}

size_t ZSTD_compressBegin_advanced_internal(ZSTD_CCtx* cctx,
                                    const void* dict, size_t dictSize,
                                    ZSTD_dictContentType_e dictContentType,
                                    ZSTD_dictTableLoadMethod_e dtlm,
                                    const ZSTD_CDict* cdict,
                                    const ZSTD_CCtx_params* params,
                                    unsigned long long pledgedSrcSize)
{
    DEBUGLOG(4, "ZSTD_compressBegin_advanced_internal: wlog=%u", params->cParams.windowLog);
    /* compression parameters verification and optimization */
    FORWARD_IF_ERROR( ZSTD_checkCParams(params->cParams) , "");
    return ZSTD_compressBegin_internal(cctx,
                                       dict, dictSize, dictContentType, dtlm,
                                       cdict,
                                       params, pledgedSrcSize,
                                       ZSTDb_not_buffered);
}

/*! ZSTD_compressBegin_advanced() :
*   @return : 0, or an error code */
size_t ZSTD_compressBegin_advanced(ZSTD_CCtx* cctx,
                             const void* dict, size_t dictSize,
                                   ZSTD_parameters params, unsigned long long pledgedSrcSize)
{
    ZSTD_CCtx_params cctxParams;
    ZSTD_CCtxParams_init_internal(&cctxParams, &params, ZSTD_NO_CLEVEL);
    return ZSTD_compressBegin_advanced_internal(cctx,
                                            dict, dictSize, ZSTD_dct_auto, ZSTD_dtlm_fast,
                                            NULL /*cdict*/,
                                            &cctxParams, pledgedSrcSize);
}

static size_t
ZSTD_compressBegin_usingDict_deprecated(ZSTD_CCtx* cctx, const void* dict, size_t dictSize, int compressionLevel)
{
    ZSTD_CCtx_params cctxParams;
    {   ZSTD_parameters const params = ZSTD_getParams_internal(compressionLevel, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_noAttachDict);
        ZSTD_CCtxParams_init_internal(&cctxParams, &params, (compressionLevel == 0) ? ZSTD_CLEVEL_DEFAULT : compressionLevel);
    }
    DEBUGLOG(4, "ZSTD_compressBegin_usingDict (dictSize=%u)", (unsigned)dictSize);
    return ZSTD_compressBegin_internal(cctx, dict, dictSize, ZSTD_dct_auto, ZSTD_dtlm_fast, NULL,
                                       &cctxParams, ZSTD_CONTENTSIZE_UNKNOWN, ZSTDb_not_buffered);
}

size_t
ZSTD_compressBegin_usingDict(ZSTD_CCtx* cctx, const void* dict, size_t dictSize, int compressionLevel)
{
    return ZSTD_compressBegin_usingDict_deprecated(cctx, dict, dictSize, compressionLevel);
}

size_t ZSTD_compressBegin(ZSTD_CCtx* cctx, int compressionLevel)
{
    return ZSTD_compressBegin_usingDict_deprecated(cctx, NULL, 0, compressionLevel);
}


/*! ZSTD_writeEpilogue() :
*   Ends a frame.
*   @return : nb of bytes written into dst (or an error code) */
static size_t ZSTD_writeEpilogue(ZSTD_CCtx* cctx, void* dst, size_t dstCapacity)
{
    BYTE* const ostart = (BYTE*)dst;
    BYTE* op = ostart;

    DEBUGLOG(4, "ZSTD_writeEpilogue");
    RETURN_ERROR_IF(cctx->stage == ZSTDcs_created, stage_wrong, "init missing");

    /* special case : empty frame */
    if (cctx->stage == ZSTDcs_init) {
        size_t fhSize = ZSTD_writeFrameHeader(dst, dstCapacity, &cctx->appliedParams, 0, 0);
        FORWARD_IF_ERROR(fhSize, "ZSTD_writeFrameHeader failed");
        dstCapacity -= fhSize;
        op += fhSize;
        cctx->stage = ZSTDcs_ongoing;
    }

    if (cctx->stage != ZSTDcs_ending) {
        /* write one last empty block, make it the "last" block */
        U32 const cBlockHeader24 = 1 /* last block */ + (((U32)bt_raw)<<1) + 0;
        ZSTD_STATIC_ASSERT(ZSTD_BLOCKHEADERSIZE == 3);
        RETURN_ERROR_IF(dstCapacity<3, dstSize_tooSmall, "no room for epilogue");
        MEM_writeLE24(op, cBlockHeader24);
        op += ZSTD_blockHeaderSize;
        dstCapacity -= ZSTD_blockHeaderSize;
    }

    if (cctx->appliedParams.fParams.checksumFlag) {
        U32 const checksum = (U32) XXH64_digest(&cctx->xxhState);
        RETURN_ERROR_IF(dstCapacity<4, dstSize_tooSmall, "no room for checksum");
        DEBUGLOG(4, "ZSTD_writeEpilogue: write checksum : %08X", (unsigned)checksum);
        MEM_writeLE32(op, checksum);
        op += 4;
    }

    cctx->stage = ZSTDcs_created;  /* return to "created but no init" status */
    return op-ostart;
}

void ZSTD_CCtx_trace(ZSTD_CCtx* cctx, size_t extraCSize)
{
#if ZSTD_TRACE
    if (cctx->traceCtx && ZSTD_trace_compress_end != NULL) {
        int const streaming = cctx->inBuffSize > 0 || cctx->outBuffSize > 0 || cctx->appliedParams.nbWorkers > 0;
        ZSTD_Trace trace;
        ZSTD_memset(&trace, 0, sizeof(trace));
        trace.version = ZSTD_VERSION_NUMBER;
        trace.streaming = streaming;
        trace.dictionaryID = cctx->dictID;
        trace.dictionarySize = cctx->dictContentSize;
        trace.uncompressedSize = cctx->consumedSrcSize;
        trace.compressedSize = cctx->producedCSize + extraCSize;
        trace.params = &cctx->appliedParams;
        trace.cctx = cctx;
        ZSTD_trace_compress_end(cctx->traceCtx, &trace);
    }
    cctx->traceCtx = 0;
#else
    (void)cctx;
    (void)extraCSize;
#endif
}

size_t ZSTD_compressEnd_public(ZSTD_CCtx* cctx,
                               void* dst, size_t dstCapacity,
                         const void* src, size_t srcSize)
{
    size_t endResult;
    size_t const cSize = ZSTD_compressContinue_internal(cctx,
                                dst, dstCapacity, src, srcSize,
                                1 /* frame mode */, 1 /* last chunk */);
    FORWARD_IF_ERROR(cSize, "ZSTD_compressContinue_internal failed");
    endResult = ZSTD_writeEpilogue(cctx, (char*)dst + cSize, dstCapacity-cSize);
    FORWARD_IF_ERROR(endResult, "ZSTD_writeEpilogue failed");
    assert(!(cctx->appliedParams.fParams.contentSizeFlag && cctx->pledgedSrcSizePlusOne == 0));
    if (cctx->pledgedSrcSizePlusOne != 0) {  /* control src size */
        ZSTD_STATIC_ASSERT(ZSTD_CONTENTSIZE_UNKNOWN == (unsigned long long)-1);
        DEBUGLOG(4, "end of frame : controlling src size");
        RETURN_ERROR_IF(
            cctx->pledgedSrcSizePlusOne != cctx->consumedSrcSize+1,
            srcSize_wrong,
             "error : pledgedSrcSize = %u, while realSrcSize = %u",
            (unsigned)cctx->pledgedSrcSizePlusOne-1,
            (unsigned)cctx->consumedSrcSize);
    }
    ZSTD_CCtx_trace(cctx, endResult);
    return cSize + endResult;
}

/* NOTE: Must just wrap ZSTD_compressEnd_public() */
size_t ZSTD_compressEnd(ZSTD_CCtx* cctx,
                        void* dst, size_t dstCapacity,
                  const void* src, size_t srcSize)
{
    return ZSTD_compressEnd_public(cctx, dst, dstCapacity, src, srcSize);
}

size_t ZSTD_compress_advanced (ZSTD_CCtx* cctx,
                               void* dst, size_t dstCapacity,
                         const void* src, size_t srcSize,
                         const void* dict,size_t dictSize,
                               ZSTD_parameters params)
{
    DEBUGLOG(4, "ZSTD_compress_advanced");
    FORWARD_IF_ERROR(ZSTD_checkCParams(params.cParams), "");
    ZSTD_CCtxParams_init_internal(&cctx->simpleApiParams, &params, ZSTD_NO_CLEVEL);
    return ZSTD_compress_advanced_internal(cctx,
                                           dst, dstCapacity,
                                           src, srcSize,
                                           dict, dictSize,
                                           &cctx->simpleApiParams);
}

/* Internal */
size_t ZSTD_compress_advanced_internal(
        ZSTD_CCtx* cctx,
        void* dst, size_t dstCapacity,
        const void* src, size_t srcSize,
        const void* dict,size_t dictSize,
        const ZSTD_CCtx_params* params)
{
    DEBUGLOG(4, "ZSTD_compress_advanced_internal (srcSize:%u)", (unsigned)srcSize);
    FORWARD_IF_ERROR( ZSTD_compressBegin_internal(cctx,
                         dict, dictSize, ZSTD_dct_auto, ZSTD_dtlm_fast, NULL,
                         params, srcSize, ZSTDb_not_buffered) , "");
    return ZSTD_compressEnd_public(cctx, dst, dstCapacity, src, srcSize);
}

size_t ZSTD_compress_usingDict(ZSTD_CCtx* cctx,
                               void* dst, size_t dstCapacity,
                         const void* src, size_t srcSize,
                         const void* dict, size_t dictSize,
                               int compressionLevel)
{
    {
        ZSTD_parameters const params = ZSTD_getParams_internal(compressionLevel, srcSize, dict ? dictSize : 0, ZSTD_cpm_noAttachDict);
        assert(params.fParams.contentSizeFlag == 1);
        ZSTD_CCtxParams_init_internal(&cctx->simpleApiParams, &params, (compressionLevel == 0) ? ZSTD_CLEVEL_DEFAULT: compressionLevel);
    }
    DEBUGLOG(4, "ZSTD_compress_usingDict (srcSize=%u)", (unsigned)srcSize);
    return ZSTD_compress_advanced_internal(cctx, dst, dstCapacity, src, srcSize, dict, dictSize, &cctx->simpleApiParams);
}

size_t ZSTD_compressCCtx(ZSTD_CCtx* cctx,
                         void* dst, size_t dstCapacity,
                   const void* src, size_t srcSize,
                         int compressionLevel)
{
    DEBUGLOG(4, "ZSTD_compressCCtx (srcSize=%u)", (unsigned)srcSize);
    assert(cctx != NULL);
    return ZSTD_compress_usingDict(cctx, dst, dstCapacity, src, srcSize, NULL, 0, compressionLevel);
}

size_t ZSTD_compress(void* dst, size_t dstCapacity,
               const void* src, size_t srcSize,
                     int compressionLevel)
{
    size_t result;
#if ZSTD_COMPRESS_HEAPMODE
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    RETURN_ERROR_IF(!cctx, memory_allocation, "ZSTD_createCCtx failed");
    result = ZSTD_compressCCtx(cctx, dst, dstCapacity, src, srcSize, compressionLevel);
    ZSTD_freeCCtx(cctx);
#else
    ZSTD_CCtx ctxBody;
    ZSTD_initCCtx(&ctxBody, ZSTD_defaultCMem);
    result = ZSTD_compressCCtx(&ctxBody, dst, dstCapacity, src, srcSize, compressionLevel);
    ZSTD_freeCCtxContent(&ctxBody);   /* can't free ctxBody itself, as it's on stack; free only heap content */
#endif
    return result;
}


/* =====  Dictionary API  ===== */

/*! ZSTD_estimateCDictSize_advanced() :
 *  Estimate amount of memory that will be needed to create a dictionary with following arguments */
size_t ZSTD_estimateCDictSize_advanced(
        size_t dictSize, ZSTD_compressionParameters cParams,
        ZSTD_dictLoadMethod_e dictLoadMethod)
{
    DEBUGLOG(5, "sizeof(ZSTD_CDict) : %u", (unsigned)sizeof(ZSTD_CDict));
    return ZSTD_cwksp_alloc_size(sizeof(ZSTD_CDict))
         + ZSTD_cwksp_alloc_size(HUF_WORKSPACE_SIZE)
         /* enableDedicatedDictSearch == 1 ensures that CDict estimation will not be too small
          * in case we are using DDS with row-hash. */
         + ZSTD_sizeof_matchState(&cParams, ZSTD_resolveRowMatchFinderMode(ZSTD_ps_auto, &cParams),
                                  /* enableDedicatedDictSearch */ 1, /* forCCtx */ 0)
         + (dictLoadMethod == ZSTD_dlm_byRef ? 0
            : ZSTD_cwksp_alloc_size(ZSTD_cwksp_align(dictSize, sizeof(void *))));
}

size_t ZSTD_estimateCDictSize(size_t dictSize, int compressionLevel)
{
    ZSTD_compressionParameters const cParams = ZSTD_getCParams_internal(compressionLevel, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_createCDict);
    return ZSTD_estimateCDictSize_advanced(dictSize, cParams, ZSTD_dlm_byCopy);
}

size_t ZSTD_sizeof_CDict(const ZSTD_CDict* cdict)
{
    if (cdict==NULL) return 0;   /* support sizeof on NULL */
    DEBUGLOG(5, "sizeof(*cdict) : %u", (unsigned)sizeof(*cdict));
    /* cdict may be in the workspace */
    return (cdict->workspace.workspace == cdict ? 0 : sizeof(*cdict))
        + ZSTD_cwksp_sizeof(&cdict->workspace);
}

static size_t ZSTD_initCDict_internal(
                    ZSTD_CDict* cdict,
              const void* dictBuffer, size_t dictSize,
                    ZSTD_dictLoadMethod_e dictLoadMethod,
                    ZSTD_dictContentType_e dictContentType,
                    ZSTD_CCtx_params params)
{
    DEBUGLOG(3, "ZSTD_initCDict_internal (dictContentType:%u)", (unsigned)dictContentType);
    assert(!ZSTD_checkCParams(params.cParams));
    cdict->matchState.cParams = params.cParams;
    cdict->matchState.dedicatedDictSearch = params.enableDedicatedDictSearch;
    if ((dictLoadMethod == ZSTD_dlm_byRef) || (!dictBuffer) || (!dictSize)) {
        cdict->dictContent = dictBuffer;
    } else {
         void *internalBuffer = ZSTD_cwksp_reserve_object(&cdict->workspace, ZSTD_cwksp_align(dictSize, sizeof(void*)));
        RETURN_ERROR_IF(!internalBuffer, memory_allocation, "NULL pointer!");
        cdict->dictContent = internalBuffer;
        ZSTD_memcpy(internalBuffer, dictBuffer, dictSize);
    }
    cdict->dictContentSize = dictSize;
    cdict->dictContentType = dictContentType;

    cdict->entropyWorkspace = (U32*)ZSTD_cwksp_reserve_object(&cdict->workspace, HUF_WORKSPACE_SIZE);


    /* Reset the state to no dictionary */
    ZSTD_reset_compressedBlockState(&cdict->cBlockState);
    FORWARD_IF_ERROR(ZSTD_reset_matchState(
        &cdict->matchState,
        &cdict->workspace,
        &params.cParams,
        params.useRowMatchFinder,
        ZSTDcrp_makeClean,
        ZSTDirp_reset,
        ZSTD_resetTarget_CDict), "");
    /* (Maybe) load the dictionary
     * Skips loading the dictionary if it is < 8 bytes.
     */
    {   params.compressionLevel = ZSTD_CLEVEL_DEFAULT;
        params.fParams.contentSizeFlag = 1;
        {   size_t const dictID = ZSTD_compress_insertDictionary(
                    &cdict->cBlockState, &cdict->matchState, NULL, &cdict->workspace,
                    &params, cdict->dictContent, cdict->dictContentSize,
                    dictContentType, ZSTD_dtlm_full, ZSTD_tfp_forCDict, cdict->entropyWorkspace);
            FORWARD_IF_ERROR(dictID, "ZSTD_compress_insertDictionary failed");
            assert(dictID <= (size_t)(U32)-1);
            cdict->dictID = (U32)dictID;
        }
    }

    return 0;
}

static ZSTD_CDict* ZSTD_createCDict_advanced_internal(size_t dictSize,
                                      ZSTD_dictLoadMethod_e dictLoadMethod,
                                      ZSTD_compressionParameters cParams,
                                      ZSTD_paramSwitch_e useRowMatchFinder,
                                      U32 enableDedicatedDictSearch,
                                      ZSTD_customMem customMem)
{
    if ((!customMem.customAlloc) ^ (!customMem.customFree)) return NULL;

    {   size_t const workspaceSize =
            ZSTD_cwksp_alloc_size(sizeof(ZSTD_CDict)) +
            ZSTD_cwksp_alloc_size(HUF_WORKSPACE_SIZE) +
            ZSTD_sizeof_matchState(&cParams, useRowMatchFinder, enableDedicatedDictSearch, /* forCCtx */ 0) +
            (dictLoadMethod == ZSTD_dlm_byRef ? 0
             : ZSTD_cwksp_alloc_size(ZSTD_cwksp_align(dictSize, sizeof(void*))));
        void* const workspace = ZSTD_customMalloc(workspaceSize, customMem);
        ZSTD_cwksp ws;
        ZSTD_CDict* cdict;

        if (!workspace) {
            ZSTD_customFree(workspace, customMem);
            return NULL;
        }

        ZSTD_cwksp_init(&ws, workspace, workspaceSize, ZSTD_cwksp_dynamic_alloc);

        cdict = (ZSTD_CDict*)ZSTD_cwksp_reserve_object(&ws, sizeof(ZSTD_CDict));
        assert(cdict != NULL);
        ZSTD_cwksp_move(&cdict->workspace, &ws);
        cdict->customMem = customMem;
        cdict->compressionLevel = ZSTD_NO_CLEVEL; /* signals advanced API usage */
        cdict->useRowMatchFinder = useRowMatchFinder;
        return cdict;
    }
}

ZSTD_CDict* ZSTD_createCDict_advanced(const void* dictBuffer, size_t dictSize,
                                      ZSTD_dictLoadMethod_e dictLoadMethod,
                                      ZSTD_dictContentType_e dictContentType,
                                      ZSTD_compressionParameters cParams,
                                      ZSTD_customMem customMem)
{
    ZSTD_CCtx_params cctxParams;
    ZSTD_memset(&cctxParams, 0, sizeof(cctxParams));
    ZSTD_CCtxParams_init(&cctxParams, 0);
    cctxParams.cParams = cParams;
    cctxParams.customMem = customMem;
    return ZSTD_createCDict_advanced2(
        dictBuffer, dictSize,
        dictLoadMethod, dictContentType,
        &cctxParams, customMem);
}

ZSTD_CDict* ZSTD_createCDict_advanced2(
        const void* dict, size_t dictSize,
        ZSTD_dictLoadMethod_e dictLoadMethod,
        ZSTD_dictContentType_e dictContentType,
        const ZSTD_CCtx_params* originalCctxParams,
        ZSTD_customMem customMem)
{
    ZSTD_CCtx_params cctxParams = *originalCctxParams;
    ZSTD_compressionParameters cParams;
    ZSTD_CDict* cdict;

    DEBUGLOG(3, "ZSTD_createCDict_advanced2, mode %u", (unsigned)dictContentType);
    if (!customMem.customAlloc ^ !customMem.customFree) return NULL;

    if (cctxParams.enableDedicatedDictSearch) {
        cParams = ZSTD_dedicatedDictSearch_getCParams(
            cctxParams.compressionLevel, dictSize);
        ZSTD_overrideCParams(&cParams, &cctxParams.cParams);
    } else {
        cParams = ZSTD_getCParamsFromCCtxParams(
            &cctxParams, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_createCDict);
    }

    if (!ZSTD_dedicatedDictSearch_isSupported(&cParams)) {
        /* Fall back to non-DDSS params */
        cctxParams.enableDedicatedDictSearch = 0;
        cParams = ZSTD_getCParamsFromCCtxParams(
            &cctxParams, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_createCDict);
    }

    DEBUGLOG(3, "ZSTD_createCDict_advanced2: DDS: %u", cctxParams.enableDedicatedDictSearch);
    cctxParams.cParams = cParams;
    cctxParams.useRowMatchFinder = ZSTD_resolveRowMatchFinderMode(cctxParams.useRowMatchFinder, &cParams);

    cdict = ZSTD_createCDict_advanced_internal(dictSize,
                        dictLoadMethod, cctxParams.cParams,
                        cctxParams.useRowMatchFinder, cctxParams.enableDedicatedDictSearch,
                        customMem);

    if (!cdict || ZSTD_isError( ZSTD_initCDict_internal(cdict,
                                    dict, dictSize,
                                    dictLoadMethod, dictContentType,
                                    cctxParams) )) {
        ZSTD_freeCDict(cdict);
        return NULL;
    }

    return cdict;
}

ZSTD_CDict* ZSTD_createCDict(const void* dict, size_t dictSize, int compressionLevel)
{
    ZSTD_compressionParameters cParams = ZSTD_getCParams_internal(compressionLevel, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_createCDict);
    ZSTD_CDict* const cdict = ZSTD_createCDict_advanced(dict, dictSize,
                                                  ZSTD_dlm_byCopy, ZSTD_dct_auto,
                                                  cParams, ZSTD_defaultCMem);
    if (cdict)
        cdict->compressionLevel = (compressionLevel == 0) ? ZSTD_CLEVEL_DEFAULT : compressionLevel;
    return cdict;
}

ZSTD_CDict* ZSTD_createCDict_byReference(const void* dict, size_t dictSize, int compressionLevel)
{
    ZSTD_compressionParameters cParams = ZSTD_getCParams_internal(compressionLevel, ZSTD_CONTENTSIZE_UNKNOWN, dictSize, ZSTD_cpm_createCDict);
    ZSTD_CDict* const cdict = ZSTD_createCDict_advanced(dict, dictSize,
                                     ZSTD_dlm_byRef, ZSTD_dct_auto,
                                     cParams, ZSTD_defaultCMem);
    if (cdict)
        cdict->compressionLevel = (compressionLevel == 0) ? ZSTD_CLEVEL_DEFAULT : compressionLevel;
    return cdict;
}

size_t ZSTD_freeCDict(ZSTD_CDict* cdict)
{
    if (cdict==NULL) return 0;   /* support free on NULL */
    {   ZSTD_customMem const cMem = cdict->customMem;
        int cdictInWorkspace = ZSTD_cwksp_owns_buffer(&cdict->workspace, cdict);
        ZSTD_cwksp_free(&cdict->workspace, cMem);
        if (!cdictInWorkspace) {
            ZSTD_customFree(cdict, cMem);
        }
        return 0;
    }
}

/*! ZSTD_initStaticCDict_advanced() :
 *  Generate a digested dictionary in provided memory area.
 *  workspace: The memory area to emplace the dictionary into.
 *             Provided pointer must 8-bytes aligned.
 *             It must outlive dictionary usage.
 *  workspaceSize: Use ZSTD_estimateCDictSize()
 *                 to determine how large workspace must be.
 *  cParams : use ZSTD_getCParams() to transform a compression level
 *            into its relevants cParams.
 * @return : pointer to ZSTD_CDict*, or NULL if error (size too small)
 *  Note : there is no corresponding "free" function.
 *         Since workspace was allocated externally, it must be freed externally.
 */
const ZSTD_CDict* ZSTD_initStaticCDict(
                                 void* workspace, size_t workspaceSize,
                           const void* dict, size_t dictSize,
                                 ZSTD_dictLoadMethod_e dictLoadMethod,
                                 ZSTD_dictContentType_e dictContentType,
                                 ZSTD_compressionParameters cParams)
{
    ZSTD_paramSwitch_e const useRowMatchFinder = ZSTD_resolveRowMatchFinderMode(ZSTD_ps_auto, &cParams);
    /* enableDedicatedDictSearch == 1 ensures matchstate is not too small in case this CDict will be used for DDS + row hash */
    size_t const matchStateSize = ZSTD_sizeof_matchState(&cParams, useRowMatchFinder, /* enableDedicatedDictSearch */ 1, /* forCCtx */ 0);
    size_t const neededSize = ZSTD_cwksp_alloc_size(sizeof(ZSTD_CDict))
                            + (dictLoadMethod == ZSTD_dlm_byRef ? 0
                               : ZSTD_cwksp_alloc_size(ZSTD_cwksp_align(dictSize, sizeof(void*))))
                            + ZSTD_cwksp_alloc_size(HUF_WORKSPACE_SIZE)
                            + matchStateSize;
    ZSTD_CDict* cdict;
    ZSTD_CCtx_params params;

    if ((size_t)workspace & 7) return NULL;  /* 8-aligned */

    {
        ZSTD_cwksp ws;
        ZSTD_cwksp_init(&ws, workspace, workspaceSize, ZSTD_cwksp_static_alloc);
        cdict = (ZSTD_CDict*)ZSTD_cwksp_reserve_object(&ws, sizeof(ZSTD_CDict));
        if (cdict == NULL) return NULL;
        ZSTD_cwksp_move(&cdict->workspace, &ws);
    }

    DEBUGLOG(4, "(workspaceSize < neededSize) : (%u < %u) => %u",
        (unsigned)workspaceSize, (unsigned)neededSize, (unsigned)(workspaceSize < neededSize));
    if (workspaceSize < neededSize) return NULL;

    ZSTD_CCtxParams_init(&params, 0);
    params.cParams = cParams;
    params.useRowMatchFinder = useRowMatchFinder;
    cdict->useRowMatchFinder = useRowMatchFinder;
    cdict->compressionLevel = ZSTD_NO_CLEVEL;

    if (ZSTD_isError( ZSTD_initCDict_internal(cdict,
                                              dict, dictSize,
                                              dictLoadMethod, dictContentType,
                                              params) ))
        return NULL;

    return cdict;
}

ZSTD_compressionParameters ZSTD_getCParamsFromCDict(const ZSTD_CDict* cdict)
{
    assert(cdict != NULL);
    return cdict->matchState.cParams;
}

/*! ZSTD_getDictID_fromCDict() :
 *  Provides the dictID of the dictionary loaded into `cdict`.
 *  If @return == 0, the dictionary is not conformant to Zstandard specification, or empty.
 *  Non-conformant dictionaries can still be loaded, but as content-only dictionaries. */
unsigned ZSTD_getDictID_fromCDict(const ZSTD_CDict* cdict)
{
    if (cdict==NULL) return 0;
    return cdict->dictID;
}

/* ZSTD_compressBegin_usingCDict_internal() :
 * Implementation of various ZSTD_compressBegin_usingCDict* functions.
 */
static size_t ZSTD_compressBegin_usingCDict_internal(
    ZSTD_CCtx* const cctx, const ZSTD_CDict* const cdict,
    ZSTD_frameParameters const fParams, unsigned long long const pledgedSrcSize)
{
    ZSTD_CCtx_params cctxParams;
    DEBUGLOG(4, "ZSTD_compressBegin_usingCDict_internal");
    RETURN_ERROR_IF(cdict==NULL, dictionary_wrong, "NULL pointer!");
    /* Initialize the cctxParams from the cdict */
    {
        ZSTD_parameters params;
        params.fParams = fParams;
        params.cParams = ( pledgedSrcSize < ZSTD_USE_CDICT_PARAMS_SRCSIZE_CUTOFF
                        || pledgedSrcSize < cdict->dictContentSize * ZSTD_USE_CDICT_PARAMS_DICTSIZE_MULTIPLIER
                        || pledgedSrcSize == ZSTD_CONTENTSIZE_UNKNOWN
                        || cdict->compressionLevel == 0 ) ?
                ZSTD_getCParamsFromCDict(cdict)
              : ZSTD_getCParams(cdict->compressionLevel,
                                pledgedSrcSize,
                                cdict->dictContentSize);
        ZSTD_CCtxParams_init_internal(&cctxParams, &params, cdict->compressionLevel);
    }
    /* Increase window log to fit the entire dictionary and source if the
     * source size is known. Limit the increase to 19, which is the
     * window log for compression level 1 with the largest source size.
     */
    if (pledgedSrcSize != ZSTD_CONTENTSIZE_UNKNOWN) {
        U32 const limitedSrcSize = (U32)MIN(pledgedSrcSize, 1U << 19);
        U32 const limitedSrcLog = limitedSrcSize > 1 ? ZSTD_highbit32(limitedSrcSize - 1) + 1 : 1;
        cctxParams.cParams.windowLog = MAX(cctxParams.cParams.windowLog, limitedSrcLog);
    }
    return ZSTD_compressBegin_internal(cctx,
                                        NULL, 0, ZSTD_dct_auto, ZSTD_dtlm_fast,
                                        cdict,
                                        &cctxParams, pledgedSrcSize,
                                        ZSTDb_not_buffered);
}


/* ZSTD_compressBegin_usingCDict_advanced() :
 * This function is DEPRECATED.
 * cdict must be != NULL */
size_t ZSTD_compressBegin_usingCDict_advanced(
    ZSTD_CCtx* const cctx, const ZSTD_CDict* const cdict,
    ZSTD_frameParameters const fParams, unsigned long long const pledgedSrcSize)
{
    return ZSTD_compressBegin_usingCDict_internal(cctx, cdict, fParams, pledgedSrcSize);
}

/* ZSTD_compressBegin_usingCDict() :
 * cdict must be != NULL */
size_t ZSTD_compressBegin_usingCDict_deprecated(ZSTD_CCtx* cctx, const ZSTD_CDict* cdict)
{
    ZSTD_frameParameters const fParams = { 0 /*content*/, 0 /*checksum*/, 0 /*noDictID*/ };
    return ZSTD_compressBegin_usingCDict_internal(cctx, cdict, fParams, ZSTD_CONTENTSIZE_UNKNOWN);
}

size_t ZSTD_compressBegin_usingCDict(ZSTD_CCtx* cctx, const ZSTD_CDict* cdict)
{
    return ZSTD_compressBegin_usingCDict_deprecated(cctx, cdict);
}

/*! ZSTD_compress_usingCDict_internal():
 * Implementation of various ZSTD_compress_usingCDict* functions.
 */
static size_t ZSTD_compress_usingCDict_internal(ZSTD_CCtx* cctx,
                                void* dst, size_t dstCapacity,
                                const void* src, size_t srcSize,
                                const ZSTD_CDict* cdict, ZSTD_frameParameters fParams)
{
    FORWARD_IF_ERROR(ZSTD_compressBegin_usingCDict_internal(cctx, cdict, fParams, srcSize), ""); /* will check if cdict != NULL */
    return ZSTD_compressEnd_public(cctx, dst, dstCapacity, src, srcSize);
}

/*! ZSTD_compress_usingCDict_advanced():
 * This function is DEPRECATED.
 */
size_t ZSTD_compress_usingCDict_advanced(ZSTD_CCtx* cctx,
                                void* dst, size_t dstCapacity,
                                const void* src, size_t srcSize,
                                const ZSTD_CDict* cdict, ZSTD_frameParameters fParams)
{
    return ZSTD_compress_usingCDict_internal(cctx, dst, dstCapacity, src, srcSize, cdict, fParams);
}

/*! ZSTD_compress_usingCDict() :
 *  Compression using a digested Dictionary.
 *  Faster startup than ZSTD_compress_usingDict(), recommended when same dictionary is used multiple times.
 *  Note that compression parameters are decided at CDict creation time
 *  while frame parameters are hardcoded */
size_t ZSTD_compress_usingCDict(ZSTD_CCtx* cctx,
                                void* dst, size_t dstCapacity,
                                const void* src, size_t srcSize,
                                const ZSTD_CDict* cdict)
{
    ZSTD_frameParameters const fParams = { 1 /*content*/, 0 /*checksum*/, 0 /*noDictID*/ };
    return ZSTD_compress_usingCDict_internal(cctx, dst, dstCapacity, src, srcSize, cdict, fParams);
}



/* ******************************************************************
*  Streaming
********************************************************************/

ZSTD_CStream* ZSTD_createCStream(void)
{
    DEBUGLOG(3, "ZSTD_createCStream");
    return ZSTD_createCStream_advanced(ZSTD_defaultCMem);
}

ZSTD_CStream* ZSTD_initStaticCStream(void *workspace, size_t workspaceSize)
{
    return ZSTD_initStaticCCtx(workspace, workspaceSize);
}

ZSTD_CStream* ZSTD_createCStream_advanced(ZSTD_customMem customMem)
{   /* CStream and CCtx are now same object */
    return ZSTD_createCCtx_advanced(customMem);
}

size_t ZSTD_freeCStream(ZSTD_CStream* zcs)
{
    return ZSTD_freeCCtx(zcs);   /* same object */
}



/*======   Initialization   ======*/

size_t ZSTD_CStreamInSize(void)  { return ZSTD_BLOCKSIZE_MAX; }

size_t ZSTD_CStreamOutSize(void)
{
    return ZSTD_compressBound(ZSTD_BLOCKSIZE_MAX) + ZSTD_blockHeaderSize + 4 /* 32-bits hash */ ;
}

static ZSTD_cParamMode_e ZSTD_getCParamMode(ZSTD_CDict const* cdict, ZSTD_CCtx_params const* params, U64 pledgedSrcSize)
{
    if (cdict != NULL && ZSTD_shouldAttachDict(cdict, params, pledgedSrcSize))
        return ZSTD_cpm_attachDict;
    else
        return ZSTD_cpm_noAttachDict;
}

/* ZSTD_resetCStream():
 * pledgedSrcSize == 0 means "unknown" */
size_t ZSTD_resetCStream(ZSTD_CStream* zcs, unsigned long long pss)
{
    /* temporary : 0 interpreted as "unknown" during transition period.
     * Users willing to specify "unknown" **must** use ZSTD_CONTENTSIZE_UNKNOWN.
     * 0 will be interpreted as "empty" in the future.
     */
    U64 const pledgedSrcSize = (pss==0) ? ZSTD_CONTENTSIZE_UNKNOWN : pss;
    DEBUGLOG(4, "ZSTD_resetCStream: pledgedSrcSize = %u", (unsigned)pledgedSrcSize);
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setPledgedSrcSize(zcs, pledgedSrcSize) , "");
    return 0;
}

/*! ZSTD_initCStream_internal() :
 *  Note : for lib/compress only. Used by zstdmt_compress.c.
 *  Assumption 1 : params are valid
 *  Assumption 2 : either dict, or cdict, is defined, not both */
size_t ZSTD_initCStream_internal(ZSTD_CStream* zcs,
                    const void* dict, size_t dictSize, const ZSTD_CDict* cdict,
                    const ZSTD_CCtx_params* params,
                    unsigned long long pledgedSrcSize)
{
    DEBUGLOG(4, "ZSTD_initCStream_internal");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setPledgedSrcSize(zcs, pledgedSrcSize) , "");
    assert(!ZSTD_isError(ZSTD_checkCParams(params->cParams)));
    zcs->requestedParams = *params;
    assert(!((dict) && (cdict)));  /* either dict or cdict, not both */
    if (dict) {
        FORWARD_IF_ERROR( ZSTD_CCtx_loadDictionary(zcs, dict, dictSize) , "");
    } else {
        /* Dictionary is cleared if !cdict */
        FORWARD_IF_ERROR( ZSTD_CCtx_refCDict(zcs, cdict) , "");
    }
    return 0;
}

/* ZSTD_initCStream_usingCDict_advanced() :
 * same as ZSTD_initCStream_usingCDict(), with control over frame parameters */
size_t ZSTD_initCStream_usingCDict_advanced(ZSTD_CStream* zcs,
                                            const ZSTD_CDict* cdict,
                                            ZSTD_frameParameters fParams,
                                            unsigned long long pledgedSrcSize)
{
    DEBUGLOG(4, "ZSTD_initCStream_usingCDict_advanced");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setPledgedSrcSize(zcs, pledgedSrcSize) , "");
    zcs->requestedParams.fParams = fParams;
    FORWARD_IF_ERROR( ZSTD_CCtx_refCDict(zcs, cdict) , "");
    return 0;
}

/* note : cdict must outlive compression session */
size_t ZSTD_initCStream_usingCDict(ZSTD_CStream* zcs, const ZSTD_CDict* cdict)
{
    DEBUGLOG(4, "ZSTD_initCStream_usingCDict");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_refCDict(zcs, cdict) , "");
    return 0;
}


/* ZSTD_initCStream_advanced() :
 * pledgedSrcSize must be exact.
 * if srcSize is not known at init time, use value ZSTD_CONTENTSIZE_UNKNOWN.
 * dict is loaded with default parameters ZSTD_dct_auto and ZSTD_dlm_byCopy. */
size_t ZSTD_initCStream_advanced(ZSTD_CStream* zcs,
                                 const void* dict, size_t dictSize,
                                 ZSTD_parameters params, unsigned long long pss)
{
    /* for compatibility with older programs relying on this behavior.
     * Users should now specify ZSTD_CONTENTSIZE_UNKNOWN.
     * This line will be removed in the future.
     */
    U64 const pledgedSrcSize = (pss==0 && params.fParams.contentSizeFlag==0) ? ZSTD_CONTENTSIZE_UNKNOWN : pss;
    DEBUGLOG(4, "ZSTD_initCStream_advanced");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setPledgedSrcSize(zcs, pledgedSrcSize) , "");
    FORWARD_IF_ERROR( ZSTD_checkCParams(params.cParams) , "");
    ZSTD_CCtxParams_setZstdParams(&zcs->requestedParams, &params);
    FORWARD_IF_ERROR( ZSTD_CCtx_loadDictionary(zcs, dict, dictSize) , "");
    return 0;
}

size_t ZSTD_initCStream_usingDict(ZSTD_CStream* zcs, const void* dict, size_t dictSize, int compressionLevel)
{
    DEBUGLOG(4, "ZSTD_initCStream_usingDict");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setParameter(zcs, ZSTD_c_compressionLevel, compressionLevel) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_loadDictionary(zcs, dict, dictSize) , "");
    return 0;
}

size_t ZSTD_initCStream_srcSize(ZSTD_CStream* zcs, int compressionLevel, unsigned long long pss)
{
    /* temporary : 0 interpreted as "unknown" during transition period.
     * Users willing to specify "unknown" **must** use ZSTD_CONTENTSIZE_UNKNOWN.
     * 0 will be interpreted as "empty" in the future.
     */
    U64 const pledgedSrcSize = (pss==0) ? ZSTD_CONTENTSIZE_UNKNOWN : pss;
    DEBUGLOG(4, "ZSTD_initCStream_srcSize");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_refCDict(zcs, NULL) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setParameter(zcs, ZSTD_c_compressionLevel, compressionLevel) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setPledgedSrcSize(zcs, pledgedSrcSize) , "");
    return 0;
}

size_t ZSTD_initCStream(ZSTD_CStream* zcs, int compressionLevel)
{
    DEBUGLOG(4, "ZSTD_initCStream");
    FORWARD_IF_ERROR( ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_refCDict(zcs, NULL) , "");
    FORWARD_IF_ERROR( ZSTD_CCtx_setParameter(zcs, ZSTD_c_compressionLevel, compressionLevel) , "");
    return 0;
}

/*======   Compression   ======*/

static size_t ZSTD_nextInputSizeHint(const ZSTD_CCtx* cctx)
{
    if (cctx->appliedParams.inBufferMode == ZSTD_bm_stable) {
        return cctx->blockSize - cctx->stableIn_notConsumed;
    }
    assert(cctx->appliedParams.inBufferMode == ZSTD_bm_buffered);
    {   size_t hintInSize = cctx->inBuffTarget - cctx->inBuffPos;
        if (hintInSize==0) hintInSize = cctx->blockSize;
        return hintInSize;
    }
}

/** ZSTD_compressStream_generic():
 *  internal function for all *compressStream*() variants
 * @return : hint size for next input to complete ongoing block */
static size_t ZSTD_compressStream_generic(ZSTD_CStream* zcs,
                                          ZSTD_outBuffer* output,
                                          ZSTD_inBuffer* input,
                                          ZSTD_EndDirective const flushMode)
{
    const char* const istart = (assert(input != NULL), (const char*)input->src);
    const char* const iend = (istart != NULL) ? istart + input->size : istart;
    const char* ip = (istart != NULL) ? istart + input->pos : istart;
    char* const ostart = (assert(output != NULL), (char*)output->dst);
    char* const oend = (ostart != NULL) ? ostart + output->size : ostart;
    char* op = (ostart != NULL) ? ostart + output->pos : ostart;
    U32 someMoreWork = 1;

    /* check expectations */
    DEBUGLOG(5, "ZSTD_compressStream_generic, flush=%i, srcSize = %zu", (int)flushMode, input->size - input->pos);
    assert(zcs != NULL);
    if (zcs->appliedParams.inBufferMode == ZSTD_bm_stable) {
        assert(input->pos >= zcs->stableIn_notConsumed);
        input->pos -= zcs->stableIn_notConsumed;
        if (ip) ip -= zcs->stableIn_notConsumed;
        zcs->stableIn_notConsumed = 0;
    }
    if (zcs->appliedParams.inBufferMode == ZSTD_bm_buffered) {
        assert(zcs->inBuff != NULL);
        assert(zcs->inBuffSize > 0);
    }
    if (zcs->appliedParams.outBufferMode == ZSTD_bm_buffered) {
        assert(zcs->outBuff !=  NULL);
        assert(zcs->outBuffSize > 0);
    }
    if (input->src == NULL) assert(input->size == 0);
    assert(input->pos <= input->size);
    if (output->dst == NULL) assert(output->size == 0);
    assert(output->pos <= output->size);
    assert((U32)flushMode <= (U32)ZSTD_e_end);

    while (someMoreWork) {
        switch(zcs->streamStage)
        {
        case zcss_init:
            RETURN_ERROR(init_missing, "call ZSTD_initCStream() first!");

        case zcss_load:
            if ( (flushMode == ZSTD_e_end)
              && ( (size_t)(oend-op) >= ZSTD_compressBound(iend-ip)     /* Enough output space */
                || zcs->appliedParams.outBufferMode == ZSTD_bm_stable)  /* OR we are allowed to return dstSizeTooSmall */
              && (zcs->inBuffPos == 0) ) {
                /* shortcut to compression pass directly into output buffer */
                size_t const cSize = ZSTD_compressEnd_public(zcs,
                                                op, oend-op, ip, iend-ip);
                DEBUGLOG(4, "ZSTD_compressEnd : cSize=%u", (unsigned)cSize);
                FORWARD_IF_ERROR(cSize, "ZSTD_compressEnd failed");
                ip = iend;
                op += cSize;
                zcs->frameEnded = 1;
                ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
                someMoreWork = 0; break;
            }
            /* complete loading into inBuffer in buffered mode */
            if (zcs->appliedParams.inBufferMode == ZSTD_bm_buffered) {
                size_t const toLoad = zcs->inBuffTarget - zcs->inBuffPos;
                size_t const loaded = ZSTD_limitCopy(
                                        zcs->inBuff + zcs->inBuffPos, toLoad,
                                        ip, iend-ip);
                zcs->inBuffPos += loaded;
                if (ip) ip += loaded;
                if ( (flushMode == ZSTD_e_continue)
                  && (zcs->inBuffPos < zcs->inBuffTarget) ) {
                    /* not enough input to fill full block : stop here */
                    someMoreWork = 0; break;
                }
                if ( (flushMode == ZSTD_e_flush)
                  && (zcs->inBuffPos == zcs->inToCompress) ) {
                    /* empty */
                    someMoreWork = 0; break;
                }
            } else {
                assert(zcs->appliedParams.inBufferMode == ZSTD_bm_stable);
                if ( (flushMode == ZSTD_e_continue)
                  && ( (size_t)(iend - ip) < zcs->blockSize) ) {
                    /* can't compress a full block : stop here */
                    zcs->stableIn_notConsumed = (size_t)(iend - ip);
                    ip = iend;  /* pretend to have consumed input */
                    someMoreWork = 0; break;
                }
                if ( (flushMode == ZSTD_e_flush)
                  && (ip == iend) ) {
                    /* empty */
                    someMoreWork = 0; break;
                }
            }
            /* compress current block (note : this stage cannot be stopped in the middle) */
            DEBUGLOG(5, "stream compression stage (flushMode==%u)", flushMode);
            {   int const inputBuffered = (zcs->appliedParams.inBufferMode == ZSTD_bm_buffered);
                void* cDst;
                size_t cSize;
                size_t oSize = oend-op;
                size_t const iSize = inputBuffered ? zcs->inBuffPos - zcs->inToCompress
                                                   : MIN((size_t)(iend - ip), zcs->blockSize);
                if (oSize >= ZSTD_compressBound(iSize) || zcs->appliedParams.outBufferMode == ZSTD_bm_stable)
                    cDst = op;   /* compress into output buffer, to skip flush stage */
                else
                    cDst = zcs->outBuff, oSize = zcs->outBuffSize;
                if (inputBuffered) {
                    unsigned const lastBlock = (flushMode == ZSTD_e_end) && (ip==iend);
                    cSize = lastBlock ?
                            ZSTD_compressEnd_public(zcs, cDst, oSize,
                                        zcs->inBuff + zcs->inToCompress, iSize) :
                            ZSTD_compressContinue_public(zcs, cDst, oSize,
                                        zcs->inBuff + zcs->inToCompress, iSize);
                    FORWARD_IF_ERROR(cSize, "%s", lastBlock ? "ZSTD_compressEnd failed" : "ZSTD_compressContinue failed");
                    zcs->frameEnded = lastBlock;
                    /* prepare next block */
                    zcs->inBuffTarget = zcs->inBuffPos + zcs->blockSize;
                    if (zcs->inBuffTarget > zcs->inBuffSize)
                        zcs->inBuffPos = 0, zcs->inBuffTarget = zcs->blockSize;
                    DEBUGLOG(5, "inBuffTarget:%u / inBuffSize:%u",
                            (unsigned)zcs->inBuffTarget, (unsigned)zcs->inBuffSize);
                    if (!lastBlock)
                        assert(zcs->inBuffTarget <= zcs->inBuffSize);
                    zcs->inToCompress = zcs->inBuffPos;
                } else { /* !inputBuffered, hence ZSTD_bm_stable */
                    unsigned const lastBlock = (flushMode == ZSTD_e_end) && (ip + iSize == iend);
                    cSize = lastBlock ?
                            ZSTD_compressEnd_public(zcs, cDst, oSize, ip, iSize) :
                            ZSTD_compressContinue_public(zcs, cDst, oSize, ip, iSize);
                    /* Consume the input prior to error checking to mirror buffered mode. */
                    if (ip) ip += iSize;
                    FORWARD_IF_ERROR(cSize, "%s", lastBlock ? "ZSTD_compressEnd failed" : "ZSTD_compressContinue failed");
                    zcs->frameEnded = lastBlock;
                    if (lastBlock) assert(ip == iend);
                }
                if (cDst == op) {  /* no need to flush */
                    op += cSize;
                    if (zcs->frameEnded) {
                        DEBUGLOG(5, "Frame completed directly in outBuffer");
                        someMoreWork = 0;
                        ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
                    }
                    break;
                }
                zcs->outBuffContentSize = cSize;
                zcs->outBuffFlushedSize = 0;
                zcs->streamStage = zcss_flush; /* pass-through to flush stage */
            }
	    ZSTD_FALLTHROUGH;
        case zcss_flush:
            DEBUGLOG(5, "flush stage");
            assert(zcs->appliedParams.outBufferMode == ZSTD_bm_buffered);
            {   size_t const toFlush = zcs->outBuffContentSize - zcs->outBuffFlushedSize;
                size_t const flushed = ZSTD_limitCopy(op, (size_t)(oend-op),
                            zcs->outBuff + zcs->outBuffFlushedSize, toFlush);
                DEBUGLOG(5, "toFlush: %u into %u ==> flushed: %u",
                            (unsigned)toFlush, (unsigned)(oend-op), (unsigned)flushed);
                if (flushed)
                    op += flushed;
                zcs->outBuffFlushedSize += flushed;
                if (toFlush!=flushed) {
                    /* flush not fully completed, presumably because dst is too small */
                    assert(op==oend);
                    someMoreWork = 0;
                    break;
                }
                zcs->outBuffContentSize = zcs->outBuffFlushedSize = 0;
                if (zcs->frameEnded) {
                    DEBUGLOG(5, "Frame completed on flush");
                    someMoreWork = 0;
                    ZSTD_CCtx_reset(zcs, ZSTD_reset_session_only);
                    break;
                }
                zcs->streamStage = zcss_load;
                break;
            }

        default: /* impossible */
            assert(0);
        }
    }

    input->pos = ip - istart;
    output->pos = op - ostart;
    if (zcs->frameEnded) return 0;
    return ZSTD_nextInputSizeHint(zcs);
}

static size_t ZSTD_nextInputSizeHint_MTorST(const ZSTD_CCtx* cctx)
{
#ifdef ZSTD_MULTITHREAD
    if (cctx->appliedParams.nbWorkers >= 1) {
        assert(cctx->mtctx != NULL);
        return ZSTDMT_nextInputSizeHint(cctx->mtctx);
    }
#endif
    return ZSTD_nextInputSizeHint(cctx);

}

size_t ZSTD_compressStream(ZSTD_CStream* zcs, ZSTD_outBuffer* output, ZSTD_inBuffer* input)
{
    FORWARD_IF_ERROR( ZSTD_compressStream2(zcs, output, input, ZSTD_e_continue) , "");
    return ZSTD_nextInputSizeHint_MTorST(zcs);
}

/* After a compression call set the expected input/output buffer.
 * This is validated at the start of the next compression call.
 */
static void
ZSTD_setBufferExpectations(ZSTD_CCtx* cctx, const ZSTD_outBuffer* output, const ZSTD_inBuffer* input)
{
    DEBUGLOG(5, "ZSTD_setBufferExpectations (for advanced stable in/out modes)");
    if (cctx->appliedParams.inBufferMode == ZSTD_bm_stable) {
        cctx->expectedInBuffer = *input;
    }
    if (cctx->appliedParams.outBufferMode == ZSTD_bm_stable) {
        cctx->expectedOutBufferSize = output->size - output->pos;
    }
}

/* Validate that the input/output buffers match the expectations set by
 * ZSTD_setBufferExpectations.
 */
static size_t ZSTD_checkBufferStability(ZSTD_CCtx const* cctx,
                                        ZSTD_outBuffer const* output,
                                        ZSTD_inBuffer const* input,
                                        ZSTD_EndDirective endOp)
{
    if (cctx->appliedParams.inBufferMode == ZSTD_bm_stable) {
        ZSTD_inBuffer const expect = cctx->expectedInBuffer;
        if (expect.src != input->src || expect.pos != input->pos)
            RETURN_ERROR(stabilityCondition_notRespected, "ZSTD_c_stableInBuffer enabled but input differs!");
    }
    (void)endOp;
    if (cctx->appliedParams.outBufferMode == ZSTD_bm_stable) {
        size_t const outBufferSize = output->size - output->pos;
        if (cctx->expectedOutBufferSize != outBufferSize)
            RETURN_ERROR(stabilityCondition_notRespected, "ZSTD_c_stableOutBuffer enabled but output size differs!");
    }
    return 0;
}

static size_t ZSTD_CCtx_init_compressStream2(ZSTD_CCtx* cctx,
                                             ZSTD_EndDirective endOp,
                                             size_t inSize)
{
    ZSTD_CCtx_params params = cctx->requestedParams;
    ZSTD_prefixDict const prefixDict = cctx->prefixDict;
    FORWARD_IF_ERROR( ZSTD_initLocalDict(cctx) , ""); /* Init the local dict if present. */
    ZSTD_memset(&cctx->prefixDict, 0, sizeof(cctx->prefixDict));   /* single usage */
    assert(prefixDict.dict==NULL || cctx->cdict==NULL);    /* only one can be set */
    if (cctx->cdict && !cctx->localDict.cdict) {
        /* Let the cdict's compression level take priority over the requested params.
         * But do not take the cdict's compression level if the "cdict" is actually a localDict
         * generated from ZSTD_initLocalDict().
         */
        params.compressionLevel = cctx->cdict->compressionLevel;
    }
    DEBUGLOG(4, "ZSTD_compressStream2 : transparent init stage");
    if (endOp == ZSTD_e_end) cctx->pledgedSrcSizePlusOne = inSize + 1;  /* auto-determine pledgedSrcSize */

    {   size_t const dictSize = prefixDict.dict
                ? prefixDict.dictSize
                : (cctx->cdict ? cctx->cdict->dictContentSize : 0);
        ZSTD_cParamMode_e const mode = ZSTD_getCParamMode(cctx->cdict, &params, cctx->pledgedSrcSizePlusOne - 1);
        params.cParams = ZSTD_getCParamsFromCCtxParams(
                &params, cctx->pledgedSrcSizePlusOne-1,
                dictSize, mode);
    }

    params.useBlockSplitter = ZSTD_resolveBlockSplitterMode(params.useBlockSplitter, &params.cParams);
    params.ldmParams.enableLdm = ZSTD_resolveEnableLdm(params.ldmParams.enableLdm, &params.cParams);
    params.useRowMatchFinder = ZSTD_resolveRowMatchFinderMode(params.useRowMatchFinder, &params.cParams);
    params.validateSequences = ZSTD_resolveExternalSequenceValidation(params.validateSequences);
    params.maxBlockSize = ZSTD_resolveMaxBlockSize(params.maxBlockSize);
    params.searchForExternalRepcodes = ZSTD_resolveExternalRepcodeSearch(params.searchForExternalRepcodes, params.compressionLevel);

#ifdef ZSTD_MULTITHREAD
    /* If external matchfinder is enabled, make sure to fail before checking job size (for consistency) */
    RETURN_ERROR_IF(
        ZSTD_hasExtSeqProd(&params) && params.nbWorkers >= 1,
        parameter_combination_unsupported,
        "External sequence producer isn't supported with nbWorkers >= 1"
    );

    if ((cctx->pledgedSrcSizePlusOne-1) <= ZSTDMT_JOBSIZE_MIN) {
        params.nbWorkers = 0; /* do not invoke multi-threading when src size is too small */
    }
    if (params.nbWorkers > 0) {
#if ZSTD_TRACE
        cctx->traceCtx = (ZSTD_trace_compress_begin != NULL) ? ZSTD_trace_compress_begin(cctx) : 0;
#endif
        /* mt context creation */
        if (cctx->mtctx == NULL) {
            DEBUGLOG(4, "ZSTD_compressStream2: creating new mtctx for nbWorkers=%u",
                        params.nbWorkers);
            cctx->mtctx = ZSTDMT_createCCtx_advanced((U32)params.nbWorkers, cctx->customMem, cctx->pool);
            RETURN_ERROR_IF(cctx->mtctx == NULL, memory_allocation, "NULL pointer!");
        }
        /* mt compression */
        DEBUGLOG(4, "call ZSTDMT_initCStream_internal as nbWorkers=%u", params.nbWorkers);
        FORWARD_IF_ERROR( ZSTDMT_initCStream_internal(
                    cctx->mtctx,
                    prefixDict.dict, prefixDict.dictSize, prefixDict.dictContentType,
                    cctx->cdict, params, cctx->pledgedSrcSizePlusOne-1) , "");
        cctx->dictID = cctx->cdict ? cctx->cdict->dictID : 0;
        cctx->dictContentSize = cctx->cdict ? cctx->cdict->dictContentSize : prefixDict.dictSize;
        cctx->consumedSrcSize = 0;
        cctx->producedCSize = 0;
        cctx->streamStage = zcss_load;
        cctx->appliedParams = params;
    } else
#endif  /* ZSTD_MULTITHREAD */
    {   U64 const pledgedSrcSize = cctx->pledgedSrcSizePlusOne - 1;
        assert(!ZSTD_isError(ZSTD_checkCParams(params.cParams)));
        FORWARD_IF_ERROR( ZSTD_compressBegin_internal(cctx,
                prefixDict.dict, prefixDict.dictSize, prefixDict.dictContentType, ZSTD_dtlm_fast,
                cctx->cdict,
                &params, pledgedSrcSize,
                ZSTDb_buffered) , "");
        assert(cctx->appliedParams.nbWorkers == 0);
        cctx->inToCompress = 0;
        cctx->inBuffPos = 0;
        if (cctx->appliedParams.inBufferMode == ZSTD_bm_buffered) {
            /* for small input: avoid automatic flush on reaching end of block, since
            * it would require to add a 3-bytes null block to end frame
            */
            cctx->inBuffTarget = cctx->blockSize + (cctx->blockSize == pledgedSrcSize);
        } else {
            cctx->inBuffTarget = 0;
        }
        cctx->outBuffContentSize = cctx->outBuffFlushedSize = 0;
        cctx->streamStage = zcss_load;
        cctx->frameEnded = 0;
    }
    return 0;
}

/* @return provides a minimum amount of data remaining to be flushed from internal buffers
 */
size_t ZSTD_compressStream2( ZSTD_CCtx* cctx,
                             ZSTD_outBuffer* output,
                             ZSTD_inBuffer* input,
                             ZSTD_EndDirective endOp)
{
    DEBUGLOG(5, "ZSTD_compressStream2, endOp=%u ", (unsigned)endOp);
    /* check conditions */
    RETURN_ERROR_IF(output->pos > output->size, dstSize_tooSmall, "invalid output buffer");
    RETURN_ERROR_IF(input->pos  > input->size, srcSize_wrong, "invalid input buffer");
    RETURN_ERROR_IF((U32)endOp > (U32)ZSTD_e_end, parameter_outOfBound, "invalid endDirective");
    assert(cctx != NULL);

    /* transparent initialization stage */
    if (cctx->streamStage == zcss_init) {
        size_t const inputSize = input->size - input->pos;  /* no obligation to start from pos==0 */
        size_t const totalInputSize = inputSize + cctx->stableIn_notConsumed;
        if ( (cctx->requestedParams.inBufferMode == ZSTD_bm_stable) /* input is presumed stable, across invocations */
          && (endOp == ZSTD_e_continue)                             /* no flush requested, more input to come */
          && (totalInputSize < ZSTD_BLOCKSIZE_MAX) ) {              /* not even reached one block yet */
            if (cctx->stableIn_notConsumed) {  /* not the first time */
                /* check stable source guarantees */
                RETURN_ERROR_IF(input->src != cctx->expectedInBuffer.src, stabilityCondition_notRespected, "stableInBuffer condition not respected: wrong src pointer");
                RETURN_ERROR_IF(input->pos != cctx->expectedInBuffer.size, stabilityCondition_notRespected, "stableInBuffer condition not respected: externally modified pos");
            }
            /* pretend input was consumed, to give a sense forward progress */
            input->pos = input->size;
            /* save stable inBuffer, for later control, and flush/end */
            cctx->expectedInBuffer = *input;
            /* but actually input wasn't consumed, so keep track of position from where compression shall resume */
            cctx->stableIn_notConsumed += inputSize;
            /* don't initialize yet, wait for the first block of flush() order, for better parameters adaptation */
            return ZSTD_FRAMEHEADERSIZE_MIN(cctx->requestedParams.format);  /* at least some header to produce */
        }
        FORWARD_IF_ERROR(ZSTD_CCtx_init_compressStream2(cctx, endOp, totalInputSize), "compressStream2 initialization failed");
        ZSTD_setBufferExpectations(cctx, output, input);   /* Set initial buffer expectations now that we've initialized */
    }
    /* end of transparent initialization stage */

    FORWARD_IF_ERROR(ZSTD_checkBufferStability(cctx, output, input, endOp), "invalid buffers");
    /* compression stage */
#ifdef ZSTD_MULTITHREAD
    if (cctx->appliedParams.nbWorkers > 0) {
        size_t flushMin;
        if (cctx->cParamsChanged) {
            ZSTDMT_updateCParams_whileCompressing(cctx->mtctx, &cctx->requestedParams);
            cctx->cParamsChanged = 0;
        }
        if (cctx->stableIn_notConsumed) {
            assert(cctx->appliedParams.inBufferMode == ZSTD_bm_stable);
            /* some early data was skipped - make it available for consumption */
            assert(input->pos >= cctx->stableIn_notConsumed);
            input->pos -= cctx->stableIn_notConsumed;
            cctx->stableIn_notConsumed = 0;
        }
        for (;;) {
            size_t const ipos = input->pos;
            size_t const opos = output->pos;
            flushMin = ZSTDMT_compressStream_generic(cctx->mtctx, output, input, endOp);
            cctx->consumedSrcSize += (U64)(input->pos - ipos);
            cctx->producedCSize += (U64)(output->pos - opos);
            if ( ZSTD_isError(flushMin)
              || (endOp == ZSTD_e_end && flushMin == 0) ) { /* compression completed */
                if (flushMin == 0)
                    ZSTD_CCtx_trace(cctx, 0);
                ZSTD_CCtx_reset(cctx, ZSTD_reset_session_only);
            }
            FORWARD_IF_ERROR(flushMin, "ZSTDMT_compressStream_generic failed");

            if (endOp == ZSTD_e_continue) {
                /* We only require some progress with ZSTD_e_continue, not maximal progress.
                 * We're done if we've consumed or produced any bytes, or either buffer is
                 * full.
                 */
                if (input->pos != ipos || output->pos != opos || input->pos == input->size || output->pos == output->size)
                    break;
            } else {
                assert(endOp == ZSTD_e_flush || endOp == ZSTD_e_end);
                /* We require maximal progress. We're done when the flush is complete or the
                 * output buffer is full.
                 */
                if (flushMin == 0 || output->pos == output->size)
                    break;
            }
        }
        DEBUGLOG(5, "completed ZSTD_compressStream2 delegating to ZSTDMT_compressStream_generic");
        /* Either we don't require maximum forward progress, we've finished the
         * flush, or we are out of output space.
         */
        assert(endOp == ZSTD_e_continue || flushMin == 0 || output->pos == output->size);
        ZSTD_setBufferExpectations(cctx, output, input);
        return flushMin;
    }
#endif /* ZSTD_MULTITHREAD */
    FORWARD_IF_ERROR( ZSTD_compressStream_generic(cctx, output, input, endOp) , "");
    DEBUGLOG(5, "completed ZSTD_compressStream2");
    ZSTD_setBufferExpectations(cctx, output, input);
    return cctx->outBuffContentSize - cctx->outBuffFlushedSize; /* remaining to flush */
}

size_t ZSTD_compressStream2_simpleArgs (
                            ZSTD_CCtx* cctx,
                            void* dst, size_t dstCapacity, size_t* dstPos,
                      const void* src, size_t srcSize, size_t* srcPos,
                            ZSTD_EndDirective endOp)
{
    ZSTD_outBuffer output;
    ZSTD_inBuffer  input;
    output.dst = dst;
    output.size = dstCapacity;
    output.pos = *dstPos;
    input.src = src;
    input.size = srcSize;
    input.pos = *srcPos;
    /* ZSTD_compressStream2() will check validity of dstPos and srcPos */
    {   size_t const cErr = ZSTD_compressStream2(cctx, &output, &input, endOp);
        *dstPos = output.pos;
        *srcPos = input.pos;
        return cErr;
    }
}

size_t ZSTD_compress2(ZSTD_CCtx* cctx,
                      void* dst, size_t dstCapacity,
                      const void* src, size_t srcSize)
{
    ZSTD_bufferMode_e const originalInBufferMode = cctx->requestedParams.inBufferMode;
    ZSTD_bufferMode_e const originalOutBufferMode = cctx->requestedParams.outBufferMode;
    DEBUGLOG(4, "ZSTD_compress2 (srcSize=%u)", (unsigned)srcSize);
    ZSTD_CCtx_reset(cctx, ZSTD_reset_session_only);
    /* Enable stable input/output buffers. */
    cctx->requestedParams.inBufferMode = ZSTD_bm_stable;
    cctx->requestedParams.outBufferMode = ZSTD_bm_stable;
    {   size_t oPos = 0;
        size_t iPos = 0;
        size_t const result = ZSTD_compressStream2_simpleArgs(cctx,
                                        dst, dstCapacity, &oPos,
                                        src, srcSize, &iPos,
                                        ZSTD_e_end);
        /* Reset to the original values. */
        cctx->requestedParams.inBufferMode = originalInBufferMode;
        cctx->requestedParams.outBufferMode = originalOutBufferMode;

        FORWARD_IF_ERROR(result, "ZSTD_compressStream2_simpleArgs failed");
        if (result != 0) {  /* compression not completed, due to lack of output space */
            assert(oPos == dstCapacity);
            RETURN_ERROR(dstSize_tooSmall, "");
        }
        assert(iPos == srcSize);   /* all input is expected consumed */
        return oPos;
    }
}

/* ZSTD_validateSequence() :
 * @offCode : is presumed to follow format required by ZSTD_storeSeq()
 * @returns a ZSTD error code if sequence is not valid
 */
static size_t
ZSTD_validateSequence(U32 offCode, U32 matchLength, U32 minMatch,
                      size_t posInSrc, U32 windowLog, size_t dictSize, int useSequenceProducer)
{
    U32 const windowSize = 1u << windowLog;
    /* posInSrc represents the amount of data the decoder would decode up to this point.
     * As long as the amount of data decoded is less than or equal to window size, offsets may be
     * larger than the total length of output decoded in order to reference the dict, even larger than
     * window size. After output surpasses windowSize, we're limited to windowSize offsets again.
     */
    size_t const offsetBound = posInSrc > windowSize ? (size_t)windowSize : posInSrc + (size_t)dictSize;
    size_t const matchLenLowerBound = (minMatch == 3 || useSequenceProducer) ? 3 : 4;
    RETURN_ERROR_IF(offCode > OFFSET_TO_OFFBASE(offsetBound), externalSequences_invalid, "Offset too large!");
    /* Validate maxNbSeq is large enough for the given matchLength and minMatch */
    RETURN_ERROR_IF(matchLength < matchLenLowerBound, externalSequences_invalid, "Matchlength too small for the minMatch");
    return 0;
}

/* Returns an offset code, given a sequence's raw offset, the ongoing repcode array, and whether litLength == 0 */
static U32 ZSTD_finalizeOffBase(U32 rawOffset, const U32 rep[ZSTD_REP_NUM], U32 ll0)
{
    U32 offBase = OFFSET_TO_OFFBASE(rawOffset);

    if (!ll0 && rawOffset == rep[0]) {
        offBase = REPCODE1_TO_OFFBASE;
    } else if (rawOffset == rep[1]) {
        offBase = REPCODE_TO_OFFBASE(2 - ll0);
    } else if (rawOffset == rep[2]) {
        offBase = REPCODE_TO_OFFBASE(3 - ll0);
    } else if (ll0 && rawOffset == rep[0] - 1) {
        offBase = REPCODE3_TO_OFFBASE;
    }
    return offBase;
}

size_t
ZSTD_copySequencesToSeqStoreExplicitBlockDelim(ZSTD_CCtx* cctx,
                                              ZSTD_sequencePosition* seqPos,
                                        const ZSTD_Sequence* const inSeqs, size_t inSeqsSize,
                                        const void* src, size_t blockSize,
                                        ZSTD_paramSwitch_e externalRepSearch)
{
    U32 idx = seqPos->idx;
    U32 const startIdx = idx;
    BYTE const* ip = (BYTE const*)(src);
    const BYTE* const iend = ip + blockSize;
    repcodes_t updatedRepcodes;
    U32 dictSize;

    DEBUGLOG(5, "ZSTD_copySequencesToSeqStoreExplicitBlockDelim (blockSize = %zu)", blockSize);

    if (cctx->cdict) {
        dictSize = (U32)cctx->cdict->dictContentSize;
    } else if (cctx->prefixDict.dict) {
        dictSize = (U32)cctx->prefixDict.dictSize;
    } else {
        dictSize = 0;
    }
    ZSTD_memcpy(updatedRepcodes.rep, cctx->blockState.prevCBlock->rep, sizeof(repcodes_t));
    for (; idx < inSeqsSize && (inSeqs[idx].matchLength != 0 || inSeqs[idx].offset != 0); ++idx) {
        U32 const litLength = inSeqs[idx].litLength;
        U32 const matchLength = inSeqs[idx].matchLength;
        U32 offBase;

        if (externalRepSearch == ZSTD_ps_disable) {
            offBase = OFFSET_TO_OFFBASE(inSeqs[idx].offset);
        } else {
            U32 const ll0 = (litLength == 0);
            offBase = ZSTD_finalizeOffBase(inSeqs[idx].offset, updatedRepcodes.rep, ll0);
            ZSTD_updateRep(updatedRepcodes.rep, offBase, ll0);
        }

        DEBUGLOG(6, "Storing sequence: (of: %u, ml: %u, ll: %u)", offBase, matchLength, litLength);
        if (cctx->appliedParams.validateSequences) {
            seqPos->posInSrc += litLength + matchLength;
            FORWARD_IF_ERROR(ZSTD_validateSequence(offBase, matchLength, cctx->appliedParams.cParams.minMatch, seqPos->posInSrc,
                                                cctx->appliedParams.cParams.windowLog, dictSize, ZSTD_hasExtSeqProd(&cctx->appliedParams)),
                                                "Sequence validation failed");
        }
        RETURN_ERROR_IF(idx - seqPos->idx >= cctx->seqStore.maxNbSeq, externalSequences_invalid,
                        "Not enough memory allocated. Try adjusting ZSTD_c_minMatch.");
        ZSTD_storeSeq(&cctx->seqStore, litLength, ip, iend, offBase, matchLength);
        ip += matchLength + litLength;
    }

    /* If we skipped repcode search while parsing, we need to update repcodes now */
    assert(externalRepSearch != ZSTD_ps_auto);
    assert(idx >= startIdx);
    if (externalRepSearch == ZSTD_ps_disable && idx != startIdx) {
        U32* const rep = updatedRepcodes.rep;
        U32 lastSeqIdx = idx - 1; /* index of last non-block-delimiter sequence */

        if (lastSeqIdx >= startIdx + 2) {
            rep[2] = inSeqs[lastSeqIdx - 2].offset;
            rep[1] = inSeqs[lastSeqIdx - 1].offset;
            rep[0] = inSeqs[lastSeqIdx].offset;
        } else if (lastSeqIdx == startIdx + 1) {
            rep[2] = rep[0];
            rep[1] = inSeqs[lastSeqIdx - 1].offset;
            rep[0] = inSeqs[lastSeqIdx].offset;
        } else {
            assert(lastSeqIdx == startIdx);
            rep[2] = rep[1];
            rep[1] = rep[0];
            rep[0] = inSeqs[lastSeqIdx].offset;
        }
    }

    ZSTD_memcpy(cctx->blockState.nextCBlock->rep, updatedRepcodes.rep, sizeof(repcodes_t));

    if (inSeqs[idx].litLength) {
        DEBUGLOG(6, "Storing last literals of size: %u", inSeqs[idx].litLength);
        ZSTD_storeLastLiterals(&cctx->seqStore, ip, inSeqs[idx].litLength);
        ip += inSeqs[idx].litLength;
        seqPos->posInSrc += inSeqs[idx].litLength;
    }
    RETURN_ERROR_IF(ip != iend, externalSequences_invalid, "Blocksize doesn't agree with block delimiter!");
    seqPos->idx = idx+1;
    return 0;
}

size_t
ZSTD_copySequencesToSeqStoreNoBlockDelim(ZSTD_CCtx* cctx, ZSTD_sequencePosition* seqPos,
                                   const ZSTD_Sequence* const inSeqs, size_t inSeqsSize,
                                   const void* src, size_t blockSize, ZSTD_paramSwitch_e externalRepSearch)
{
    U32 idx = seqPos->idx;
    U32 startPosInSequence = seqPos->posInSequence;
    U32 endPosInSequence = seqPos->posInSequence + (U32)blockSize;
    size_t dictSize;
    BYTE const* ip = (BYTE const*)(src);
    BYTE const* iend = ip + blockSize;  /* May be adjusted if we decide to process fewer than blockSize bytes */
    repcodes_t updatedRepcodes;
    U32 bytesAdjustment = 0;
    U32 finalMatchSplit = 0;

    /* TODO(embg) support fast parsing mode in noBlockDelim mode */
    (void)externalRepSearch;

    if (cctx->cdict) {
        dictSize = cctx->cdict->dictContentSize;
    } else if (cctx->prefixDict.dict) {
        dictSize = cctx->prefixDict.dictSize;
    } else {
        dictSize = 0;
    }
    DEBUGLOG(5, "ZSTD_copySequencesToSeqStoreNoBlockDelim: idx: %u PIS: %u blockSize: %zu", idx, startPosInSequence, blockSize);
    DEBUGLOG(5, "Start seq: idx: %u (of: %u ml: %u ll: %u)", idx, inSeqs[idx].offset, inSeqs[idx].matchLength, inSeqs[idx].litLength);
    ZSTD_memcpy(updatedRepcodes.rep, cctx->blockState.prevCBlock->rep, sizeof(repcodes_t));
    while (endPosInSequence && idx < inSeqsSize && !finalMatchSplit) {
        const ZSTD_Sequence currSeq = inSeqs[idx];
        U32 litLength = currSeq.litLength;
        U32 matchLength = currSeq.matchLength;
        U32 const rawOffset = currSeq.offset;
        U32 offBase;

        /* Modify the sequence depending on where endPosInSequence lies */
        if (endPosInSequence >= currSeq.litLength + currSeq.matchLength) {
            if (startPosInSequence >= litLength) {
                startPosInSequence -= litLength;
                litLength = 0;
                matchLength -= startPosInSequence;
            } else {
                litLength -= startPosInSequence;
            }
            /* Move to the next sequence */
            endPosInSequence -= currSeq.litLength + currSeq.matchLength;
            startPosInSequence = 0;
        } else {
            /* This is the final (partial) sequence we're adding from inSeqs, and endPosInSequence
               does not reach the end of the match. So, we have to split the sequence */
            DEBUGLOG(6, "Require a split: diff: %u, idx: %u PIS: %u",
                     currSeq.litLength + currSeq.matchLength - endPosInSequence, idx, endPosInSequence);
            if (endPosInSequence > litLength) {
                U32 firstHalfMatchLength;
                litLength = startPosInSequence >= litLength ? 0 : litLength - startPosInSequence;
                firstHalfMatchLength = endPosInSequence - startPosInSequence - litLength;
                if (matchLength > blockSize && firstHalfMatchLength >= cctx->appliedParams.cParams.minMatch) {
                    /* Only ever split the match if it is larger than the block size */
                    U32 secondHalfMatchLength = currSeq.matchLength + currSeq.litLength - endPosInSequence;
                    if (secondHalfMatchLength < cctx->appliedParams.cParams.minMatch) {
                        /* Move the endPosInSequence backward so that it creates match of minMatch length */
                        endPosInSequence -= cctx->appliedParams.cParams.minMatch - secondHalfMatchLength;
                        bytesAdjustment = cctx->appliedParams.cParams.minMatch - secondHalfMatchLength;
                        firstHalfMatchLength -= bytesAdjustment;
                    }
                    matchLength = firstHalfMatchLength;
                    /* Flag that we split the last match - after storing the sequence, exit the loop,
                       but keep the value of endPosInSequence */
                    finalMatchSplit = 1;
                } else {
                    /* Move the position in sequence backwards so that we don't split match, and break to store
                     * the last literals. We use the original currSeq.litLength as a marker for where endPosInSequence
                     * should go. We prefer to do this whenever it is not necessary to split the match, or if doing so
                     * would cause the first half of the match to be too small
                     */
                    bytesAdjustment = endPosInSequence - currSeq.litLength;
                    endPosInSequence = currSeq.litLength;
                    break;
                }
            } else {
                /* This sequence ends inside the literals, break to store the last literals */
                break;
            }
        }
        /* Check if this offset can be represented with a repcode */
        {   U32 const ll0 = (litLength == 0);
            offBase = ZSTD_finalizeOffBase(rawOffset, updatedRepcodes.rep, ll0);
            ZSTD_updateRep(updatedRepcodes.rep, offBase, ll0);
        }

        if (cctx->appliedParams.validateSequences) {
            seqPos->posInSrc += litLength + matchLength;
            FORWARD_IF_ERROR(ZSTD_validateSequence(offBase, matchLength, cctx->appliedParams.cParams.minMatch, seqPos->posInSrc,
                                                   cctx->appliedParams.cParams.windowLog, dictSize, ZSTD_hasExtSeqProd(&cctx->appliedParams)),
                                                   "Sequence validation failed");
        }
        DEBUGLOG(6, "Storing sequence: (of: %u, ml: %u, ll: %u)", offBase, matchLength, litLength);
        RETURN_ERROR_IF(idx - seqPos->idx >= cctx->seqStore.maxNbSeq, externalSequences_invalid,
                        "Not enough memory allocated. Try adjusting ZSTD_c_minMatch.");
        ZSTD_storeSeq(&cctx->seqStore, litLength, ip, iend, offBase, matchLength);
        ip += matchLength + litLength;
        if (!finalMatchSplit)
            idx++; /* Next Sequence */
    }
    DEBUGLOG(5, "Ending seq: idx: %u (of: %u ml: %u ll: %u)", idx, inSeqs[idx].offset, inSeqs[idx].matchLength, inSeqs[idx].litLength);
    assert(idx == inSeqsSize || endPosInSequence <= inSeqs[idx].litLength + inSeqs[idx].matchLength);
    seqPos->idx = idx;
    seqPos->posInSequence = endPosInSequence;
    ZSTD_memcpy(cctx->blockState.nextCBlock->rep, updatedRepcodes.rep, sizeof(repcodes_t));

    iend -= bytesAdjustment;
    if (ip != iend) {
        /* Store any last literals */
        U32 lastLLSize = (U32)(iend - ip);
        assert(ip <= iend);
        DEBUGLOG(6, "Storing last literals of size: %u", lastLLSize);
        ZSTD_storeLastLiterals(&cctx->seqStore, ip, lastLLSize);
        seqPos->posInSrc += lastLLSize;
    }

    return bytesAdjustment;
}

typedef size_t (*ZSTD_sequenceCopier) (ZSTD_CCtx* cctx, ZSTD_sequencePosition* seqPos,
                                       const ZSTD_Sequence* const inSeqs, size_t inSeqsSize,
                                       const void* src, size_t blockSize, ZSTD_paramSwitch_e externalRepSearch);
static ZSTD_sequenceCopier ZSTD_selectSequenceCopier(ZSTD_sequenceFormat_e mode)
{
    ZSTD_sequenceCopier sequenceCopier = NULL;
    assert(ZSTD_cParam_withinBounds(ZSTD_c_blockDelimiters, mode));
    if (mode == ZSTD_sf_explicitBlockDelimiters) {
        return ZSTD_copySequencesToSeqStoreExplicitBlockDelim;
    } else if (mode == ZSTD_sf_noBlockDelimiters) {
        return ZSTD_copySequencesToSeqStoreNoBlockDelim;
    }
    assert(sequenceCopier != NULL);
    return sequenceCopier;
}

/* Discover the size of next block by searching for the delimiter.
 * Note that a block delimiter **must** exist in this mode,
 * otherwise it's an input error.
 * The block size retrieved will be later compared to ensure it remains within bounds */
static size_t
blockSize_explicitDelimiter(const ZSTD_Sequence* inSeqs, size_t inSeqsSize, ZSTD_sequencePosition seqPos)
{
    int end = 0;
    size_t blockSize = 0;
    size_t spos = seqPos.idx;
    DEBUGLOG(6, "blockSize_explicitDelimiter : seq %zu / %zu", spos, inSeqsSize);
    assert(spos <= inSeqsSize);
    while (spos < inSeqsSize) {
        end = (inSeqs[spos].offset == 0);
        blockSize += inSeqs[spos].litLength + inSeqs[spos].matchLength;
        if (end) {
            if (inSeqs[spos].matchLength != 0)
                RETURN_ERROR(externalSequences_invalid, "delimiter format error : both matchlength and offset must be == 0");
            break;
        }
        spos++;
    }
    if (!end)
        RETURN_ERROR(externalSequences_invalid, "Reached end of sequences without finding a block delimiter");
    return blockSize;
}

/* More a "target" block size */
static size_t blockSize_noDelimiter(size_t blockSize, size_t remaining)
{
    int const lastBlock = (remaining <= blockSize);
    return lastBlock ? remaining : blockSize;
}

static size_t determine_blockSize(ZSTD_sequenceFormat_e mode,
                           size_t blockSize, size_t remaining,
                     const ZSTD_Sequence* inSeqs, size_t inSeqsSize, ZSTD_sequencePosition seqPos)
{
    DEBUGLOG(6, "determine_blockSize : remainingSize = %zu", remaining);
    if (mode == ZSTD_sf_noBlockDelimiters)
        return blockSize_noDelimiter(blockSize, remaining);
    {   size_t const explicitBlockSize = blockSize_explicitDelimiter(inSeqs, inSeqsSize, seqPos);
        FORWARD_IF_ERROR(explicitBlockSize, "Error while determining block size with explicit delimiters");
        if (explicitBlockSize > blockSize)
            RETURN_ERROR(externalSequences_invalid, "sequences incorrectly define a too large block");
        if (explicitBlockSize > remaining)
            RETURN_ERROR(externalSequences_invalid, "sequences define a frame longer than source");
        return explicitBlockSize;
    }
}

/* Compress, block-by-block, all of the sequences given.
 *
 * Returns the cumulative size of all compressed blocks (including their headers),
 * otherwise a ZSTD error.
 */
static size_t
ZSTD_compressSequences_internal(ZSTD_CCtx* cctx,
                                void* dst, size_t dstCapacity,
                          const ZSTD_Sequence* inSeqs, size_t inSeqsSize,
                          const void* src, size_t srcSize)
{
    size_t cSize = 0;
    size_t remaining = srcSize;
    ZSTD_sequencePosition seqPos = {0, 0, 0};

    BYTE const* ip = (BYTE const*)src;
    BYTE* op = (BYTE*)dst;
    ZSTD_sequenceCopier const sequenceCopier = ZSTD_selectSequenceCopier(cctx->appliedParams.blockDelimiters);

    DEBUGLOG(4, "ZSTD_compressSequences_internal srcSize: %zu, inSeqsSize: %zu", srcSize, inSeqsSize);
    /* Special case: empty frame */
    if (remaining == 0) {
        U32 const cBlockHeader24 = 1 /* last block */ + (((U32)bt_raw)<<1);
        RETURN_ERROR_IF(dstCapacity<4, dstSize_tooSmall, "No room for empty frame block header");
        MEM_writeLE32(op, cBlockHeader24);
        op += ZSTD_blockHeaderSize;
        dstCapacity -= ZSTD_blockHeaderSize;
        cSize += ZSTD_blockHeaderSize;
    }

    while (remaining) {
        size_t compressedSeqsSize;
        size_t cBlockSize;
        size_t additionalByteAdjustment;
        size_t blockSize = determine_blockSize(cctx->appliedParams.blockDelimiters,
                                        cctx->blockSize, remaining,
                                        inSeqs, inSeqsSize, seqPos);
        U32 const lastBlock = (blockSize == remaining);
        FORWARD_IF_ERROR(blockSize, "Error while trying to determine block size");
        assert(blockSize <= remaining);
        ZSTD_resetSeqStore(&cctx->seqStore);
        DEBUGLOG(5, "Working on new block. Blocksize: %zu (total:%zu)", blockSize, (ip - (const BYTE*)src) + blockSize);

        additionalByteAdjustment = sequenceCopier(cctx, &seqPos, inSeqs, inSeqsSize, ip, blockSize, cctx->appliedParams.searchForExternalRepcodes);
        FORWARD_IF_ERROR(additionalByteAdjustment, "Bad sequence copy");
        blockSize -= additionalByteAdjustment;

        /* If blocks are too small, emit as a nocompress block */
        /* TODO: See 3090. We reduced MIN_CBLOCK_SIZE from 3 to 2 so to compensate we are adding
         * additional 1. We need to revisit and change this logic to be more consistent */
        if (blockSize < MIN_CBLOCK_SIZE+ZSTD_blockHeaderSize+1+1) {
            cBlockSize = ZSTD_noCompressBlock(op, dstCapacity, ip, blockSize, lastBlock);
            FORWARD_IF_ERROR(cBlockSize, "Nocompress block failed");
            DEBUGLOG(5, "Block too small, writing out nocompress block: cSize: %zu", cBlockSize);
            cSize += cBlockSize;
            ip += blockSize;
            op += cBlockSize;
            remaining -= blockSize;
            dstCapacity -= cBlockSize;
            continue;
        }

        RETURN_ERROR_IF(dstCapacity < ZSTD_blockHeaderSize, dstSize_tooSmall, "not enough dstCapacity to write a new compressed block");
        compressedSeqsSize = ZSTD_entropyCompressSeqStore(&cctx->seqStore,
                                &cctx->blockState.prevCBlock->entropy, &cctx->blockState.nextCBlock->entropy,
                                &cctx->appliedParams,
                                op + ZSTD_blockHeaderSize /* Leave space for block header */, dstCapacity - ZSTD_blockHeaderSize,
                                blockSize,
                                cctx->entropyWorkspace, ENTROPY_WORKSPACE_SIZE /* statically allocated in resetCCtx */,
                                cctx->bmi2);
        FORWARD_IF_ERROR(compressedSeqsSize, "Compressing sequences of block failed");
        DEBUGLOG(5, "Compressed sequences size: %zu", compressedSeqsSize);

        if (!cctx->isFirstBlock &&
            ZSTD_maybeRLE(&cctx->seqStore) &&
            ZSTD_isRLE(ip, blockSize)) {
            /* We don't want to emit our first block as a RLE even if it qualifies because
            * doing so will cause the decoder (cli only) to throw a "should consume all input error."
            * This is only an issue for zstd <= v1.4.3
            */
            compressedSeqsSize = 1;
        }

        if (compressedSeqsSize == 0) {
            /* ZSTD_noCompressBlock writes the block header as well */
            cBlockSize = ZSTD_noCompressBlock(op, dstCapacity, ip, blockSize, lastBlock);
            FORWARD_IF_ERROR(cBlockSize, "ZSTD_noCompressBlock failed");
            DEBUGLOG(5, "Writing out nocompress block, size: %zu", cBlockSize);
        } else if (compressedSeqsSize == 1) {
            cBlockSize = ZSTD_rleCompressBlock(op, dstCapacity, *ip, blockSize, lastBlock);
            FORWARD_IF_ERROR(cBlockSize, "ZSTD_rleCompressBlock failed");
            DEBUGLOG(5, "Writing out RLE block, size: %zu", cBlockSize);
        } else {
            U32 cBlockHeader;
            /* Error checking and repcodes update */
            ZSTD_blockState_confirmRepcodesAndEntropyTables(&cctx->blockState);
            if (cctx->blockState.prevCBlock->entropy.fse.offcode_repeatMode == FSE_repeat_valid)
                cctx->blockState.prevCBlock->entropy.fse.offcode_repeatMode = FSE_repeat_check;

            /* Write block header into beginning of block*/
            cBlockHeader = lastBlock + (((U32)bt_compressed)<<1) + (U32)(compressedSeqsSize << 3);
            MEM_writeLE24(op, cBlockHeader);
            cBlockSize = ZSTD_blockHeaderSize + compressedSeqsSize;
            DEBUGLOG(5, "Writing out compressed block, size: %zu", cBlockSize);
        }

        cSize += cBlockSize;

        if (lastBlock) {
            break;
        } else {
            ip += blockSize;
            op += cBlockSize;
            remaining -= blockSize;
            dstCapacity -= cBlockSize;
            cctx->isFirstBlock = 0;
        }
        DEBUGLOG(5, "cSize running total: %zu (remaining dstCapacity=%zu)", cSize, dstCapacity);
    }

    DEBUGLOG(4, "cSize final total: %zu", cSize);
    return cSize;
}

size_t ZSTD_compressSequences(ZSTD_CCtx* cctx,
                              void* dst, size_t dstCapacity,
                              const ZSTD_Sequence* inSeqs, size_t inSeqsSize,
                              const void* src, size_t srcSize)
{
    BYTE* op = (BYTE*)dst;
    size_t cSize = 0;
    size_t compressedBlocksSize = 0;
    size_t frameHeaderSize = 0;

    /* Transparent initialization stage, same as compressStream2() */
    DEBUGLOG(4, "ZSTD_compressSequences (dstCapacity=%zu)", dstCapacity);
    assert(cctx != NULL);
    FORWARD_IF_ERROR(ZSTD_CCtx_init_compressStream2(cctx, ZSTD_e_end, srcSize), "CCtx initialization failed");
    /* Begin writing output, starting with frame header */
    frameHeaderSize = ZSTD_writeFrameHeader(op, dstCapacity, &cctx->appliedParams, srcSize, cctx->dictID);
    op += frameHeaderSize;
    dstCapacity -= frameHeaderSize;
    cSize += frameHeaderSize;
    if (cctx->appliedParams.fParams.checksumFlag && srcSize) {
        XXH64_update(&cctx->xxhState, src, srcSize);
    }
    /* cSize includes block header size and compressed sequences size */
    compressedBlocksSize = ZSTD_compressSequences_internal(cctx,
                                                           op, dstCapacity,
                                                           inSeqs, inSeqsSize,
                                                           src, srcSize);
    FORWARD_IF_ERROR(compressedBlocksSize, "Compressing blocks failed!");
    cSize += compressedBlocksSize;
    dstCapacity -= compressedBlocksSize;

    if (cctx->appliedParams.fParams.checksumFlag) {
        U32 const checksum = (U32) XXH64_digest(&cctx->xxhState);
        RETURN_ERROR_IF(dstCapacity<4, dstSize_tooSmall, "no room for checksum");
        DEBUGLOG(4, "Write checksum : %08X", (unsigned)checksum);
        MEM_writeLE32((char*)dst + cSize, checksum);
        cSize += 4;
    }

    DEBUGLOG(4, "Final compressed size: %zu", cSize);
    return cSize;
}

/*======   Finalize   ======*/

static ZSTD_inBuffer inBuffer_forEndFlush(const ZSTD_CStream* zcs)
{
    const ZSTD_inBuffer nullInput = { NULL, 0, 0 };
    const int stableInput = (zcs->appliedParams.inBufferMode == ZSTD_bm_stable);
    return stableInput ? zcs->expectedInBuffer : nullInput;
}

/*! ZSTD_flushStream() :
 * @return : amount of data remaining to flush */
size_t ZSTD_flushStream(ZSTD_CStream* zcs, ZSTD_outBuffer* output)
{
    ZSTD_inBuffer input = inBuffer_forEndFlush(zcs);
    input.size = input.pos; /* do not ingest more input during flush */
    return ZSTD_compressStream2(zcs, output, &input, ZSTD_e_flush);
}


size_t ZSTD_endStream(ZSTD_CStream* zcs, ZSTD_outBuffer* output)
{
    ZSTD_inBuffer input = inBuffer_forEndFlush(zcs);
    size_t const remainingToFlush = ZSTD_compressStream2(zcs, output, &input, ZSTD_e_end);
    FORWARD_IF_ERROR(remainingToFlush , "ZSTD_compressStream2(,,ZSTD_e_end) failed");
    if (zcs->appliedParams.nbWorkers > 0) return remainingToFlush;   /* minimal estimation */
    /* single thread mode : attempt to calculate remaining to flush more precisely */
    {   size_t const lastBlockSize = zcs->frameEnded ? 0 : ZSTD_BLOCKHEADERSIZE;
        size_t const checksumSize = (size_t)(zcs->frameEnded ? 0 : zcs->appliedParams.fParams.checksumFlag * 4);
        size_t const toFlush = remainingToFlush + lastBlockSize + checksumSize;
        DEBUGLOG(4, "ZSTD_endStream : remaining to flush : %u", (unsigned)toFlush);
        return toFlush;
    }
}


/*-=====  Pre-defined compression levels  =====-*/
#include "clevels.h"

int ZSTD_maxCLevel(void) { return ZSTD_MAX_CLEVEL; }
int ZSTD_minCLevel(void) { return (int)-ZSTD_TARGETLENGTH_MAX; }
int ZSTD_defaultCLevel(void) { return ZSTD_CLEVEL_DEFAULT; }

static ZSTD_compressionParameters ZSTD_dedicatedDictSearch_getCParams(int const compressionLevel, size_t const dictSize)
{
    ZSTD_compressionParameters cParams = ZSTD_getCParams_internal(compressionLevel, 0, dictSize, ZSTD_cpm_createCDict);
    switch (cParams.strategy) {
        case ZSTD_fast:
        case ZSTD_dfast:
            break;
        case ZSTD_greedy:
        case ZSTD_lazy:
        case ZSTD_lazy2:
            cParams.hashLog += ZSTD_LAZY_DDSS_BUCKET_LOG;
            break;
        case ZSTD_btlazy2:
        case ZSTD_btopt:
        case ZSTD_btultra:
        case ZSTD_btultra2:
            break;
    }
    return cParams;
}

static int ZSTD_dedicatedDictSearch_isSupported(
        ZSTD_compressionParameters const* cParams)
{
    return (cParams->strategy >= ZSTD_greedy)
        && (cParams->strategy <= ZSTD_lazy2)
        && (cParams->hashLog > cParams->chainLog)
        && (cParams->chainLog <= 24);
}

/**
 * Reverses the adjustment applied to cparams when enabling dedicated dict
 * search. This is used to recover the params set to be used in the working
 * context. (Otherwise, those tables would also grow.)
 */
static void ZSTD_dedicatedDictSearch_revertCParams(
        ZSTD_compressionParameters* cParams) {
    switch (cParams->strategy) {
        case ZSTD_fast:
        case ZSTD_dfast:
            break;
        case ZSTD_greedy:
        case ZSTD_lazy:
        case ZSTD_lazy2:
            cParams->hashLog -= ZSTD_LAZY_DDSS_BUCKET_LOG;
            if (cParams->hashLog < ZSTD_HASHLOG_MIN) {
                cParams->hashLog = ZSTD_HASHLOG_MIN;
            }
            break;
        case ZSTD_btlazy2:
        case ZSTD_btopt:
        case ZSTD_btultra:
        case ZSTD_btultra2:
            break;
    }
}

static U64 ZSTD_getCParamRowSize(U64 srcSizeHint, size_t dictSize, ZSTD_cParamMode_e mode)
{
    switch (mode) {
    case ZSTD_cpm_unknown:
    case ZSTD_cpm_noAttachDict:
    case ZSTD_cpm_createCDict:
        break;
    case ZSTD_cpm_attachDict:
        dictSize = 0;
        break;
    default:
        assert(0);
        break;
    }
    {   int const unknown = srcSizeHint == ZSTD_CONTENTSIZE_UNKNOWN;
        size_t const addedSize = unknown && dictSize > 0 ? 500 : 0;
        return unknown && dictSize == 0 ? ZSTD_CONTENTSIZE_UNKNOWN : srcSizeHint+dictSize+addedSize;
    }
}

/*! ZSTD_getCParams_internal() :
 * @return ZSTD_compressionParameters structure for a selected compression level, srcSize and dictSize.
 *  Note: srcSizeHint 0 means 0, use ZSTD_CONTENTSIZE_UNKNOWN for unknown.
 *        Use dictSize == 0 for unknown or unused.
 *  Note: `mode` controls how we treat the `dictSize`. See docs for `ZSTD_cParamMode_e`. */
static ZSTD_compressionParameters ZSTD_getCParams_internal(int compressionLevel, unsigned long long srcSizeHint, size_t dictSize, ZSTD_cParamMode_e mode)
{
    U64 const rSize = ZSTD_getCParamRowSize(srcSizeHint, dictSize, mode);
    U32 const tableID = (rSize <= 256 KB) + (rSize <= 128 KB) + (rSize <= 16 KB);
    int row;
    DEBUGLOG(5, "ZSTD_getCParams_internal (cLevel=%i)", compressionLevel);

    /* row */
    if (compressionLevel == 0) row = ZSTD_CLEVEL_DEFAULT;   /* 0 == default */
    else if (compressionLevel < 0) row = 0;   /* entry 0 is baseline for fast mode */
    else if (compressionLevel > ZSTD_MAX_CLEVEL) row = ZSTD_MAX_CLEVEL;
    else row = compressionLevel;

    {   ZSTD_compressionParameters cp = ZSTD_defaultCParameters[tableID][row];
        DEBUGLOG(5, "ZSTD_getCParams_internal selected tableID: %u row: %u strat: %u", tableID, row, (U32)cp.strategy);
        /* acceleration factor */
        if (compressionLevel < 0) {
            int const clampedCompressionLevel = MAX(ZSTD_minCLevel(), compressionLevel);
            cp.targetLength = (unsigned)(-clampedCompressionLevel);
        }
        /* refine parameters based on srcSize & dictSize */
        return ZSTD_adjustCParams_internal(cp, srcSizeHint, dictSize, mode, ZSTD_ps_auto);
    }
}

/*! ZSTD_getCParams() :
 * @return ZSTD_compressionParameters structure for a selected compression level, srcSize and dictSize.
 *  Size values are optional, provide 0 if not known or unused */
ZSTD_compressionParameters ZSTD_getCParams(int compressionLevel, unsigned long long srcSizeHint, size_t dictSize)
{
    if (srcSizeHint == 0) srcSizeHint = ZSTD_CONTENTSIZE_UNKNOWN;
    return ZSTD_getCParams_internal(compressionLevel, srcSizeHint, dictSize, ZSTD_cpm_unknown);
}

/*! ZSTD_getParams() :
 *  same idea as ZSTD_getCParams()
 * @return a `ZSTD_parameters` structure (instead of `ZSTD_compressionParameters`).
 *  Fields of `ZSTD_frameParameters` are set to default values */
static ZSTD_parameters ZSTD_getParams_internal(int compressionLevel, unsigned long long srcSizeHint, size_t dictSize, ZSTD_cParamMode_e mode) {
    ZSTD_parameters params;
    ZSTD_compressionParameters const cParams = ZSTD_getCParams_internal(compressionLevel, srcSizeHint, dictSize, mode);
    DEBUGLOG(5, "ZSTD_getParams (cLevel=%i)", compressionLevel);
    ZSTD_memset(&params, 0, sizeof(params));
    params.cParams = cParams;
    params.fParams.contentSizeFlag = 1;
    return params;
}

/*! ZSTD_getParams() :
 *  same idea as ZSTD_getCParams()
 * @return a `ZSTD_parameters` structure (instead of `ZSTD_compressionParameters`).
 *  Fields of `ZSTD_frameParameters` are set to default values */
ZSTD_parameters ZSTD_getParams(int compressionLevel, unsigned long long srcSizeHint, size_t dictSize) {
    if (srcSizeHint == 0) srcSizeHint = ZSTD_CONTENTSIZE_UNKNOWN;
    return ZSTD_getParams_internal(compressionLevel, srcSizeHint, dictSize, ZSTD_cpm_unknown);
}

void ZSTD_registerSequenceProducer(
    ZSTD_CCtx* zc,
    void* extSeqProdState,
    ZSTD_sequenceProducer_F extSeqProdFunc
) {
    assert(zc != NULL);
    ZSTD_CCtxParams_registerSequenceProducer(
        &zc->requestedParams, extSeqProdState, extSeqProdFunc
    );
}

void ZSTD_CCtxParams_registerSequenceProducer(
  ZSTD_CCtx_params* params,
  void* extSeqProdState,
  ZSTD_sequenceProducer_F extSeqProdFunc
) {
    assert(params != NULL);
    if (extSeqProdFunc != NULL) {
        params->extSeqProdFunc = extSeqProdFunc;
        params->extSeqProdState = extSeqProdState;
    } else {
        params->extSeqProdFunc = NULL;
        params->extSeqProdState = NULL;
    }
}
