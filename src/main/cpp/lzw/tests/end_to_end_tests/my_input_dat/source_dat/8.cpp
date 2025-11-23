/* ******************************************************************
 * Huffman encoder, part of New Generation Entropy library
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 *  You can contact the author at :
 *  - FSE+HUF source repository : https://github.com/Cyan4973/FiniteStateEntropy
 *  - Public forum : https://groups.google.com/forum/#!forum/lz4c
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
****************************************************************** */

/* **************************************************************
*  Compiler specifics
****************************************************************/
#ifdef _MSC_VER    /* Visual Studio */
#  pragma warning(disable : 4127)        /* disable: C4127: conditional expression is constant */
#endif


/* **************************************************************
*  Includes
****************************************************************/
#include "../common/zstd_deps.h"     /* ZSTD_memcpy, ZSTD_memset */
#include "../common/compiler.h"
#include "../common/bitstream.h"
#include "hist.h"
#define FSE_STATIC_LINKING_ONLY   /* FSE_optimalTableLog_internal */
#include "../common/fse.h"        /* header compression */
#include "../common/huf.h"
#include "../common/error_private.h"
#include "../common/bits.h"       /* ZSTD_highbit32 */


/* **************************************************************
*  Error Management
****************************************************************/
#define HUF_isError ERR_isError
#define HUF_STATIC_ASSERT(c) DEBUG_STATIC_ASSERT(c)   /* use only *after* variable declarations */


/* **************************************************************
*  Required declarations
****************************************************************/
typedef struct nodeElt_s {
    U32 count;
    U16 parent;
    BYTE byte;
    BYTE nbBits;
} nodeElt;


/* **************************************************************
*  Debug Traces
****************************************************************/

#if DEBUGLEVEL >= 2

static size_t showU32(const U32* arr, size_t size)
{
    size_t u;
    for (u=0; u<size; u++) {
        RAWLOG(6, " %u", arr[u]); (void)arr;
    }
    RAWLOG(6, " \n");
    return size;
}

static size_t HUF_getNbBits(HUF_CElt elt);

static size_t showCTableBits(const HUF_CElt* ctable, size_t size)
{
    size_t u;
    for (u=0; u<size; u++) {
        RAWLOG(6, " %zu", HUF_getNbBits(ctable[u])); (void)ctable;
    }
    RAWLOG(6, " \n");
    return size;

}

static size_t showHNodeSymbols(const nodeElt* hnode, size_t size)
{
    size_t u;
    for (u=0; u<size; u++) {
        RAWLOG(6, " %u", hnode[u].byte); (void)hnode;
    }
    RAWLOG(6, " \n");
    return size;
}

static size_t showHNodeBits(const nodeElt* hnode, size_t size)
{
    size_t u;
    for (u=0; u<size; u++) {
        RAWLOG(6, " %u", hnode[u].nbBits); (void)hnode;
    }
    RAWLOG(6, " \n");
    return size;
}

#endif


/* *******************************************************
*  HUF : Huffman block compression
*********************************************************/
#define HUF_WORKSPACE_MAX_ALIGNMENT 8

static void* HUF_alignUpWorkspace(void* workspace, size_t* workspaceSizePtr, size_t align)
{
    size_t const mask = align - 1;
    size_t const rem = (size_t)workspace & mask;
    size_t const add = (align - rem) & mask;
    BYTE* const aligned = (BYTE*)workspace + add;
    assert((align & (align - 1)) == 0); /* pow 2 */
    assert(align <= HUF_WORKSPACE_MAX_ALIGNMENT);
    if (*workspaceSizePtr >= add) {
        assert(add < align);
        assert(((size_t)aligned & mask) == 0);
        *workspaceSizePtr -= add;
        return aligned;
    } else {
        *workspaceSizePtr = 0;
        return NULL;
    }
}


/* HUF_compressWeights() :
 * Same as FSE_compress(), but dedicated to huff0's weights compression.
 * The use case needs much less stack memory.
 * Note : all elements within weightTable are supposed to be <= HUF_TABLELOG_MAX.
 */
#define MAX_FSE_TABLELOG_FOR_HUFF_HEADER 6

typedef struct {
    FSE_CTable CTable[FSE_CTABLE_SIZE_U32(MAX_FSE_TABLELOG_FOR_HUFF_HEADER, HUF_TABLELOG_MAX)];
    U32 scratchBuffer[FSE_BUILD_CTABLE_WORKSPACE_SIZE_U32(HUF_TABLELOG_MAX, MAX_FSE_TABLELOG_FOR_HUFF_HEADER)];
    unsigned count[HUF_TABLELOG_MAX+1];
    S16 norm[HUF_TABLELOG_MAX+1];
} HUF_CompressWeightsWksp;

static size_t
HUF_compressWeights(void* dst, size_t dstSize,
              const void* weightTable, size_t wtSize,
                    void* workspace, size_t workspaceSize)
{
    BYTE* const ostart = (BYTE*) dst;
    BYTE* op = ostart;
    BYTE* const oend = ostart + dstSize;

    unsigned maxSymbolValue = HUF_TABLELOG_MAX;
    U32 tableLog = MAX_FSE_TABLELOG_FOR_HUFF_HEADER;
    HUF_CompressWeightsWksp* wksp = (HUF_CompressWeightsWksp*)HUF_alignUpWorkspace(workspace, &workspaceSize, ZSTD_ALIGNOF(U32));

    if (workspaceSize < sizeof(HUF_CompressWeightsWksp)) return ERROR(GENERIC);

    /* init conditions */
    if (wtSize <= 1) return 0;  /* Not compressible */

    /* Scan input and build symbol stats */
    {   unsigned const maxCount = HIST_count_simple(wksp->count, &maxSymbolValue, weightTable, wtSize);   /* never fails */
        if (maxCount == wtSize) return 1;   /* only a single symbol in src : rle */
        if (maxCount == 1) return 0;        /* each symbol present maximum once => not compressible */
    }

    tableLog = FSE_optimalTableLog(tableLog, wtSize, maxSymbolValue);
    CHECK_F( FSE_normalizeCount(wksp->norm, tableLog, wksp->count, wtSize, maxSymbolValue, /* useLowProbCount */ 0) );

    /* Write table description header */
    {   CHECK_V_F(hSize, FSE_writeNCount(op, (size_t)(oend-op), wksp->norm, maxSymbolValue, tableLog) );
        op += hSize;
    }

    /* Compress */
    CHECK_F( FSE_buildCTable_wksp(wksp->CTable, wksp->norm, maxSymbolValue, tableLog, wksp->scratchBuffer, sizeof(wksp->scratchBuffer)) );
    {   CHECK_V_F(cSize, FSE_compress_usingCTable(op, (size_t)(oend - op), weightTable, wtSize, wksp->CTable) );
        if (cSize == 0) return 0;   /* not enough space for compressed data */
        op += cSize;
    }

    return (size_t)(op-ostart);
}

static size_t HUF_getNbBits(HUF_CElt elt)
{
    return elt & 0xFF;
}

static size_t HUF_getNbBitsFast(HUF_CElt elt)
{
    return elt;
}

static size_t HUF_getValue(HUF_CElt elt)
{
    return elt & ~(size_t)0xFF;
}

static size_t HUF_getValueFast(HUF_CElt elt)
{
    return elt;
}

static void HUF_setNbBits(HUF_CElt* elt, size_t nbBits)
{
    assert(nbBits <= HUF_TABLELOG_ABSOLUTEMAX);
    *elt = nbBits;
}

static void HUF_setValue(HUF_CElt* elt, size_t value)
{
    size_t const nbBits = HUF_getNbBits(*elt);
    if (nbBits > 0) {
        assert((value >> nbBits) == 0);
        *elt |= value << (sizeof(HUF_CElt) * 8 - nbBits);
    }
}

HUF_CTableHeader HUF_readCTableHeader(HUF_CElt const* ctable)
{
    HUF_CTableHeader header;
    ZSTD_memcpy(&header, ctable, sizeof(header));
    return header;
}

static void HUF_writeCTableHeader(HUF_CElt* ctable, U32 tableLog, U32 maxSymbolValue)
{
    HUF_CTableHeader header;
    HUF_STATIC_ASSERT(sizeof(ctable[0]) == sizeof(header));
    ZSTD_memset(&header, 0, sizeof(header));
    assert(tableLog < 256);
    header.tableLog = (BYTE)tableLog;
    assert(maxSymbolValue < 256);
    header.maxSymbolValue = (BYTE)maxSymbolValue;
    ZSTD_memcpy(ctable, &header, sizeof(header));
}

typedef struct {
    HUF_CompressWeightsWksp wksp;
    BYTE bitsToWeight[HUF_TABLELOG_MAX + 1];   /* precomputed conversion table */
    BYTE huffWeight[HUF_SYMBOLVALUE_MAX];
} HUF_WriteCTableWksp;

size_t HUF_writeCTable_wksp(void* dst, size_t maxDstSize,
                            const HUF_CElt* CTable, unsigned maxSymbolValue, unsigned huffLog,
                            void* workspace, size_t workspaceSize)
{
    HUF_CElt const* const ct = CTable + 1;
    BYTE* op = (BYTE*)dst;
    U32 n;
    HUF_WriteCTableWksp* wksp = (HUF_WriteCTableWksp*)HUF_alignUpWorkspace(workspace, &workspaceSize, ZSTD_ALIGNOF(U32));

    HUF_STATIC_ASSERT(HUF_CTABLE_WORKSPACE_SIZE >= sizeof(HUF_WriteCTableWksp));

    assert(HUF_readCTableHeader(CTable).maxSymbolValue == maxSymbolValue);
    assert(HUF_readCTableHeader(CTable).tableLog == huffLog);

    /* check conditions */
    if (workspaceSize < sizeof(HUF_WriteCTableWksp)) return ERROR(GENERIC);
    if (maxSymbolValue > HUF_SYMBOLVALUE_MAX) return ERROR(maxSymbolValue_tooLarge);

    /* convert to weight */
    wksp->bitsToWeight[0] = 0;
    for (n=1; n<huffLog+1; n++)
        wksp->bitsToWeight[n] = (BYTE)(huffLog + 1 - n);
    for (n=0; n<maxSymbolValue; n++)
        wksp->huffWeight[n] = wksp->bitsToWeight[HUF_getNbBits(ct[n])];

    /* attempt weights compression by FSE */
    if (maxDstSize < 1) return ERROR(dstSize_tooSmall);
    {   CHECK_V_F(hSize, HUF_compressWeights(op+1, maxDstSize-1, wksp->huffWeight, maxSymbolValue, &wksp->wksp, sizeof(wksp->wksp)) );
        if ((hSize>1) & (hSize < maxSymbolValue/2)) {   /* FSE compressed */
            op[0] = (BYTE)hSize;
            return hSize+1;
    }   }

    /* write raw values as 4-bits (max : 15) */
    if (maxSymbolValue > (256-128)) return ERROR(GENERIC);   /* should not happen : likely means source cannot be compressed */
    if (((maxSymbolValue+1)/2) + 1 > maxDstSize) return ERROR(dstSize_tooSmall);   /* not enough space within dst buffer */
    op[0] = (BYTE)(128 /*special case*/ + (maxSymbolValue-1));
    wksp->huffWeight[maxSymbolValue] = 0;   /* to be sure it doesn't cause msan issue in final combination */
    for (n=0; n<maxSymbolValue; n+=2)
        op[(n/2)+1] = (BYTE)((wksp->huffWeight[n] << 4) + wksp->huffWeight[n+1]);
    return ((maxSymbolValue+1)/2) + 1;
}


size_t HUF_readCTable (HUF_CElt* CTable, unsigned* maxSymbolValuePtr, const void* src, size_t srcSize, unsigned* hasZeroWeights)
{
    BYTE huffWeight[HUF_SYMBOLVALUE_MAX + 1];   /* init not required, even though some static analyzer may complain */
    U32 rankVal[HUF_TABLELOG_ABSOLUTEMAX + 1];   /* large enough for values from 0 to 16 */
    U32 tableLog = 0;
    U32 nbSymbols = 0;
    HUF_CElt* const ct = CTable + 1;

    /* get symbol weights */
    CHECK_V_F(readSize, HUF_readStats(huffWeight, HUF_SYMBOLVALUE_MAX+1, rankVal, &nbSymbols, &tableLog, src, srcSize));
    *hasZeroWeights = (rankVal[0] > 0);

    /* check result */
    if (tableLog > HUF_TABLELOG_MAX) return ERROR(tableLog_tooLarge);
    if (nbSymbols > *maxSymbolValuePtr+1) return ERROR(maxSymbolValue_tooSmall);

    *maxSymbolValuePtr = nbSymbols - 1;

    HUF_writeCTableHeader(CTable, tableLog, *maxSymbolValuePtr);

    /* Prepare base value per rank */
    {   U32 n, nextRankStart = 0;
        for (n=1; n<=tableLog; n++) {
            U32 curr = nextRankStart;
            nextRankStart += (rankVal[n] << (n-1));
            rankVal[n] = curr;
    }   }

    /* fill nbBits */
    {   U32 n; for (n=0; n<nbSymbols; n++) {
            const U32 w = huffWeight[n];
            HUF_setNbBits(ct + n, (BYTE)(tableLog + 1 - w) & -(w != 0));
    }   }

    /* fill val */
    {   U16 nbPerRank[HUF_TABLELOG_MAX+2]  = {0};  /* support w=0=>n=tableLog+1 */
        U16 valPerRank[HUF_TABLELOG_MAX+2] = {0};
        { U32 n; for (n=0; n<nbSymbols; n++) nbPerRank[HUF_getNbBits(ct[n])]++; }
        /* determine stating value per rank */
        valPerRank[tableLog+1] = 0;   /* for w==0 */
        {   U16 min = 0;
            U32 n; for (n=tableLog; n>0; n--) {  /* start at n=tablelog <-> w=1 */
                valPerRank[n] = min;     /* get starting value within each rank */
                min += nbPerRank[n];
                min >>= 1;
        }   }
        /* assign value within rank, symbol order */
        { U32 n; for (n=0; n<nbSymbols; n++) HUF_setValue(ct + n, valPerRank[HUF_getNbBits(ct[n])]++); }
    }

    return readSize;
}

U32 HUF_getNbBitsFromCTable(HUF_CElt const* CTable, U32 symbolValue)
{
    const HUF_CElt* const ct = CTable + 1;
    assert(symbolValue <= HUF_SYMBOLVALUE_MAX);
    if (symbolValue > HUF_readCTableHeader(CTable).maxSymbolValue)
        return 0;
    return (U32)HUF_getNbBits(ct[symbolValue]);
}


/**
 * HUF_setMaxHeight():
 * Try to enforce @targetNbBits on the Huffman tree described in @huffNode.
 *
 * It attempts to convert all nodes with nbBits > @targetNbBits
 * to employ @targetNbBits instead. Then it adjusts the tree
 * so that it remains a valid canonical Huffman tree.
 *
 * @pre               The sum of the ranks of each symbol == 2^largestBits,
 *                    where largestBits == huffNode[lastNonNull].nbBits.
 * @post              The sum of the ranks of each symbol == 2^largestBits,
 *                    where largestBits is the return value (expected <= targetNbBits).
 *
 * @param huffNode    The Huffman tree modified in place to enforce targetNbBits.
 *                    It's presumed sorted, from most frequent to rarest symbol.
 * @param lastNonNull The symbol with the lowest count in the Huffman tree.
 * @param targetNbBits  The allowed number of bits, which the Huffman tree
 *                    may not respect. After this function the Huffman tree will
 *                    respect targetNbBits.
 * @return            The maximum number of bits of the Huffman tree after adjustment.
 */
static U32 HUF_setMaxHeight(nodeElt* huffNode, U32 lastNonNull, U32 targetNbBits)
{
    const U32 largestBits = huffNode[lastNonNull].nbBits;
    /* early exit : no elt > targetNbBits, so the tree is already valid. */
    if (largestBits <= targetNbBits) return largestBits;

    DEBUGLOG(5, "HUF_setMaxHeight (targetNbBits = %u)", targetNbBits);

    /* there are several too large elements (at least >= 2) */
    {   int totalCost = 0;
        const U32 baseCost = 1 << (largestBits - targetNbBits);
        int n = (int)lastNonNull;

        /* Adjust any ranks > targetNbBits to targetNbBits.
         * Compute totalCost, which is how far the sum of the ranks is
         * we are over 2^largestBits after adjust the offending ranks.
         */
        while (huffNode[n].nbBits > targetNbBits) {
            totalCost += baseCost - (1 << (largestBits - huffNode[n].nbBits));
            huffNode[n].nbBits = (BYTE)targetNbBits;
            n--;
        }
        /* n stops at huffNode[n].nbBits <= targetNbBits */
        assert(huffNode[n].nbBits <= targetNbBits);
        /* n end at index of smallest symbol using < targetNbBits */
        while (huffNode[n].nbBits == targetNbBits) --n;

        /* renorm totalCost from 2^largestBits to 2^targetNbBits
         * note : totalCost is necessarily a multiple of baseCost */
        assert(((U32)totalCost & (baseCost - 1)) == 0);
        totalCost >>= (largestBits - targetNbBits);
        assert(totalCost > 0);

        /* repay normalized cost */
        {   U32 const noSymbol = 0xF0F0F0F0;
            U32 rankLast[HUF_TABLELOG_MAX+2];

            /* Get pos of last (smallest = lowest cum. count) symbol per rank */
            ZSTD_memset(rankLast, 0xF0, sizeof(rankLast));
            {   U32 currentNbBits = targetNbBits;
                int pos;
                for (pos=n ; pos >= 0; pos--) {
                    if (huffNode[pos].nbBits >= currentNbBits) continue;
                    currentNbBits = huffNode[pos].nbBits;   /* < targetNbBits */
                    rankLast[targetNbBits-currentNbBits] = (U32)pos;
            }   }

            while (totalCost > 0) {
                /* Try to reduce the next power of 2 above totalCost because we
                 * gain back half the rank.
                 */
                U32 nBitsToDecrease = ZSTD_highbit32((U32)totalCost) + 1;
                for ( ; nBitsToDecrease > 1; nBitsToDecrease--) {
                    U32 const highPos = rankLast[nBitsToDecrease];
                    U32 const lowPos = rankLast[nBitsToDecrease-1];
                    if (highPos == noSymbol) continue;
                    /* Decrease highPos if no symbols of lowPos or if it is
                     * not cheaper to remove 2 lowPos than highPos.
                     */
                    if (lowPos == noSymbol) break;
                    {   U32 const highTotal = huffNode[highPos].count;
                        U32 const lowTotal = 2 * huffNode[lowPos].count;
                        if (highTotal <= lowTotal) break;
                }   }
                /* only triggered when no more rank 1 symbol left => find closest one (note : there is necessarily at least one !) */
                assert(rankLast[nBitsToDecrease] != noSymbol || nBitsToDecrease == 1);
                /* HUF_MAX_TABLELOG test just to please gcc 5+; but it should not be necessary */
                while ((nBitsToDecrease<=HUF_TABLELOG_MAX) && (rankLast[nBitsToDecrease] == noSymbol))
                    nBitsToDecrease++;
                assert(rankLast[nBitsToDecrease] != noSymbol);
                /* Increase the number of bits to gain back half the rank cost. */
                totalCost -= 1 << (nBitsToDecrease-1);
                huffNode[rankLast[nBitsToDecrease]].nbBits++;

                /* Fix up the new rank.
                 * If the new rank was empty, this symbol is now its smallest.
                 * Otherwise, this symbol will be the largest in the new rank so no adjustment.
                 */
                if (rankLast[nBitsToDecrease-1] == noSymbol)
                    rankLast[nBitsToDecrease-1] = rankLast[nBitsToDecrease];
                /* Fix up the old rank.
                 * If the symbol was at position 0, meaning it was the highest weight symbol in the tree,
                 * it must be the only symbol in its rank, so the old rank now has no symbols.
                 * Otherwise, since the Huffman nodes are sorted by count, the previous position is now
                 * the smallest node in the rank. If the previous position belongs to a different rank,
                 * then the rank is now empty.
                 */
                if (rankLast[nBitsToDecrease] == 0)    /* special case, reached largest symbol */
                    rankLast[nBitsToDecrease] = noSymbol;
                else {
                    rankLast[nBitsToDecrease]--;
                    if (huffNode[rankLast[nBitsToDecrease]].nbBits != targetNbBits-nBitsToDecrease)
                        rankLast[nBitsToDecrease] = noSymbol;   /* this rank is now empty */
                }
            }   /* while (totalCost > 0) */

            /* If we've removed too much weight, then we have to add it back.
             * To avoid overshooting again, we only adjust the smallest rank.
             * We take the largest nodes from the lowest rank 0 and move them
             * to rank 1. There's guaranteed to be enough rank 0 symbols because
             * TODO.
             */
            while (totalCost < 0) {  /* Sometimes, cost correction overshoot */
                /* special case : no rank 1 symbol (using targetNbBits-1);
                 * let's create one from largest rank 0 (using targetNbBits).
                 */
                if (rankLast[1] == noSymbol) {
                    while (huffNode[n].nbBits == targetNbBits) n--;
                    huffNode[n+1].nbBits--;
                    assert(n >= 0);
                    rankLast[1] = (U32)(n+1);
                    totalCost++;
                    continue;
                }
                huffNode[ rankLast[1] + 1 ].nbBits--;
                rankLast[1]++;
                totalCost ++;
            }
        }   /* repay normalized cost */
    }   /* there are several too large elements (at least >= 2) */

    return targetNbBits;
}

typedef struct {
    U16 base;
    U16 curr;
} rankPos;

typedef nodeElt huffNodeTable[2 * (HUF_SYMBOLVALUE_MAX + 1)];

/* Number of buckets available for HUF_sort() */
#define RANK_POSITION_TABLE_SIZE 192

typedef struct {
  huffNodeTable huffNodeTbl;
  rankPos rankPosition[RANK_POSITION_TABLE_SIZE];
} HUF_buildCTable_wksp_tables;

/* RANK_POSITION_DISTINCT_COUNT_CUTOFF == Cutoff point in HUF_sort() buckets for which we use log2 bucketing.
 * Strategy is to use as many buckets as possible for representing distinct
 * counts while using the remainder to represent all "large" counts.
 *
 * To satisfy this requirement for 192 buckets, we can do the following:
 * Let buckets 0-166 represent distinct counts of [0, 166]
 * Let buckets 166 to 192 represent all remaining counts up to RANK_POSITION_MAX_COUNT_LOG using log2 bucketing.
 */
#define RANK_POSITION_MAX_COUNT_LOG 32
#define RANK_POSITION_LOG_BUCKETS_BEGIN ((RANK_POSITION_TABLE_SIZE - 1) - RANK_POSITION_MAX_COUNT_LOG - 1 /* == 158 */)
#define RANK_POSITION_DISTINCT_COUNT_CUTOFF (RANK_POSITION_LOG_BUCKETS_BEGIN + ZSTD_highbit32(RANK_POSITION_LOG_BUCKETS_BEGIN) /* == 166 */)

/* Return the appropriate bucket index for a given count. See definition of
 * RANK_POSITION_DISTINCT_COUNT_CUTOFF for explanation of bucketing strategy.
 */
static U32 HUF_getIndex(U32 const count) {
    return (count < RANK_POSITION_DISTINCT_COUNT_CUTOFF)
        ? count
        : ZSTD_highbit32(count) + RANK_POSITION_LOG_BUCKETS_BEGIN;
}

/* Helper swap function for HUF_quickSortPartition() */
static void HUF_swapNodes(nodeElt* a, nodeElt* b) {
	nodeElt tmp = *a;
	*a = *b;
	*b = tmp;
}

/* Returns 0 if the huffNode array is not sorted by descending count */
MEM_STATIC int HUF_isSorted(nodeElt huffNode[], U32 const maxSymbolValue1) {
    U32 i;
    for (i = 1; i < maxSymbolValue1; ++i) {
        if (huffNode[i].count > huffNode[i-1].count) {
            return 0;
        }
    }
    return 1;
}

/* Insertion sort by descending order */
HINT_INLINE void HUF_insertionSort(nodeElt huffNode[], int const low, int const high) {
    int i;
    int const size = high-low+1;
    huffNode += low;
    for (i = 1; i < size; ++i) {
        nodeElt const key = huffNode[i];
        int j = i - 1;
        while (j >= 0 && huffNode[j].count < key.count) {
            huffNode[j + 1] = huffNode[j];
            j--;
        }
        huffNode[j + 1] = key;
    }
}

/* Pivot helper function for quicksort. */
static int HUF_quickSortPartition(nodeElt arr[], int const low, int const high) {
    /* Simply select rightmost element as pivot. "Better" selectors like
     * median-of-three don't experimentally appear to have any benefit.
     */
    U32 const pivot = arr[high].count;
    int i = low - 1;
    int j = low;
    for ( ; j < high; j++) {
        if (arr[j].count > pivot) {
            i++;
            HUF_swapNodes(&arr[i], &arr[j]);
        }
    }
    HUF_swapNodes(&arr[i + 1], &arr[high]);
    return i + 1;
}

/* Classic quicksort by descending with partially iterative calls
 * to reduce worst case callstack size.
 */
static void HUF_simpleQuickSort(nodeElt arr[], int low, int high) {
    int const kInsertionSortThreshold = 8;
    if (high - low < kInsertionSortThreshold) {
        HUF_insertionSort(arr, low, high);
        return;
    }
    while (low < high) {
        int const idx = HUF_quickSortPartition(arr, low, high);
        if (idx - low < high - idx) {
            HUF_simpleQuickSort(arr, low, idx - 1);
            low = idx + 1;
        } else {
            HUF_simpleQuickSort(arr, idx + 1, high);
            high = idx - 1;
        }
    }
}

/**
 * HUF_sort():
 * Sorts the symbols [0, maxSymbolValue] by count[symbol] in decreasing order.
 * This is a typical bucket sorting strategy that uses either quicksort or insertion sort to sort each bucket.
 *
 * @param[out] huffNode       Sorted symbols by decreasing count. Only members `.count` and `.byte` are filled.
 *                            Must have (maxSymbolValue + 1) entries.
 * @param[in]  count          Histogram of the symbols.
 * @param[in]  maxSymbolValue Maximum symbol value.
 * @param      rankPosition   This is a scratch workspace. Must have RANK_POSITION_TABLE_SIZE entries.
 */
static void HUF_sort(nodeElt huffNode[], const unsigned count[], U32 const maxSymbolValue, rankPos rankPosition[]) {
    U32 n;
    U32 const maxSymbolValue1 = maxSymbolValue+1;

    /* Compute base and set curr to base.
     * For symbol s let lowerRank = HUF_getIndex(count[n]) and rank = lowerRank + 1.
     * See HUF_getIndex to see bucketing strategy.
     * We attribute each symbol to lowerRank's base value, because we want to know where
     * each rank begins in the output, so for rank R we want to count ranks R+1 and above.
     */
    ZSTD_memset(rankPosition, 0, sizeof(*rankPosition) * RANK_POSITION_TABLE_SIZE);
    for (n = 0; n < maxSymbolValue1; ++n) {
        U32 lowerRank = HUF_getIndex(count[n]);
        assert(lowerRank < RANK_POSITION_TABLE_SIZE - 1);
        rankPosition[lowerRank].base++;
    }

    assert(rankPosition[RANK_POSITION_TABLE_SIZE - 1].base == 0);
    /* Set up the rankPosition table */
    for (n = RANK_POSITION_TABLE_SIZE - 1; n > 0; --n) {
        rankPosition[n-1].base += rankPosition[n].base;
        rankPosition[n-1].curr = rankPosition[n-1].base;
    }

    /* Insert each symbol into their appropriate bucket, setting up rankPosition table. */
    for (n = 0; n < maxSymbolValue1; ++n) {
        U32 const c = count[n];
        U32 const r = HUF_getIndex(c) + 1;
        U32 const pos = rankPosition[r].curr++;
        assert(pos < maxSymbolValue1);
        huffNode[pos].count = c;
        huffNode[pos].byte  = (BYTE)n;
    }

    /* Sort each bucket. */
    for (n = RANK_POSITION_DISTINCT_COUNT_CUTOFF; n < RANK_POSITION_TABLE_SIZE - 1; ++n) {
        int const bucketSize = rankPosition[n].curr - rankPosition[n].base;
        U32 const bucketStartIdx = rankPosition[n].base;
        if (bucketSize > 1) {
            assert(bucketStartIdx < maxSymbolValue1);
            HUF_simpleQuickSort(huffNode + bucketStartIdx, 0, bucketSize-1);
        }
    }

    assert(HUF_isSorted(huffNode, maxSymbolValue1));
}


/** HUF_buildCTable_wksp() :
 *  Same as HUF_buildCTable(), but using externally allocated scratch buffer.
 *  `workSpace` must be aligned on 4-bytes boundaries, and be at least as large as sizeof(HUF_buildCTable_wksp_tables).
 */
#define STARTNODE (HUF_SYMBOLVALUE_MAX+1)

/* HUF_buildTree():
 * Takes the huffNode array sorted by HUF_sort() and builds an unlimited-depth Huffman tree.
 *
 * @param huffNode        The array sorted by HUF_sort(). Builds the Huffman tree in this array.
 * @param maxSymbolValue  The maximum symbol value.
 * @return                The smallest node in the Huffman tree (by count).
 */
static int HUF_buildTree(nodeElt* huffNode, U32 maxSymbolValue)
{
    nodeElt* const huffNode0 = huffNode - 1;
    int nonNullRank;
    int lowS, lowN;
    int nodeNb = STARTNODE;
    int n, nodeRoot;
    DEBUGLOG(5, "HUF_buildTree (alphabet size = %u)", maxSymbolValue + 1);
    /* init for parents */
    nonNullRank = (int)maxSymbolValue;
    while(huffNode[nonNullRank].count == 0) nonNullRank--;
    lowS = nonNullRank; nodeRoot = nodeNb + lowS - 1; lowN = nodeNb;
    huffNode[nodeNb].count = huffNode[lowS].count + huffNode[lowS-1].count;
    huffNode[lowS].parent = huffNode[lowS-1].parent = (U16)nodeNb;
    nodeNb++; lowS-=2;
    for (n=nodeNb; n<=nodeRoot; n++) huffNode[n].count = (U32)(1U<<30);
    huffNode0[0].count = (U32)(1U<<31);  /* fake entry, strong barrier */

    /* create parents */
    while (nodeNb <= nodeRoot) {
        int const n1 = (huffNode[lowS].count < huffNode[lowN].count) ? lowS-- : lowN++;
        int const n2 = (huffNode[lowS].count < huffNode[lowN].count) ? lowS-- : lowN++;
        huffNode[nodeNb].count = huffNode[n1].count + huffNode[n2].count;
        huffNode[n1].parent = huffNode[n2].parent = (U16)nodeNb;
        nodeNb++;
    }

    /* distribute weights (unlimited tree height) */
    huffNode[nodeRoot].nbBits = 0;
    for (n=nodeRoot-1; n>=STARTNODE; n--)
        huffNode[n].nbBits = huffNode[ huffNode[n].parent ].nbBits + 1;
    for (n=0; n<=nonNullRank; n++)
        huffNode[n].nbBits = huffNode[ huffNode[n].parent ].nbBits + 1;

    DEBUGLOG(6, "Initial distribution of bits completed (%zu sorted symbols)", showHNodeBits(huffNode, maxSymbolValue+1));

    return nonNullRank;
}

/**
 * HUF_buildCTableFromTree():
 * Build the CTable given the Huffman tree in huffNode.
 *
 * @param[out] CTable         The output Huffman CTable.
 * @param      huffNode       The Huffman tree.
 * @param      nonNullRank    The last and smallest node in the Huffman tree.
 * @param      maxSymbolValue The maximum symbol value.
 * @param      maxNbBits      The exact maximum number of bits used in the Huffman tree.
 */
static void HUF_buildCTableFromTree(HUF_CElt* CTable, nodeElt const* huffNode, int nonNullRank, U32 maxSymbolValue, U32 maxNbBits)
{
    HUF_CElt* const ct = CTable + 1;
    /* fill result into ctable (val, nbBits) */
    int n;
    U16 nbPerRank[HUF_TABLELOG_MAX+1] = {0};
    U16 valPerRank[HUF_TABLELOG_MAX+1] = {0};
    int const alphabetSize = (int)(maxSymbolValue + 1);
    for (n=0; n<=nonNullRank; n++)
        nbPerRank[huffNode[n].nbBits]++;
    /* determine starting value per rank */
    {   U16 min = 0;
        for (n=(int)maxNbBits; n>0; n--) {
            valPerRank[n] = min;      /* get starting value within each rank */
            min += nbPerRank[n];
            min >>= 1;
    }   }
    for (n=0; n<alphabetSize; n++)
        HUF_setNbBits(ct + huffNode[n].byte, huffNode[n].nbBits);   /* push nbBits per symbol, symbol order */
    for (n=0; n<alphabetSize; n++)
        HUF_setValue(ct + n, valPerRank[HUF_getNbBits(ct[n])]++);   /* assign value within rank, symbol order */

    HUF_writeCTableHeader(CTable, maxNbBits, maxSymbolValue);
}

size_t
HUF_buildCTable_wksp(HUF_CElt* CTable, const unsigned* count, U32 maxSymbolValue, U32 maxNbBits,
                     void* workSpace, size_t wkspSize)
{
    HUF_buildCTable_wksp_tables* const wksp_tables =
        (HUF_buildCTable_wksp_tables*)HUF_alignUpWorkspace(workSpace, &wkspSize, ZSTD_ALIGNOF(U32));
    nodeElt* const huffNode0 = wksp_tables->huffNodeTbl;
    nodeElt* const huffNode = huffNode0+1;
    int nonNullRank;

    HUF_STATIC_ASSERT(HUF_CTABLE_WORKSPACE_SIZE == sizeof(HUF_buildCTable_wksp_tables));

    DEBUGLOG(5, "HUF_buildCTable_wksp (alphabet size = %u)", maxSymbolValue+1);

    /* safety checks */
    if (wkspSize < sizeof(HUF_buildCTable_wksp_tables))
        return ERROR(workSpace_tooSmall);
    if (maxNbBits == 0) maxNbBits = HUF_TABLELOG_DEFAULT;
    if (maxSymbolValue > HUF_SYMBOLVALUE_MAX)
        return ERROR(maxSymbolValue_tooLarge);
    ZSTD_memset(huffNode0, 0, sizeof(huffNodeTable));

    /* sort, decreasing order */
    HUF_sort(huffNode, count, maxSymbolValue, wksp_tables->rankPosition);
    DEBUGLOG(6, "sorted symbols completed (%zu symbols)", showHNodeSymbols(huffNode, maxSymbolValue+1));

    /* build tree */
    nonNullRank = HUF_buildTree(huffNode, maxSymbolValue);

    /* determine and enforce maxTableLog */
    maxNbBits = HUF_setMaxHeight(huffNode, (U32)nonNullRank, maxNbBits);
    if (maxNbBits > HUF_TABLELOG_MAX) return ERROR(GENERIC);   /* check fit into table */

    HUF_buildCTableFromTree(CTable, huffNode, nonNullRank, maxSymbolValue, maxNbBits);

    return maxNbBits;
}

size_t HUF_estimateCompressedSize(const HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue)
{
    HUF_CElt const* ct = CTable + 1;
    size_t nbBits = 0;
    int s;
    for (s = 0; s <= (int)maxSymbolValue; ++s) {
        nbBits += HUF_getNbBits(ct[s]) * count[s];
    }
    return nbBits >> 3;
}

int HUF_validateCTable(const HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue) {
    HUF_CTableHeader header = HUF_readCTableHeader(CTable);
    HUF_CElt const* ct = CTable + 1;
    int bad = 0;
    int s;

    assert(header.tableLog <= HUF_TABLELOG_ABSOLUTEMAX);

    if (header.maxSymbolValue < maxSymbolValue)
        return 0;

    for (s = 0; s <= (int)maxSymbolValue; ++s) {
        bad |= (count[s] != 0) & (HUF_getNbBits(ct[s]) == 0);
    }
    return !bad;
}

size_t HUF_compressBound(size_t size) { return HUF_COMPRESSBOUND(size); }

/** HUF_CStream_t:
 * Huffman uses its own BIT_CStream_t implementation.
 * There are three major differences from BIT_CStream_t:
 *   1. HUF_addBits() takes a HUF_CElt (size_t) which is
 *      the pair (nbBits, value) in the format:
 *      format:
 *        - Bits [0, 4)            = nbBits
 *        - Bits [4, 64 - nbBits)  = 0
 *        - Bits [64 - nbBits, 64) = value
 *   2. The bitContainer is built from the upper bits and
 *      right shifted. E.g. to add a new value of N bits
 *      you right shift the bitContainer by N, then or in
 *      the new value into the N upper bits.
 *   3. The bitstream has two bit containers. You can add
 *      bits to the second container and merge them into
 *      the first container.
 */

#define HUF_BITS_IN_CONTAINER (sizeof(size_t) * 8)

typedef struct {
    size_t bitContainer[2];
    size_t bitPos[2];

    BYTE* startPtr;
    BYTE* ptr;
    BYTE* endPtr;
} HUF_CStream_t;

/**! HUF_initCStream():
 * Initializes the bitstream.
 * @returns 0 or an error code.
 */
static size_t HUF_initCStream(HUF_CStream_t* bitC,
                                  void* startPtr, size_t dstCapacity)
{
    ZSTD_memset(bitC, 0, sizeof(*bitC));
    bitC->startPtr = (BYTE*)startPtr;
    bitC->ptr = bitC->startPtr;
    bitC->endPtr = bitC->startPtr + dstCapacity - sizeof(bitC->bitContainer[0]);
    if (dstCapacity <= sizeof(bitC->bitContainer[0])) return ERROR(dstSize_tooSmall);
    return 0;
}

/*! HUF_addBits():
 * Adds the symbol stored in HUF_CElt elt to the bitstream.
 *
 * @param elt   The element we're adding. This is a (nbBits, value) pair.
 *              See the HUF_CStream_t docs for the format.
 * @param idx   Insert into the bitstream at this idx.
 * @param kFast This is a template parameter. If the bitstream is guaranteed
 *              to have at least 4 unused bits after this call it may be 1,
 *              otherwise it must be 0. HUF_addBits() is faster when fast is set.
 */
FORCE_INLINE_TEMPLATE void HUF_addBits(HUF_CStream_t* bitC, HUF_CElt elt, int idx, int kFast)
{
    assert(idx <= 1);
    assert(HUF_getNbBits(elt) <= HUF_TABLELOG_ABSOLUTEMAX);
    /* This is efficient on x86-64 with BMI2 because shrx
     * only reads the low 6 bits of the register. The compiler
     * knows this and elides the mask. When fast is set,
     * every operation can use the same value loaded from elt.
     */
    bitC->bitContainer[idx] >>= HUF_getNbBits(elt);
    bitC->bitContainer[idx] |= kFast ? HUF_getValueFast(elt) : HUF_getValue(elt);
    /* We only read the low 8 bits of bitC->bitPos[idx] so it
     * doesn't matter that the high bits have noise from the value.
     */
    bitC->bitPos[idx] += HUF_getNbBitsFast(elt);
    assert((bitC->bitPos[idx] & 0xFF) <= HUF_BITS_IN_CONTAINER);
    /* The last 4-bits of elt are dirty if fast is set,
     * so we must not be overwriting bits that have already been
     * inserted into the bit container.
     */
#if DEBUGLEVEL >= 1
    {
        size_t const nbBits = HUF_getNbBits(elt);
        size_t const dirtyBits = nbBits == 0 ? 0 : ZSTD_highbit32((U32)nbBits) + 1;
        (void)dirtyBits;
        /* Middle bits are 0. */
        assert(((elt >> dirtyBits) << (dirtyBits + nbBits)) == 0);
        /* We didn't overwrite any bits in the bit container. */
        assert(!kFast || (bitC->bitPos[idx] & 0xFF) <= HUF_BITS_IN_CONTAINER);
        (void)dirtyBits;
    }
#endif
}

FORCE_INLINE_TEMPLATE void HUF_zeroIndex1(HUF_CStream_t* bitC)
{
    bitC->bitContainer[1] = 0;
    bitC->bitPos[1] = 0;
}

/*! HUF_mergeIndex1() :
 * Merges the bit container @ index 1 into the bit container @ index 0
 * and zeros the bit container @ index 1.
 */
FORCE_INLINE_TEMPLATE void HUF_mergeIndex1(HUF_CStream_t* bitC)
{
    assert((bitC->bitPos[1] & 0xFF) < HUF_BITS_IN_CONTAINER);
    bitC->bitContainer[0] >>= (bitC->bitPos[1] & 0xFF);
    bitC->bitContainer[0] |= bitC->bitContainer[1];
    bitC->bitPos[0] += bitC->bitPos[1];
    assert((bitC->bitPos[0] & 0xFF) <= HUF_BITS_IN_CONTAINER);
}

/*! HUF_flushBits() :
* Flushes the bits in the bit container @ index 0.
*
* @post bitPos will be < 8.
* @param kFast If kFast is set then we must know a-priori that
*              the bit container will not overflow.
*/
FORCE_INLINE_TEMPLATE void HUF_flushBits(HUF_CStream_t* bitC, int kFast)
{
    /* The upper bits of bitPos are noisy, so we must mask by 0xFF. */
    size_t const nbBits = bitC->bitPos[0] & 0xFF;
    size_t const nbBytes = nbBits >> 3;
    /* The top nbBits bits of bitContainer are the ones we need. */
    size_t const bitContainer = bitC->bitContainer[0] >> (HUF_BITS_IN_CONTAINER - nbBits);
    /* Mask bitPos to account for the bytes we consumed. */
    bitC->bitPos[0] &= 7;
    assert(nbBits > 0);
    assert(nbBits <= sizeof(bitC->bitContainer[0]) * 8);
    assert(bitC->ptr <= bitC->endPtr);
    MEM_writeLEST(bitC->ptr, bitContainer);
    bitC->ptr += nbBytes;
    assert(!kFast || bitC->ptr <= bitC->endPtr);
    if (!kFast && bitC->ptr > bitC->endPtr) bitC->ptr = bitC->endPtr;
    /* bitContainer doesn't need to be modified because the leftover
     * bits are already the top bitPos bits. And we don't care about
     * noise in the lower values.
     */
}

/*! HUF_endMark()
 * @returns The Huffman stream end mark: A 1-bit value = 1.
 */
static HUF_CElt HUF_endMark(void)
{
    HUF_CElt endMark;
    HUF_setNbBits(&endMark, 1);
    HUF_setValue(&endMark, 1);
    return endMark;
}

/*! HUF_closeCStream() :
 *  @return Size of CStream, in bytes,
 *          or 0 if it could not fit into dstBuffer */
static size_t HUF_closeCStream(HUF_CStream_t* bitC)
{
    HUF_addBits(bitC, HUF_endMark(), /* idx */ 0, /* kFast */ 0);
    HUF_flushBits(bitC, /* kFast */ 0);
    {
        size_t const nbBits = bitC->bitPos[0] & 0xFF;
        if (bitC->ptr >= bitC->endPtr) return 0; /* overflow detected */
        return (size_t)(bitC->ptr - bitC->startPtr) + (nbBits > 0);
    }
}

FORCE_INLINE_TEMPLATE void
HUF_encodeSymbol(HUF_CStream_t* bitCPtr, U32 symbol, const HUF_CElt* CTable, int idx, int fast)
{
    HUF_addBits(bitCPtr, CTable[symbol], idx, fast);
}

FORCE_INLINE_TEMPLATE void
HUF_compress1X_usingCTable_internal_body_loop(HUF_CStream_t* bitC,
                                   const BYTE* ip, size_t srcSize,
                                   const HUF_CElt* ct,
                                   int kUnroll, int kFastFlush, int kLastFast)
{
    /* Join to kUnroll */
    int n = (int)srcSize;
    int rem = n % kUnroll;
    if (rem > 0) {
        for (; rem > 0; --rem) {
            HUF_encodeSymbol(bitC, ip[--n], ct, 0, /* fast */ 0);
        }
        HUF_flushBits(bitC, kFastFlush);
    }
    assert(n % kUnroll == 0);

    /* Join to 2 * kUnroll */
    if (n % (2 * kUnroll)) {
        int u;
        for (u = 1; u < kUnroll; ++u) {
            HUF_encodeSymbol(bitC, ip[n - u], ct, 0, 1);
        }
        HUF_encodeSymbol(bitC, ip[n - kUnroll], ct, 0, kLastFast);
        HUF_flushBits(bitC, kFastFlush);
        n -= kUnroll;
    }
    assert(n % (2 * kUnroll) == 0);

    for (; n>0; n-= 2 * kUnroll) {
        /* Encode kUnroll symbols into the bitstream @ index 0. */
        int u;
        for (u = 1; u < kUnroll; ++u) {
            HUF_encodeSymbol(bitC, ip[n - u], ct, /* idx */ 0, /* fast */ 1);
        }
        HUF_encodeSymbol(bitC, ip[n - kUnroll], ct, /* idx */ 0, /* fast */ kLastFast);
        HUF_flushBits(bitC, kFastFlush);
        /* Encode kUnroll symbols into the bitstream @ index 1.
         * This allows us to start filling the bit container
         * without any data dependencies.
         */
        HUF_zeroIndex1(bitC);
        for (u = 1; u < kUnroll; ++u) {
            HUF_encodeSymbol(bitC, ip[n - kUnroll - u], ct, /* idx */ 1, /* fast */ 1);
        }
        HUF_encodeSymbol(bitC, ip[n - kUnroll - kUnroll], ct, /* idx */ 1, /* fast */ kLastFast);
        /* Merge bitstream @ index 1 into the bitstream @ index 0 */
        HUF_mergeIndex1(bitC);
        HUF_flushBits(bitC, kFastFlush);
    }
    assert(n == 0);

}

/**
 * Returns a tight upper bound on the output space needed by Huffman
 * with 8 bytes buffer to handle over-writes. If the output is at least
 * this large we don't need to do bounds checks during Huffman encoding.
 */
static size_t HUF_tightCompressBound(size_t srcSize, size_t tableLog)
{
    return ((srcSize * tableLog) >> 3) + 8;
}


FORCE_INLINE_TEMPLATE size_t
HUF_compress1X_usingCTable_internal_body(void* dst, size_t dstSize,
                                   const void* src, size_t srcSize,
                                   const HUF_CElt* CTable)
{
    U32 const tableLog = HUF_readCTableHeader(CTable).tableLog;
    HUF_CElt const* ct = CTable + 1;
    const BYTE* ip = (const BYTE*) src;
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstSize;
    HUF_CStream_t bitC;

    /* init */
    if (dstSize < 8) return 0;   /* not enough space to compress */
    { BYTE* op = ostart;
      size_t const initErr = HUF_initCStream(&bitC, op, (size_t)(oend-op));
      if (HUF_isError(initErr)) return 0; }

    if (dstSize < HUF_tightCompressBound(srcSize, (size_t)tableLog) || tableLog > 11)
        HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ MEM_32bits() ? 2 : 4, /* kFast */ 0, /* kLastFast */ 0);
    else {
        if (MEM_32bits()) {
            switch (tableLog) {
            case 11:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 2, /* kFastFlush */ 1, /* kLastFast */ 0);
                break;
            case 10: ZSTD_FALLTHROUGH;
            case 9: ZSTD_FALLTHROUGH;
            case 8:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 2, /* kFastFlush */ 1, /* kLastFast */ 1);
                break;
            case 7: ZSTD_FALLTHROUGH;
            default:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 3, /* kFastFlush */ 1, /* kLastFast */ 1);
                break;
            }
        } else {
            switch (tableLog) {
            case 11:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 5, /* kFastFlush */ 1, /* kLastFast */ 0);
                break;
            case 10:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 5, /* kFastFlush */ 1, /* kLastFast */ 1);
                break;
            case 9:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 6, /* kFastFlush */ 1, /* kLastFast */ 0);
                break;
            case 8:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 7, /* kFastFlush */ 1, /* kLastFast */ 0);
                break;
            case 7:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 8, /* kFastFlush */ 1, /* kLastFast */ 0);
                break;
            case 6: ZSTD_FALLTHROUGH;
            default:
                HUF_compress1X_usingCTable_internal_body_loop(&bitC, ip, srcSize, ct, /* kUnroll */ 9, /* kFastFlush */ 1, /* kLastFast */ 1);
                break;
            }
        }
    }
    assert(bitC.ptr <= bitC.endPtr);

    return HUF_closeCStream(&bitC);
}

#if DYNAMIC_BMI2

static BMI2_TARGET_ATTRIBUTE size_t
HUF_compress1X_usingCTable_internal_bmi2(void* dst, size_t dstSize,
                                   const void* src, size_t srcSize,
                                   const HUF_CElt* CTable)
{
    return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

static size_t
HUF_compress1X_usingCTable_internal_default(void* dst, size_t dstSize,
                                      const void* src, size_t srcSize,
                                      const HUF_CElt* CTable)
{
    return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

static size_t
HUF_compress1X_usingCTable_internal(void* dst, size_t dstSize,
                              const void* src, size_t srcSize,
                              const HUF_CElt* CTable, const int flags)
{
    if (flags & HUF_flags_bmi2) {
        return HUF_compress1X_usingCTable_internal_bmi2(dst, dstSize, src, srcSize, CTable);
    }
    return HUF_compress1X_usingCTable_internal_default(dst, dstSize, src, srcSize, CTable);
}

#else

static size_t
HUF_compress1X_usingCTable_internal(void* dst, size_t dstSize,
                              const void* src, size_t srcSize,
                              const HUF_CElt* CTable, const int flags)
{
    (void)flags;
    return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

#endif

size_t HUF_compress1X_usingCTable(void* dst, size_t dstSize, const void* src, size_t srcSize, const HUF_CElt* CTable, int flags)
{
    return HUF_compress1X_usingCTable_internal(dst, dstSize, src, srcSize, CTable, flags);
}

static size_t
HUF_compress4X_usingCTable_internal(void* dst, size_t dstSize,
                              const void* src, size_t srcSize,
                              const HUF_CElt* CTable, int flags)
{
    size_t const segmentSize = (srcSize+3)/4;   /* first 3 segments */
    const BYTE* ip = (const BYTE*) src;
    const BYTE* const iend = ip + srcSize;
    BYTE* const ostart = (BYTE*) dst;
    BYTE* const oend = ostart + dstSize;
    BYTE* op = ostart;

    if (dstSize < 6 + 1 + 1 + 1 + 8) return 0;   /* minimum space to compress successfully */
    if (srcSize < 12) return 0;   /* no saving possible : too small input */
    op += 6;   /* jumpTable */

    assert(op <= oend);
    {   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, (size_t)(oend-op), ip, segmentSize, CTable, flags) );
        if (cSize == 0 || cSize > 65535) return 0;
        MEM_writeLE16(ostart, (U16)cSize);
        op += cSize;
    }

    ip += segmentSize;
    assert(op <= oend);
    {   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, (size_t)(oend-op), ip, segmentSize, CTable, flags) );
        if (cSize == 0 || cSize > 65535) return 0;
        MEM_writeLE16(ostart+2, (U16)cSize);
        op += cSize;
    }

    ip += segmentSize;
    assert(op <= oend);
    {   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, (size_t)(oend-op), ip, segmentSize, CTable, flags) );
        if (cSize == 0 || cSize > 65535) return 0;
        MEM_writeLE16(ostart+4, (U16)cSize);
        op += cSize;
    }

    ip += segmentSize;
    assert(op <= oend);
    assert(ip <= iend);
    {   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, (size_t)(oend-op), ip, (size_t)(iend-ip), CTable, flags) );
        if (cSize == 0 || cSize > 65535) return 0;
        op += cSize;
    }

    return (size_t)(op-ostart);
}

size_t HUF_compress4X_usingCTable(void* dst, size_t dstSize, const void* src, size_t srcSize, const HUF_CElt* CTable, int flags)
{
    return HUF_compress4X_usingCTable_internal(dst, dstSize, src, srcSize, CTable, flags);
}

typedef enum { HUF_singleStream, HUF_fourStreams } HUF_nbStreams_e;

static size_t HUF_compressCTable_internal(
                BYTE* const ostart, BYTE* op, BYTE* const oend,
                const void* src, size_t srcSize,
                HUF_nbStreams_e nbStreams, const HUF_CElt* CTable, const int flags)
{
    size_t const cSize = (nbStreams==HUF_singleStream) ?
                         HUF_compress1X_usingCTable_internal(op, (size_t)(oend - op), src, srcSize, CTable, flags) :
                         HUF_compress4X_usingCTable_internal(op, (size_t)(oend - op), src, srcSize, CTable, flags);
    if (HUF_isError(cSize)) { return cSize; }
    if (cSize==0) { return 0; }   /* uncompressible */
    op += cSize;
    /* check compressibility */
    assert(op >= ostart);
    if ((size_t)(op-ostart) >= srcSize-1) { return 0; }
    return (size_t)(op-ostart);
}

typedef struct {
    unsigned count[HUF_SYMBOLVALUE_MAX + 1];
    HUF_CElt CTable[HUF_CTABLE_SIZE_ST(HUF_SYMBOLVALUE_MAX)];
    union {
        HUF_buildCTable_wksp_tables buildCTable_wksp;
        HUF_WriteCTableWksp writeCTable_wksp;
        U32 hist_wksp[HIST_WKSP_SIZE_U32];
    } wksps;
} HUF_compress_tables_t;

#define SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE 4096
#define SUSPECT_INCOMPRESSIBLE_SAMPLE_RATIO 10  /* Must be >= 2 */

unsigned HUF_cardinality(const unsigned* count, unsigned maxSymbolValue)
{
    unsigned cardinality = 0;
    unsigned i;

    for (i = 0; i < maxSymbolValue + 1; i++) {
        if (count[i] != 0) cardinality += 1;
    }

    return cardinality;
}

unsigned HUF_minTableLog(unsigned symbolCardinality)
{
    U32 minBitsSymbols = ZSTD_highbit32(symbolCardinality) + 1;
    return minBitsSymbols;
}

unsigned HUF_optimalTableLog(
            unsigned maxTableLog,
            size_t srcSize,
            unsigned maxSymbolValue,
            void* workSpace, size_t wkspSize,
            HUF_CElt* table,
      const unsigned* count,
            int flags)
{
    assert(srcSize > 1); /* Not supported, RLE should be used instead */
    assert(wkspSize >= sizeof(HUF_buildCTable_wksp_tables));

    if (!(flags & HUF_flags_optimalDepth)) {
        /* cheap evaluation, based on FSE */
        return FSE_optimalTableLog_internal(maxTableLog, srcSize, maxSymbolValue, 1);
    }

    {   BYTE* dst = (BYTE*)workSpace + sizeof(HUF_WriteCTableWksp);
        size_t dstSize = wkspSize - sizeof(HUF_WriteCTableWksp);
        size_t hSize, newSize;
        const unsigned symbolCardinality = HUF_cardinality(count, maxSymbolValue);
        const unsigned minTableLog = HUF_minTableLog(symbolCardinality);
        size_t optSize = ((size_t) ~0) - 1;
        unsigned optLog = maxTableLog, optLogGuess;

        DEBUGLOG(6, "HUF_optimalTableLog: probing huf depth (srcSize=%zu)", srcSize);

        /* Search until size increases */
        for (optLogGuess = minTableLog; optLogGuess <= maxTableLog; optLogGuess++) {
            DEBUGLOG(7, "checking for huffLog=%u", optLogGuess);

            {   size_t maxBits = HUF_buildCTable_wksp(table, count, maxSymbolValue, optLogGuess, workSpace, wkspSize);
                if (ERR_isError(maxBits)) continue;

                if (maxBits < optLogGuess && optLogGuess > minTableLog) break;

                hSize = HUF_writeCTable_wksp(dst, dstSize, table, maxSymbolValue, (U32)maxBits, workSpace, wkspSize);
            }

            if (ERR_isError(hSize)) continue;

            newSize = HUF_estimateCompressedSize(table, count, maxSymbolValue) + hSize;

            if (newSize > optSize + 1) {
                break;
            }

            if (newSize < optSize) {
                optSize = newSize;
                optLog = optLogGuess;
            }
        }
        assert(optLog <= HUF_TABLELOG_MAX);
        return optLog;
    }
}

/* HUF_compress_internal() :
 * `workSpace_align4` must be aligned on 4-bytes boundaries,
 * and occupies the same space as a table of HUF_WORKSPACE_SIZE_U64 unsigned */
static size_t
HUF_compress_internal (void* dst, size_t dstSize,
                 const void* src, size_t srcSize,
                       unsigned maxSymbolValue, unsigned huffLog,
                       HUF_nbStreams_e nbStreams,
                       void* workSpace, size_t wkspSize,
                       HUF_CElt* oldHufTable, HUF_repeat* repeat, int flags)
{
    HUF_compress_tables_t* const table = (HUF_compress_tables_t*)HUF_alignUpWorkspace(workSpace, &wkspSize, ZSTD_ALIGNOF(size_t));
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstSize;
    BYTE* op = ostart;

    DEBUGLOG(5, "HUF_compress_internal (srcSize=%zu)", srcSize);
    HUF_STATIC_ASSERT(sizeof(*table) + HUF_WORKSPACE_MAX_ALIGNMENT <= HUF_WORKSPACE_SIZE);

    /* checks & inits */
    if (wkspSize < sizeof(*table)) return ERROR(workSpace_tooSmall);
    if (!srcSize) return 0;  /* Uncompressed */
    if (!dstSize) return 0;  /* cannot fit anything within dst budget */
    if (srcSize > HUF_BLOCKSIZE_MAX) return ERROR(srcSize_wrong);   /* current block size limit */
    if (huffLog > HUF_TABLELOG_MAX) return ERROR(tableLog_tooLarge);
    if (maxSymbolValue > HUF_SYMBOLVALUE_MAX) return ERROR(maxSymbolValue_tooLarge);
    if (!maxSymbolValue) maxSymbolValue = HUF_SYMBOLVALUE_MAX;
    if (!huffLog) huffLog = HUF_TABLELOG_DEFAULT;

    /* Heuristic : If old table is valid, use it for small inputs */
    if ((flags & HUF_flags_preferRepeat) && repeat && *repeat == HUF_repeat_valid) {
        return HUF_compressCTable_internal(ostart, op, oend,
                                           src, srcSize,
                                           nbStreams, oldHufTable, flags);
    }

    /* If uncompressible data is suspected, do a smaller sampling first */
    DEBUG_STATIC_ASSERT(SUSPECT_INCOMPRESSIBLE_SAMPLE_RATIO >= 2);
    if ((flags & HUF_flags_suspectUncompressible) && srcSize >= (SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE * SUSPECT_INCOMPRESSIBLE_SAMPLE_RATIO)) {
        size_t largestTotal = 0;
        DEBUGLOG(5, "input suspected incompressible : sampling to check");
        {   unsigned maxSymbolValueBegin = maxSymbolValue;
            CHECK_V_F(largestBegin, HIST_count_simple (table->count, &maxSymbolValueBegin, (const BYTE*)src, SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE) );
            largestTotal += largestBegin;
        }
        {   unsigned maxSymbolValueEnd = maxSymbolValue;
            CHECK_V_F(largestEnd, HIST_count_simple (table->count, &maxSymbolValueEnd, (const BYTE*)src + srcSize - SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE, SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE) );
            largestTotal += largestEnd;
        }
        if (largestTotal <= ((2 * SUSPECT_INCOMPRESSIBLE_SAMPLE_SIZE) >> 7)+4) return 0;   /* heuristic : probably not compressible enough */
    }

    /* Scan input and build symbol stats */
    {   CHECK_V_F(largest, HIST_count_wksp (table->count, &maxSymbolValue, (const BYTE*)src, srcSize, table->wksps.hist_wksp, sizeof(table->wksps.hist_wksp)) );
        if (largest == srcSize) { *ostart = ((const BYTE*)src)[0]; return 1; }   /* single symbol, rle */
        if (largest <= (srcSize >> 7)+4) return 0;   /* heuristic : probably not compressible enough */
    }
    DEBUGLOG(6, "histogram detail completed (%zu symbols)", showU32(table->count, maxSymbolValue+1));

    /* Check validity of previous table */
    if ( repeat
      && *repeat == HUF_repeat_check
      && !HUF_validateCTable(oldHufTable, table->count, maxSymbolValue)) {
        *repeat = HUF_repeat_none;
    }
    /* Heuristic : use existing table for small inputs */
    if ((flags & HUF_flags_preferRepeat) && repeat && *repeat != HUF_repeat_none) {
        return HUF_compressCTable_internal(ostart, op, oend,
                                           src, srcSize,
                                           nbStreams, oldHufTable, flags);
    }

    /* Build Huffman Tree */
    huffLog = HUF_optimalTableLog(huffLog, srcSize, maxSymbolValue, &table->wksps, sizeof(table->wksps), table->CTable, table->count, flags);
    {   size_t const maxBits = HUF_buildCTable_wksp(table->CTable, table->count,
                                            maxSymbolValue, huffLog,
                                            &table->wksps.buildCTable_wksp, sizeof(table->wksps.buildCTable_wksp));
        CHECK_F(maxBits);
        huffLog = (U32)maxBits;
        DEBUGLOG(6, "bit distribution completed (%zu symbols)", showCTableBits(table->CTable + 1, maxSymbolValue+1));
    }

    /* Write table description header */
    {   CHECK_V_F(hSize, HUF_writeCTable_wksp(op, dstSize, table->CTable, maxSymbolValue, huffLog,
                                              &table->wksps.writeCTable_wksp, sizeof(table->wksps.writeCTable_wksp)) );
        /* Check if using previous huffman table is beneficial */
        if (repeat && *repeat != HUF_repeat_none) {
            size_t const oldSize = HUF_estimateCompressedSize(oldHufTable, table->count, maxSymbolValue);
            size_t const newSize = HUF_estimateCompressedSize(table->CTable, table->count, maxSymbolValue);
            if (oldSize <= hSize + newSize || hSize + 12 >= srcSize) {
                return HUF_compressCTable_internal(ostart, op, oend,
                                                   src, srcSize,
                                                   nbStreams, oldHufTable, flags);
        }   }

        /* Use the new huffman table */
        if (hSize + 12ul >= srcSize) { return 0; }
        op += hSize;
        if (repeat) { *repeat = HUF_repeat_none; }
        if (oldHufTable)
            ZSTD_memcpy(oldHufTable, table->CTable, sizeof(table->CTable));  /* Save new table */
    }
    return HUF_compressCTable_internal(ostart, op, oend,
                                       src, srcSize,
                                       nbStreams, table->CTable, flags);
}

size_t HUF_compress1X_repeat (void* dst, size_t dstSize,
                      const void* src, size_t srcSize,
                      unsigned maxSymbolValue, unsigned huffLog,
                      void* workSpace, size_t wkspSize,
                      HUF_CElt* hufTable, HUF_repeat* repeat, int flags)
{
    DEBUGLOG(5, "HUF_compress1X_repeat (srcSize = %zu)", srcSize);
    return HUF_compress_internal(dst, dstSize, src, srcSize,
                                 maxSymbolValue, huffLog, HUF_singleStream,
                                 workSpace, wkspSize, hufTable,
                                 repeat, flags);
}

/* HUF_compress4X_repeat():
 * compress input using 4 streams.
 * consider skipping quickly
 * reuse an existing huffman compression table */
size_t HUF_compress4X_repeat (void* dst, size_t dstSize,
                      const void* src, size_t srcSize,
                      unsigned maxSymbolValue, unsigned huffLog,
                      void* workSpace, size_t wkspSize,
                      HUF_CElt* hufTable, HUF_repeat* repeat, int flags)
{
    DEBUGLOG(5, "HUF_compress4X_repeat (srcSize = %zu)", srcSize);
    return HUF_compress_internal(dst, dstSize, src, srcSize,
                                 maxSymbolValue, huffLog, HUF_fourStreams,
                                 workSpace, wkspSize,
                                 hufTable, repeat, flags);
}


/* ******************************************************************
 * FSE : Finite State Entropy encoder
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 *  You can contact the author at :
 *  - FSE source repository : https://github.com/Cyan4973/FiniteStateEntropy
 *  - Public forum : https://groups.google.com/forum/#!forum/lz4c
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
****************************************************************** */

/* **************************************************************
*  Includes
****************************************************************/
#include "../common/compiler.h"
#include "../common/mem.h"        /* U32, U16, etc. */
#include "../common/debug.h"      /* assert, DEBUGLOG */
#include "hist.h"       /* HIST_count_wksp */
#include "../common/bitstream.h"
#define FSE_STATIC_LINKING_ONLY
#include "../common/fse.h"
#include "../common/error_private.h"
#define ZSTD_DEPS_NEED_MALLOC
#define ZSTD_DEPS_NEED_MATH64
#include "../common/zstd_deps.h"  /* ZSTD_malloc, ZSTD_free, ZSTD_memcpy, ZSTD_memset */
#include "../common/bits.h" /* ZSTD_highbit32 */


/* **************************************************************
*  Error Management
****************************************************************/
#define FSE_isError ERR_isError


/* **************************************************************
*  Templates
****************************************************************/
/*
  designed to be included
  for type-specific functions (template emulation in C)
  Objective is to write these functions only once, for improved maintenance
*/

/* safety checks */
#ifndef FSE_FUNCTION_EXTENSION
#  error "FSE_FUNCTION_EXTENSION must be defined"
#endif
#ifndef FSE_FUNCTION_TYPE
#  error "FSE_FUNCTION_TYPE must be defined"
#endif

/* Function names */
#define FSE_CAT(X,Y) X##Y
#define FSE_FUNCTION_NAME(X,Y) FSE_CAT(X,Y)
#define FSE_TYPE_NAME(X,Y) FSE_CAT(X,Y)


/* Function templates */

/* FSE_buildCTable_wksp() :
 * Same as FSE_buildCTable(), but using an externally allocated scratch buffer (`workSpace`).
 * wkspSize should be sized to handle worst case situation, which is `1<<max_tableLog * sizeof(FSE_FUNCTION_TYPE)`
 * workSpace must also be properly aligned with FSE_FUNCTION_TYPE requirements
 */
size_t FSE_buildCTable_wksp(FSE_CTable* ct,
                      const short* normalizedCounter, unsigned maxSymbolValue, unsigned tableLog,
                            void* workSpace, size_t wkspSize)
{
    U32 const tableSize = 1 << tableLog;
    U32 const tableMask = tableSize - 1;
    void* const ptr = ct;
    U16* const tableU16 = ( (U16*) ptr) + 2;
    void* const FSCT = ((U32*)ptr) + 1 /* header */ + (tableLog ? tableSize>>1 : 1) ;
    FSE_symbolCompressionTransform* const symbolTT = (FSE_symbolCompressionTransform*) (FSCT);
    U32 const step = FSE_TABLESTEP(tableSize);
    U32 const maxSV1 = maxSymbolValue+1;

    U16* cumul = (U16*)workSpace;   /* size = maxSV1 */
    FSE_FUNCTION_TYPE* const tableSymbol = (FSE_FUNCTION_TYPE*)(cumul + (maxSV1+1));  /* size = tableSize */

    U32 highThreshold = tableSize-1;

    assert(((size_t)workSpace & 1) == 0);  /* Must be 2 bytes-aligned */
    if (FSE_BUILD_CTABLE_WORKSPACE_SIZE(maxSymbolValue, tableLog) > wkspSize) return ERROR(tableLog_tooLarge);
    /* CTable header */
    tableU16[-2] = (U16) tableLog;
    tableU16[-1] = (U16) maxSymbolValue;
    assert(tableLog < 16);   /* required for threshold strategy to work */

    /* For explanations on how to distribute symbol values over the table :
     * https://fastcompression.blogspot.fr/2014/02/fse-distributing-symbol-values.html */

     #ifdef __clang_analyzer__
     ZSTD_memset(tableSymbol, 0, sizeof(*tableSymbol) * tableSize);   /* useless initialization, just to keep scan-build happy */
     #endif

    /* symbol start positions */
    {   U32 u;
        cumul[0] = 0;
        for (u=1; u <= maxSV1; u++) {
            if (normalizedCounter[u-1]==-1) {  /* Low proba symbol */
                cumul[u] = cumul[u-1] + 1;
                tableSymbol[highThreshold--] = (FSE_FUNCTION_TYPE)(u-1);
            } else {
                assert(normalizedCounter[u-1] >= 0);
                cumul[u] = cumul[u-1] + (U16)normalizedCounter[u-1];
                assert(cumul[u] >= cumul[u-1]);  /* no overflow */
        }   }
        cumul[maxSV1] = (U16)(tableSize+1);
    }

    /* Spread symbols */
    if (highThreshold == tableSize - 1) {
        /* Case for no low prob count symbols. Lay down 8 bytes at a time
         * to reduce branch misses since we are operating on a small block
         */
        BYTE* const spread = tableSymbol + tableSize; /* size = tableSize + 8 (may write beyond tableSize) */
        {   U64 const add = 0x0101010101010101ull;
            size_t pos = 0;
            U64 sv = 0;
            U32 s;
            for (s=0; s<maxSV1; ++s, sv += add) {
                int i;
                int const n = normalizedCounter[s];
                MEM_write64(spread + pos, sv);
                for (i = 8; i < n; i += 8) {
                    MEM_write64(spread + pos + i, sv);
                }
                assert(n>=0);
                pos += (size_t)n;
            }
        }
        /* Spread symbols across the table. Lack of lowprob symbols means that
         * we don't need variable sized inner loop, so we can unroll the loop and
         * reduce branch misses.
         */
        {   size_t position = 0;
            size_t s;
            size_t const unroll = 2; /* Experimentally determined optimal unroll */
            assert(tableSize % unroll == 0); /* FSE_MIN_TABLELOG is 5 */
            for (s = 0; s < (size_t)tableSize; s += unroll) {
                size_t u;
                for (u = 0; u < unroll; ++u) {
                    size_t const uPosition = (position + (u * step)) & tableMask;
                    tableSymbol[uPosition] = spread[s + u];
                }
                position = (position + (unroll * step)) & tableMask;
            }
            assert(position == 0);   /* Must have initialized all positions */
        }
    } else {
        U32 position = 0;
        U32 symbol;
        for (symbol=0; symbol<maxSV1; symbol++) {
            int nbOccurrences;
            int const freq = normalizedCounter[symbol];
            for (nbOccurrences=0; nbOccurrences<freq; nbOccurrences++) {
                tableSymbol[position] = (FSE_FUNCTION_TYPE)symbol;
                position = (position + step) & tableMask;
                while (position > highThreshold)
                    position = (position + step) & tableMask;   /* Low proba area */
        }   }
        assert(position==0);  /* Must have initialized all positions */
    }

    /* Build table */
    {   U32 u; for (u=0; u<tableSize; u++) {
        FSE_FUNCTION_TYPE s = tableSymbol[u];   /* note : static analyzer may not understand tableSymbol is properly initialized */
        tableU16[cumul[s]++] = (U16) (tableSize+u);   /* TableU16 : sorted by symbol order; gives next state value */
    }   }

    /* Build Symbol Transformation Table */
    {   unsigned total = 0;
        unsigned s;
        for (s=0; s<=maxSymbolValue; s++) {
            switch (normalizedCounter[s])
            {
            case  0:
                /* filling nonetheless, for compatibility with FSE_getMaxNbBits() */
                symbolTT[s].deltaNbBits = ((tableLog+1) << 16) - (1<<tableLog);
                break;

            case -1:
            case  1:
                symbolTT[s].deltaNbBits = (tableLog << 16) - (1<<tableLog);
                assert(total <= INT_MAX);
                symbolTT[s].deltaFindState = (int)(total - 1);
                total ++;
                break;
            default :
                assert(normalizedCounter[s] > 1);
                {   U32 const maxBitsOut = tableLog - ZSTD_highbit32 ((U32)normalizedCounter[s]-1);
                    U32 const minStatePlus = (U32)normalizedCounter[s] << maxBitsOut;
                    symbolTT[s].deltaNbBits = (maxBitsOut << 16) - minStatePlus;
                    symbolTT[s].deltaFindState = (int)(total - (unsigned)normalizedCounter[s]);
                    total +=  (unsigned)normalizedCounter[s];
    }   }   }   }

#if 0  /* debug : symbol costs */
    DEBUGLOG(5, "\n --- table statistics : ");
    {   U32 symbol;
        for (symbol=0; symbol<=maxSymbolValue; symbol++) {
            DEBUGLOG(5, "%3u: w=%3i,   maxBits=%u, fracBits=%.2f",
                symbol, normalizedCounter[symbol],
                FSE_getMaxNbBits(symbolTT, symbol),
                (double)FSE_bitCost(symbolTT, tableLog, symbol, 8) / 256);
    }   }
#endif

    return 0;
}



#ifndef FSE_COMMONDEFS_ONLY

/*-**************************************************************
*  FSE NCount encoding
****************************************************************/
size_t FSE_NCountWriteBound(unsigned maxSymbolValue, unsigned tableLog)
{
    size_t const maxHeaderSize = (((maxSymbolValue+1) * tableLog
                                   + 4 /* bitCount initialized at 4 */
                                   + 2 /* first two symbols may use one additional bit each */) / 8)
                                   + 1 /* round up to whole nb bytes */
                                   + 2 /* additional two bytes for bitstream flush */;
    return maxSymbolValue ? maxHeaderSize : FSE_NCOUNTBOUND;  /* maxSymbolValue==0 ? use default */
}

static size_t
FSE_writeNCount_generic (void* header, size_t headerBufferSize,
                   const short* normalizedCounter, unsigned maxSymbolValue, unsigned tableLog,
                         unsigned writeIsSafe)
{
    BYTE* const ostart = (BYTE*) header;
    BYTE* out = ostart;
    BYTE* const oend = ostart + headerBufferSize;
    int nbBits;
    const int tableSize = 1 << tableLog;
    int remaining;
    int threshold;
    U32 bitStream = 0;
    int bitCount = 0;
    unsigned symbol = 0;
    unsigned const alphabetSize = maxSymbolValue + 1;
    int previousIs0 = 0;

    /* Table Size */
    bitStream += (tableLog-FSE_MIN_TABLELOG) << bitCount;
    bitCount  += 4;

    /* Init */
    remaining = tableSize+1;   /* +1 for extra accuracy */
    threshold = tableSize;
    nbBits = (int)tableLog+1;

    while ((symbol < alphabetSize) && (remaining>1)) {  /* stops at 1 */
        if (previousIs0) {
            unsigned start = symbol;
            while ((symbol < alphabetSize) && !normalizedCounter[symbol]) symbol++;
            if (symbol == alphabetSize) break;   /* incorrect distribution */
            while (symbol >= start+24) {
                start+=24;
                bitStream += 0xFFFFU << bitCount;
                if ((!writeIsSafe) && (out > oend-2))
                    return ERROR(dstSize_tooSmall);   /* Buffer overflow */
                out[0] = (BYTE) bitStream;
                out[1] = (BYTE)(bitStream>>8);
                out+=2;
                bitStream>>=16;
            }
            while (symbol >= start+3) {
                start+=3;
                bitStream += 3U << bitCount;
                bitCount += 2;
            }
            bitStream += (symbol-start) << bitCount;
            bitCount += 2;
            if (bitCount>16) {
                if ((!writeIsSafe) && (out > oend - 2))
                    return ERROR(dstSize_tooSmall);   /* Buffer overflow */
                out[0] = (BYTE)bitStream;
                out[1] = (BYTE)(bitStream>>8);
                out += 2;
                bitStream >>= 16;
                bitCount -= 16;
        }   }
        {   int count = normalizedCounter[symbol++];
            int const max = (2*threshold-1) - remaining;
            remaining -= count < 0 ? -count : count;
            count++;   /* +1 for extra accuracy */
            if (count>=threshold)
                count += max;   /* [0..max[ [max..threshold[ (...) [threshold+max 2*threshold[ */
            bitStream += (U32)count << bitCount;
            bitCount  += nbBits;
            bitCount  -= (count<max);
            previousIs0  = (count==1);
            if (remaining<1) return ERROR(GENERIC);
            while (remaining<threshold) { nbBits--; threshold>>=1; }
        }
        if (bitCount>16) {
            if ((!writeIsSafe) && (out > oend - 2))
                return ERROR(dstSize_tooSmall);   /* Buffer overflow */
            out[0] = (BYTE)bitStream;
            out[1] = (BYTE)(bitStream>>8);
            out += 2;
            bitStream >>= 16;
            bitCount -= 16;
    }   }

    if (remaining != 1)
        return ERROR(GENERIC);  /* incorrect normalized distribution */
    assert(symbol <= alphabetSize);

    /* flush remaining bitStream */
    if ((!writeIsSafe) && (out > oend - 2))
        return ERROR(dstSize_tooSmall);   /* Buffer overflow */
    out[0] = (BYTE)bitStream;
    out[1] = (BYTE)(bitStream>>8);
    out+= (bitCount+7) /8;

    assert(out >= ostart);
    return (size_t)(out-ostart);
}


size_t FSE_writeNCount (void* buffer, size_t bufferSize,
                  const short* normalizedCounter, unsigned maxSymbolValue, unsigned tableLog)
{
    if (tableLog > FSE_MAX_TABLELOG) return ERROR(tableLog_tooLarge);   /* Unsupported */
    if (tableLog < FSE_MIN_TABLELOG) return ERROR(GENERIC);   /* Unsupported */

    if (bufferSize < FSE_NCountWriteBound(maxSymbolValue, tableLog))
        return FSE_writeNCount_generic(buffer, bufferSize, normalizedCounter, maxSymbolValue, tableLog, 0);

    return FSE_writeNCount_generic(buffer, bufferSize, normalizedCounter, maxSymbolValue, tableLog, 1 /* write in buffer is safe */);
}


/*-**************************************************************
*  FSE Compression Code
****************************************************************/

/* provides the minimum logSize to safely represent a distribution */
static unsigned FSE_minTableLog(size_t srcSize, unsigned maxSymbolValue)
{
    U32 minBitsSrc = ZSTD_highbit32((U32)(srcSize)) + 1;
    U32 minBitsSymbols = ZSTD_highbit32(maxSymbolValue) + 2;
    U32 minBits = minBitsSrc < minBitsSymbols ? minBitsSrc : minBitsSymbols;
    assert(srcSize > 1); /* Not supported, RLE should be used instead */
    return minBits;
}

unsigned FSE_optimalTableLog_internal(unsigned maxTableLog, size_t srcSize, unsigned maxSymbolValue, unsigned minus)
{
    U32 maxBitsSrc = ZSTD_highbit32((U32)(srcSize - 1)) - minus;
    U32 tableLog = maxTableLog;
    U32 minBits = FSE_minTableLog(srcSize, maxSymbolValue);
    assert(srcSize > 1); /* Not supported, RLE should be used instead */
    if (tableLog==0) tableLog = FSE_DEFAULT_TABLELOG;
    if (maxBitsSrc < tableLog) tableLog = maxBitsSrc;   /* Accuracy can be reduced */
    if (minBits > tableLog) tableLog = minBits;   /* Need a minimum to safely represent all symbol values */
    if (tableLog < FSE_MIN_TABLELOG) tableLog = FSE_MIN_TABLELOG;
    if (tableLog > FSE_MAX_TABLELOG) tableLog = FSE_MAX_TABLELOG;
    return tableLog;
}

unsigned FSE_optimalTableLog(unsigned maxTableLog, size_t srcSize, unsigned maxSymbolValue)
{
    return FSE_optimalTableLog_internal(maxTableLog, srcSize, maxSymbolValue, 2);
}

/* Secondary normalization method.
   To be used when primary method fails. */

static size_t FSE_normalizeM2(short* norm, U32 tableLog, const unsigned* count, size_t total, U32 maxSymbolValue, short lowProbCount)
{
    short const NOT_YET_ASSIGNED = -2;
    U32 s;
    U32 distributed = 0;
    U32 ToDistribute;

    /* Init */
    U32 const lowThreshold = (U32)(total >> tableLog);
    U32 lowOne = (U32)((total * 3) >> (tableLog + 1));

    for (s=0; s<=maxSymbolValue; s++) {
        if (count[s] == 0) {
            norm[s]=0;
            continue;
        }
        if (count[s] <= lowThreshold) {
            norm[s] = lowProbCount;
            distributed++;
            total -= count[s];
            continue;
        }
        if (count[s] <= lowOne) {
            norm[s] = 1;
            distributed++;
            total -= count[s];
            continue;
        }

        norm[s]=NOT_YET_ASSIGNED;
    }
    ToDistribute = (1 << tableLog) - distributed;

    if (ToDistribute == 0)
        return 0;

    if ((total / ToDistribute) > lowOne) {
        /* risk of rounding to zero */
        lowOne = (U32)((total * 3) / (ToDistribute * 2));
        for (s=0; s<=maxSymbolValue; s++) {
            if ((norm[s] == NOT_YET_ASSIGNED) && (count[s] <= lowOne)) {
                norm[s] = 1;
                distributed++;
                total -= count[s];
                continue;
        }   }
        ToDistribute = (1 << tableLog) - distributed;
    }

    if (distributed == maxSymbolValue+1) {
        /* all values are pretty poor;
           probably incompressible data (should have already been detected);
           find max, then give all remaining points to max */
        U32 maxV = 0, maxC = 0;
        for (s=0; s<=maxSymbolValue; s++)
            if (count[s] > maxC) { maxV=s; maxC=count[s]; }
        norm[maxV] += (short)ToDistribute;
        return 0;
    }

    if (total == 0) {
        /* all of the symbols were low enough for the lowOne or lowThreshold */
        for (s=0; ToDistribute > 0; s = (s+1)%(maxSymbolValue+1))
            if (norm[s] > 0) { ToDistribute--; norm[s]++; }
        return 0;
    }

    {   U64 const vStepLog = 62 - tableLog;
        U64 const mid = (1ULL << (vStepLog-1)) - 1;
        U64 const rStep = ZSTD_div64((((U64)1<<vStepLog) * ToDistribute) + mid, (U32)total);   /* scale on remaining */
        U64 tmpTotal = mid;
        for (s=0; s<=maxSymbolValue; s++) {
            if (norm[s]==NOT_YET_ASSIGNED) {
                U64 const end = tmpTotal + (count[s] * rStep);
                U32 const sStart = (U32)(tmpTotal >> vStepLog);
                U32 const sEnd = (U32)(end >> vStepLog);
                U32 const weight = sEnd - sStart;
                if (weight < 1)
                    return ERROR(GENERIC);
                norm[s] = (short)weight;
                tmpTotal = end;
    }   }   }

    return 0;
}

size_t FSE_normalizeCount (short* normalizedCounter, unsigned tableLog,
                           const unsigned* count, size_t total,
                           unsigned maxSymbolValue, unsigned useLowProbCount)
{
    /* Sanity checks */
    if (tableLog==0) tableLog = FSE_DEFAULT_TABLELOG;
    if (tableLog < FSE_MIN_TABLELOG) return ERROR(GENERIC);   /* Unsupported size */
    if (tableLog > FSE_MAX_TABLELOG) return ERROR(tableLog_tooLarge);   /* Unsupported size */
    if (tableLog < FSE_minTableLog(total, maxSymbolValue)) return ERROR(GENERIC);   /* Too small tableLog, compression potentially impossible */

    {   static U32 const rtbTable[] = {     0, 473195, 504333, 520860, 550000, 700000, 750000, 830000 };
        short const lowProbCount = useLowProbCount ? -1 : 1;
        U64 const scale = 62 - tableLog;
        U64 const step = ZSTD_div64((U64)1<<62, (U32)total);   /* <== here, one division ! */
        U64 const vStep = 1ULL<<(scale-20);
        int stillToDistribute = 1<<tableLog;
        unsigned s;
        unsigned largest=0;
        short largestP=0;
        U32 lowThreshold = (U32)(total >> tableLog);

        for (s=0; s<=maxSymbolValue; s++) {
            if (count[s] == total) return 0;   /* rle special case */
            if (count[s] == 0) { normalizedCounter[s]=0; continue; }
            if (count[s] <= lowThreshold) {
                normalizedCounter[s] = lowProbCount;
                stillToDistribute--;
            } else {
                short proba = (short)((count[s]*step) >> scale);
                if (proba<8) {
                    U64 restToBeat = vStep * rtbTable[proba];
                    proba += (count[s]*step) - ((U64)proba<<scale) > restToBeat;
                }
                if (proba > largestP) { largestP=proba; largest=s; }
                normalizedCounter[s] = proba;
                stillToDistribute -= proba;
        }   }
        if (-stillToDistribute >= (normalizedCounter[largest] >> 1)) {
            /* corner case, need another normalization method */
            size_t const errorCode = FSE_normalizeM2(normalizedCounter, tableLog, count, total, maxSymbolValue, lowProbCount);
            if (FSE_isError(errorCode)) return errorCode;
        }
        else normalizedCounter[largest] += (short)stillToDistribute;
    }

#if 0
    {   /* Print Table (debug) */
        U32 s;
        U32 nTotal = 0;
        for (s=0; s<=maxSymbolValue; s++)
            RAWLOG(2, "%3i: %4i \n", s, normalizedCounter[s]);
        for (s=0; s<=maxSymbolValue; s++)
            nTotal += abs(normalizedCounter[s]);
        if (nTotal != (1U<<tableLog))
            RAWLOG(2, "Warning !!! Total == %u != %u !!!", nTotal, 1U<<tableLog);
        getchar();
    }
#endif

    return tableLog;
}

/* fake FSE_CTable, for rle input (always same symbol) */
size_t FSE_buildCTable_rle (FSE_CTable* ct, BYTE symbolValue)
{
    void* ptr = ct;
    U16* tableU16 = ( (U16*) ptr) + 2;
    void* FSCTptr = (U32*)ptr + 2;
    FSE_symbolCompressionTransform* symbolTT = (FSE_symbolCompressionTransform*) FSCTptr;

    /* header */
    tableU16[-2] = (U16) 0;
    tableU16[-1] = (U16) symbolValue;

    /* Build table */
    tableU16[0] = 0;
    tableU16[1] = 0;   /* just in case */

    /* Build Symbol Transformation Table */
    symbolTT[symbolValue].deltaNbBits = 0;
    symbolTT[symbolValue].deltaFindState = 0;

    return 0;
}


static size_t FSE_compress_usingCTable_generic (void* dst, size_t dstSize,
                           const void* src, size_t srcSize,
                           const FSE_CTable* ct, const unsigned fast)
{
    const BYTE* const istart = (const BYTE*) src;
    const BYTE* const iend = istart + srcSize;
    const BYTE* ip=iend;

    BIT_CStream_t bitC;
    FSE_CState_t CState1, CState2;

    /* init */
    if (srcSize <= 2) return 0;
    { size_t const initError = BIT_initCStream(&bitC, dst, dstSize);
      if (FSE_isError(initError)) return 0; /* not enough space available to write a bitstream */ }

#define FSE_FLUSHBITS(s)  (fast ? BIT_flushBitsFast(s) : BIT_flushBits(s))

    if (srcSize & 1) {
        FSE_initCState2(&CState1, ct, *--ip);
        FSE_initCState2(&CState2, ct, *--ip);
        FSE_encodeSymbol(&bitC, &CState1, *--ip);
        FSE_FLUSHBITS(&bitC);
    } else {
        FSE_initCState2(&CState2, ct, *--ip);
        FSE_initCState2(&CState1, ct, *--ip);
    }

    /* join to mod 4 */
    srcSize -= 2;
    if ((sizeof(bitC.bitContainer)*8 > FSE_MAX_TABLELOG*4+7 ) && (srcSize & 2)) {  /* test bit 2 */
        FSE_encodeSymbol(&bitC, &CState2, *--ip);
        FSE_encodeSymbol(&bitC, &CState1, *--ip);
        FSE_FLUSHBITS(&bitC);
    }

    /* 2 or 4 encoding per loop */
    while ( ip>istart ) {

        FSE_encodeSymbol(&bitC, &CState2, *--ip);

        if (sizeof(bitC.bitContainer)*8 < FSE_MAX_TABLELOG*2+7 )   /* this test must be static */
            FSE_FLUSHBITS(&bitC);

        FSE_encodeSymbol(&bitC, &CState1, *--ip);

        if (sizeof(bitC.bitContainer)*8 > FSE_MAX_TABLELOG*4+7 ) {  /* this test must be static */
            FSE_encodeSymbol(&bitC, &CState2, *--ip);
            FSE_encodeSymbol(&bitC, &CState1, *--ip);
        }

        FSE_FLUSHBITS(&bitC);
    }

    FSE_flushCState(&bitC, &CState2);
    FSE_flushCState(&bitC, &CState1);
    return BIT_closeCStream(&bitC);
}

size_t FSE_compress_usingCTable (void* dst, size_t dstSize,
                           const void* src, size_t srcSize,
                           const FSE_CTable* ct)
{
    unsigned const fast = (dstSize >= FSE_BLOCKBOUND(srcSize));

    if (fast)
        return FSE_compress_usingCTable_generic(dst, dstSize, src, srcSize, ct, 1);
    else
        return FSE_compress_usingCTable_generic(dst, dstSize, src, srcSize, ct, 0);
}


size_t FSE_compressBound(size_t size) { return FSE_COMPRESSBOUND(size); }

#endif   /* FSE_COMMONDEFS_ONLY */
