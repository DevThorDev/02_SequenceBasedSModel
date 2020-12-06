# -*- coding: utf-8 -*-
################################################################################
# GenericFunctions.py #
################################################################################

import os, shutil, random, math, zipfile
import numpy

from scipy.stats import binom

import InputVariables
from Constants import *

# ### FUNCTIONS ################################################################
def drawFromPList(lProb):
    if len(lProb) == 0:
        print('Error: Probability vector has size 0!')
        return ERR_CODE
    drIdx = 0
    curSum = 0.0
    drNum = random.random()
    while curSum <= drNum and drIdx < len(lProb):
        curSum += lProb[drIdx]
        drIdx += 1
    return drIdx - 1

def addToCountDict(theDict, drEl):
    if drEl in theDict:
        theDict[drEl] += 1
    else:
        theDict[drEl] = 1

def addListValToDict(theDict, theKey, theVal):
    if theKey in theDict:
        if len(theVal) > 1:
            for idx in range(len(theVal)):
                theDict[theKey][idx] += theVal[idx]
        else:
            theDict[theKey] += theVal
    else:
        theDict[theKey] = theVal

def addToDictOfLists(theDict, theKey, theEl):
    if theKey in theDict:
        theDict[theKey].append(theEl)
    else:
        theDict[theKey] = [theEl]

def convSubSetToBinChain(lA, nE):
    dA = {}
    # convert to chains of 1s and 0s
    for lIdx in range(len(lA)):
        lRespPat = [0]*nE
        n1s = 0
        for k in lA[lIdx]:
            lRespPat[k] = 1
            n1s += 1
        dA[lIdx] = [lRespPat, n1s/nE]
    return dA

def complementOldModDict(theDic):
    for oneKey, oneVal in theDic.items():
        newVal = [oneVal, oneVal]
        theDic[oneKey] = newVal

def complementNewModDict(theDic):
    for oneKey, oneVal in theDic.items():
        n1s = 0
        for el in oneVal[0]:
            if el == 1:
                n1s += 1
        oneVal.append(n1s/len(oneVal[0]))
        theDic[oneKey] = oneVal

def randPartitionRange(numParts, rLen = 1.0):
    lShares = [0.0]*numParts
    for i in range(numParts):
        lShares[i] = random.random()
    sumShares = math.fsum(lShares)
    for i in range(numParts):
        lShares[i] *= rLen/sumShares
    return lShares

def calcRelFitness(l1, l2):
    lenGT = 0
    ovGT = 0
    for i in range(len(l1)):
        if l1[i] == 1 or l2[i] == 1:
            lenGT += 1
            if l1[i] == 1 and l2[i] == 1:
                ovGT += 1
    return [lenGT/len(l1), ovGT/len(l1)]

def calcMean(lNums):
    if len(lNums) > 0:
        return sum(lNums)/float(len(lNums))
    else:
        print('Error: list is empty - cannot calculate arithmetic mean.')
        return 0.0;

def calcRunMean(curVal, meanBef, numValsCur):
    if numValsCur > 0:
        return meanBef + (curVal - meanBef)/numValsCur
    else:
        return 0.0

def calcBinProbs(n, pSuccess):
    lBinomProbs = [0.0]*(n + 1)
    for k in range(n + 1):
        lBinomProbs[k] = binom.pmf(k, n, pSuccess)
    return lBinomProbs

def calcCondBinProbs(n, pSuccess):
    p_Xis0 = (1.0-pSuccess)**n
    lCondBinPrs = calcBinProbs(n, pSuccess)
    lCondBinPrs[0] = 0.0
    for idx in range(1, len(lCondBinPrs)):
        lCondBinPrs[idx] /= (1.0 - p_Xis0)
    return lCondBinPrs

def SUB_G2Fit(seqA1, seqA2, fAC1, fAC2, minGTF, maxGTF, domTp, hetAdM, dltI,
              lmbdDB, fHetAd, limFit):
    fGT = minGTF    # fitness of genotype
    if domTp in [D_SEQ_B_AVPD_W_HA, D_SEQ_B_MAXPD_W_HA]:
        numSeqDiff = compare2Seq(seqA1, seqA2)
        if hetAdM in [HAM_CONST_DELTA, HAM_POSDEP_DELTA]:
            fGT = numSeqDiff*dltI    # overdominance contribution
        if hetAdM == HAM_BAD_ALL_DELTA:
            fGT = numSeqDiff*dltI*min(fAC1, fAC2)   # weaker allele has effect
        if domTp == D_SEQ_B_AVPD_W_HA:
            fGT += round(calcMean([fAC1, fAC2]), ROUND_PREC)
        elif domTp == D_SEQ_B_MAXPD_W_HA:
            fGT += max(fAC1, fAC2)
        else: pass
    else:
        if domTp == D_SEQ_B_AVPD_NO_HA:
            fGT = round(calcMean([fAC1, fAC2]), ROUND_PREC)
        elif domTp == D_SEQ_B_MAXPD_NO_HA:
            fGT = max(fAC1, fAC2)
        elif domTp == D_HET_ADV_DEBOER:
            if seqA1 == seqA2:  # homozygote
                fGT = fAC1
            else:               # heterozygote
                fGT = fAC1 + (1.0 - lmbdDB*fAC1)*fAC2
        elif domTp == D_HET_ADV_MAXPD:
            if seqA1 == seqA2:  # homozygote
                fGT = fAC1
            else:               # heterozygote
                fGT = max(fAC1, fAC2) + fHetAd
        elif domTp == D_HET_ADV_AVPD:
            if seqA1 == seqA2:  # homozygote
                fGT = fAC1
            else:               # heterozygote
                fGT = calcMean([fAC1, fAC2]) + fHetAd
        elif domTp == D_HET_ADV_TAKNEI:
            if seqA1 == seqA2:  # homozygote
                fGT = fAC1
            else:               # heterozygote
                fGT = maxGTF
        else: pass
    fGT = capValue(fGT, minGTF, maxGTF, limFit in [L_FIT_1, L_FIT_MAX])
    return fGT

def getDictKey(theDict, givenVal, lIdx = -1):
    curVal = givenVal
    foundKey = ''
    for oneKey in theDict:
        if lIdx >= 0:
            curVal = theDict[oneKey][lIdx]
        else:
            curVal = theDict[oneKey]
        if curVal == givenVal:
            foundKey = oneKey
            break
    if foundKey == '':
        print('Error: key of value', givenVal, 'not found.')
    return foundKey

def getDictVal(theDict, theKey, lIdx = -1):
    if lIdx >= 0:
        try:
            theVal = theDict[theKey][lIdx]
        except KeyError:
            print('ERROR: theKey =', theKey, '; lIdx =', lIdx)
            theVal = ERR_CODE
    else:
        theVal = theDict[theKey]
    return theVal

def multDictVal(theDict, theKey, lIdx = -1, theMult = 1):
    if lIdx >= 0:
        theDict[theKey][lIdx] *= theMult
    else:
        theDict[theKey] *= theMult

def getListOfDictVals(theDict, lIdx = -1, sortList = False):
    lVals = []
    for oneKey in theDict:
        lVals.append(getDictVal(theDict, oneKey, lIdx))
    if sortList:
        lVals.sort()
    return lVals

def calcRoundValSum(theDict, lIdx = -1, theMult = 1):
    rdValSum = 0
    for oneKey in theDict:
        rdValSum += round(getDictVal(theDict, oneKey, lIdx)*theMult)
    return rdValSum

def getDictNumEvents(theDict, lIdx = -1, theMult = 1, theProb = 0.0):
    dKeyChanges = {}
    elIdx = 0
    cElThr = 0
    roundNVals = calcRoundValSum(theDict, lIdx, theMult)
    numEv = numpy.random.binomial(theMult, theProb)
    if numEv > 0:   # otherwise the dictionary stays empty - no mutation
        lElChange = sorted(random.sample(range(roundNVals), numEv))
        for oneKey in sorted(theDict):
            numEvCurKey = 0
            cElThr += round(getDictVal(theDict, oneKey, lIdx)*theMult)
            while lElChange[elIdx] < cElThr:
                numEvCurKey += 1
                elIdx += 1
                if elIdx == len(lElChange):
                    break
            if numEvCurKey > 0:
                dKeyChanges[oneKey] = numEvCurKey
                if numEvCurKey > getDictVal(theDict, oneKey, lIdx)*theMult:
                    print('Empty Mutations! Number of mutations:', numEvCurKey,
                          'number of alleles only',
                          getDictVal(theDict, oneKey, lIdx)*theMult)
            if elIdx == len(lElChange):
                break
    return dKeyChanges

def drawFromBinom(nVal, pVal, numDraws = None):
    drNum = numpy.random.binomial(nVal, pVal, numDraws)
    return drNum 

def drawFromMultinom(nVal, pList, numDraws = None):
    lDrNum = numpy.random.multinomial(nVal, pList, numDraws).tolist()
    return lDrNum

def drawFromPoisson(lambVal, numDraws = None):
    drNum = numpy.random.poisson(lambVal, numDraws)
    return drNum 

def correctProportions(theDict, lIdx = -1, extraProp = 0.0):
    for oneKey in theDict:
        if extraProp > 0:
            multDictVal(theDict, oneKey, lIdx, 1.0 - extraProp)
        else:
            multDictVal(theDict, oneKey, lIdx, 1.0/(1.0 + extraProp))

def safeRemoveEl(theDict, theKey, lIdx = -1, elSz = 0.0):
    del theDict[theKey]
    if elSz != 0.0:
        correctProportions(theDict, lIdx, -elSz)

def compare2Seq(seq1, seq2):
    numDiff = 0
    for idx in range(min(len(seq1), len(seq2))):
        if seq1[idx] != seq2[idx]:
            numDiff += 1
    return numDiff

def capValue(curVal, minVal, maxVal, capAtMax = True, capAtMin = True):
    newVal = curVal
    if capAtMin:
        newVal = max(newVal, minVal)
    if capAtMax:
        newVal = min(newVal, maxVal)
    return newVal

def SUB_NewExtVal(theDict, curK, kExtV, extV, lIdx = -1, extMax = True):
    curV = getDictVal(theDict, curK, lIdx)
    if extMax:
        if curV >= extV:
            kExtV = curK
            extV = curV
    else:
        if curV <= extV:
            kExtV = curK
            extV = curV
    return kExtV, extV

def getStartKeyVal(theDict, lKeys, lIdx = -1, extMax = True):
    if len(lKeys) > 0:
        startKey, startVal = lKeys[0], getDictVal(theDict, lKeys[0], lIdx)
    else:
        print('ERROR!')
        return ERR_CODE, ERR_CODE
    for oneKey in lKeys[1:]:
        startKey, startVal = SUB_NewExtVal(theDict, oneKey, startKey, startVal,
                                           lIdx, extMax)
    return startKey, startVal

def getKeyOfExtVal(theDict, lKFull, lKUsed, stKey, stVal, lIdx = -1,
                   extMax = True, errStr = ''):
    if len(lKFull) > 0:
        keyExtVal = stKey
        extVal = stVal
        for oneKey in lKFull:
            if oneKey not in lKUsed:
                keyExtVal, extVal = SUB_NewExtVal(theDict, oneKey, keyExtVal,
                                                  extVal, lIdx, extMax)
        lKUsed.append(keyExtVal)
        return keyExtVal
    else:
        print(errStr)
        return ERR_CODE

def addEntryToDVList(theDict, theKey, valToAdd):
    if theKey in theDict:
        valList = theDict[theKey]
        valList.append(valToAdd)
        theDict[theKey] = valList
    else:
        theDict[theKey] = [valToAdd]

def weightedMeanAndVarSimple(vVs, vWs):
    meanVal = numpy.average(vVs, weights = vWs)
    varVal = numpy.dot(vWs, (vVs - meanVal)*(vVs - meanVal))/vWs.sum()
    return meanVal, varVal

def weightedMeanAndSD(vVals, vWts, totNumEls = 0, sampleEst = True):
    theMean, theVar = 0, 0
    if len(vVals) == len(vWts):
        if len(vVals) > 0:
            theMean, theVar = weightedMeanAndVarSimple(vVals, vWts)
            if sampleEst:
                if totNumEls > 1:
                    theVar *= (totNumEls/(totNumEls - 1))
                else:
                    print('Total number of elements only', totNumEls,
                          '- sample estimate for variance not possible.')
    else:
        print('Error: length of values array (', len(vVals),
              ') unequal length of weights array (', len(vWts), ').')       
    return theMean, math.sqrt(theVar)

def getTimeStepFromFileName(nF, p1, p2, p3):
    tS = nF.split('_')[p1]
    return int(tS.split('.')[p2][p3:])

def sortListOfFileNames(lNFiles):
    dTSt = {}
    for nFile in lNFiles:
        dTSt[getTimeStepFromFileName(nFile, -1, 0, 8)] = nFile
    lNFilesSort = ['']*len(lNFiles)
    if len(lNFiles) == len(dTSt):
        ct = 0
        for curTSt in sorted(dTSt):
            lNFilesSort[ct] = dTSt[curTSt]
            ct += 1
    else:
        print('Error: Number of time steps (', len(dTSt),
              ') and files (', len(lNFiles), ') unequal.', sep ='')
    return lNFilesSort

def writeListToFile(pathToFile, tableHeader, dataList):
    outFile = open(pathToFile, 'w')
    outFile.write(tableHeader)
    for oneLine in dataList:
        for oneEl in oneLine[:-1]:
            outFile.write(str(oneEl) + ' ')
        outFile.write(str(oneLine[-1]) + '\n')
    outFile.close()

def zipFolder(nFolder, targetDir, delDir = False):            
    zipObj = zipfile.ZipFile(nFolder + '.zip', 'w', zipfile.ZIP_DEFLATED)
    lRoot = len(targetDir) + 1
    for base, dirs, files in os.walk(targetDir):
        for file in files:
            fn = os.path.join(base, file)
            zipObj.write(fn, fn[lRoot:])
    if delDir:
        shutil.rmtree(targetDir)
