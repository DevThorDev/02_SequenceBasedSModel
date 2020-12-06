# -*- coding: utf-8 -*-
################################################################################
# AnalyseRawData.py #
################################################################################
#!/opt/local/Library/Frameworks/Python.framework/Versions/3.2

import os, sys

import numpy
from scipy.stats import hmean

from math import sqrt
#from string import split

lSubDirs = ['SBMSupport']
for oneSubDir in lSubDirs:
    sys.path.append(sys.path[0] + '/' + oneSubDir)

from Constants import ROUND_PREC
from InputVariables import (iniAllFDist, iniAllProp, newAllFDist, limFitness,
                            dominType, deltaIn, gaussMuIn, gaussSigIn, minFitn,
                            maxFitn, minGTFitn, maxGTFitn, dictDefFit, lostThr,
                            disregThr, numAllThr, numRepetitions,
                            numIterations, numARIter, popSize, numAllI,
                            numAmAcPerAll, hetAdvM, lambdaDeB, fixedHA,
                            modDisp, modRate, modFile, dAllSortM, dFlags,
                            dMutInfo, dSpecAna, dClOutc, dOutput, ND_RESULTS,
                            NF_RATERES, NF_ANARES, NF_SPECRES, NF_AGERES,
                            NF_ALLRES, NF_PHYL, NF_GPRES, NF_SGPRES, NF_TGPRES,
                            NF_HA, NF_DIVRES, NF_REPRES, NF_END)
from SpecialFunctions import (genDictInp, getSeqDiffMatrix, sortAlleles,
                              getNumberOfSeqDiffs)
from GenericFunctions import (addListValToDict, SUB_G2Fit, getListOfDictVals,
                              weightedMeanAndSD, getTimeStepFromFileName,
                              sortListOfFileNames, writeListToFile)

print('AnalyseRawData.py')

def getListOfResultFileNames(cDir, cRp, nDirRes, nFGPRes):
    repOfFile = 0
    lNFRes = []
    lNFResRaw = os.listdir(cDir + '/' + nDirRes)
    for nFRes in lNFResRaw:
        if len(nFRes) >= len(nFGPRes):
            if nFRes[:len(nFGPRes)] == nFGPRes:
                try:
                    repOfFile = int(nFRes.split('_')[1][3:])
                except:
                    print('Error: Cannot extract repetition of file', nFRes)
                if cRp == repOfFile:
                    lNFRes.append(nFRes)
    return lNFRes

def writeHeaderOfRawDataAnalysisFile(dI, cRp, cDir):
    fRes = open(cDir + '/' + dI['nD_Results'] + '/' + dI['nF_AnaRes'] +
                '_Rep' + str(cRp) + dI['nF_End'], 'w')
    fRes.write('TimeStep NumberAllelesCurrent AvNumberOfPairDiffs' +
               ' AvPercPairDiffHet PercHomozygotes PercHeterozygotes' +
               ' AvFitnessHomozygotes AvFitnessHeterozygotes' +
               ' AvHeterozygoteAdvantage AvHeterozygoteAdvantagePerc' +
               ' AvOverdominanceAll AvOverdominanceHet' +
               ' AvOverdominanceHetPerc FitnessRangeAbs FitnessRangePerc' +
               ' deBoerFitnessThreshold UnfitAlleles UnfitAllelesPerc' +
               ' AllFitSDUnweighted AllFitSDWeighted AvSeqDiffAbsPerc' +
               ' AvSeqDiffPropWtPerc PopulationFitness\n')
    fRes.close()

def writeHeaderOfSpecialAnalysisFile(dI, cRp, cDir):
    fRes = open(cDir + '/' + dI['nD_Results'] + '/' + dI['nF_SpecRes'] +
                '_Rep' + str(cRp) + dI['nF_End'], 'w')
    fRes.write('TimeStep NumberAllelesTot' +
               ' NumAlleles1E6 NumAlleles1E5 NumAlleles1E4' +
               ' NumAlleles1E3 NumAlleles1E2' +
               ' NumAlleles10cp NumAlleles25cp NumAlleles100cp' +
               ' NumAlleles250cp NumAlleles1000cp' +
               ' NumAlleles10G NumAlleles25G NumAlleles100G NumAlleles250G' +
               ' NumAlleles1000G' +
               ' TotRangePerc TotSD TotSDWt TotDeBThr TotNUnfAll' +
               ' RangePerc1E6T SD1E6T SDWt1E6T DeBThr1E6T NUnfAll1E6T' +
               ' RangePerc1E5T SD1E5T SDWt1E5T DeBThr1E5T NUnfAll1E5T' +
               ' RangePerc1E4T SD1E4T SDWt1E4T DeBThr1E4T NUnfAll1E4T' +
               ' RangePerc1E3T SD1E3T SDWt1E3T DeBThr1E3T NUnfAll1E3T' +
               ' RangePerc1E2T SD1E2T SDWt1E2T DeBThr1E2T NUnfAll1E2T' +
               ' RangePerc10cpT SD10cpT SDWt10cpT DeBThr10cpT NUnfAll10cpT' +
               ' RangePerc25cpT SD25cpT SDWt25cpT DeBThr25cpT NUnfAll25cpT' +
               ' RangePerc100cpT SD100cpT SDWt100cpT DeBThr100cpT NUnfAll100cpT' +
               ' RangePerc250cpT SD250cpT SDWt250cpT DeBThr250cpT NUnfAll250cpT' +
               ' RangePerc1000cpT SD1000cpT SDWt1000cpT DeBThr1000cpT NUnfAll1000cpT' +
               ' RangePerc10GT SD10GT SDWt10GT DeBThr10GT NUnfAll10GT' +
               ' RangePerc25GT SD25GT SDWt25GT DeBThr25GT NUnfAll25GT' +
               ' RangePerc100GT SD100GT SDWt100GT DeBThr100GT NUnfAll100GT' +
               ' RangePerc250GT SD250GT SDWt250GT DeBThr250GT NUnfAll250GT' +
               ' RangePerc1000GT SD1000GT SDWt1000GT DeBThr1000GT NUnfAll1000GT' +
               '\n')
    fRes.close()

def writeHeaderOfAgeClassAnalysisFile(dI, cRp, cDir):
    fRes = open(cDir + '/' + dI['nD_Results'] + '/' + dI['nF_AgeRes'] +
                '_Rep' + str(cRp) + dI['nF_End'], 'w')
    fRes.write('class numAllelesAgeCl shareAllelesAgeCl numAllelesFitCl' +
               ' shareAllelesFitCl numAllelesFreqCl shareAllelesFreqCl\n')
    fRes.close()

def writeHeaderOfAlleleComparisonFile(dI, cRp, cDir):
    fRes = open(cDir + '/' + dI['nD_Results'] + '/' + dI['nF_AllRes'] +
                '_Rep' + str(cRp) + dI['nF_End'], 'w')
    fRes.write('TimeStep AlleleStatus Allele AlleleFitness AlleleProportion' +
               ' MarginalFitness Emerged MutationType SequenceDiffAbsPerc' +
               ' SequenceDiffPropWtPerc\n')
    fRes.close()

def addSeqDiffEntriesToDictGP(dI, dGP, nFGP, cDir):
    print('Adding sequence difference entries to genepool dictionary', nFGP)
    dSeqDiffGP = getSeqDiffMatrix(dGP)
    dVals = {}
    for curSeq in dGP:
        numSeqDiffPerc, wtSeqDiffPerc = \
            getNumberOfSeqDiffs(dI, dGP, dSeqDiffGP, list(dGP.keys()), curSeq)
        dVals[curSeq] = [numSeqDiffPerc, wtSeqDiffPerc]
    fGP = open(cDir + '/' + dI['nD_Results'] + '/' + nFGP, 'r')
    lFGP = fGP.readlines()
    fGP.close()
    fGP = open(cDir + '/' + dI['nD_Results'] + '/' + nFGP, 'w')
    fGP.write(lFGP[0][:-1] + ' SequenceDiffAbsPerc SequenceDiffPropWtPerc' +
              lFGP[0][-1:])
    for oneLine in lFGP[1:]:
        curSeq = oneLine.split()[2]
        oneLine = oneLine[:-1] + ' ' + str(dVals[curSeq][0]) + ' ' + \
            str(dVals[curSeq][1]) + oneLine[-1:]
        fGP.write(oneLine)
    fGP.close()

def parseGPResFile(dI, nGPFile, cDir):
    seqDiffInfo = True
    dGP = {}
    fGP = open(cDir + '/' + dI['nD_Results'] + '/' + nGPFile, 'r')
    lFGP = fGP.readlines()
    for oneLine in lFGP[1:]:
        lVals = oneLine.split()
        if len(lVals) in [6, 8, 10]:
            if len(lVals) == 6: # file does not contain seq. difference info
                dGP[lVals[2]] = [int(float(lVals[1])), float(lVals[3]),
                                 float(lVals[4]), float(lVals[5])]
                seqDiffInfo = False
            elif len(lVals) == 8:   # file contains add. seq. difference info
                dGP[lVals[2]] = [int(float(lVals[1])), float(lVals[3]),
                                 float(lVals[4]), float(lVals[5]),
                                 float(lVals[6]), float(lVals[7])]
            else:               # file also contains "time step emerged" info
                dGP[lVals[2]] = [int(float(lVals[1])), float(lVals[3]),
                                 float(lVals[4]), float(lVals[5]),
                                 int(float(lVals[6])), int(float(lVals[7])),
                                 float(lVals[8]), float(lVals[9])]
        else:
            print('Error in parsing file ', nGPFile, '.', sep ='')
            print('Line ', lFGP.index(oneLine), ' contains ', len(lVals),
                  ' entries.', sep ='')
    if not seqDiffInfo:
        addSeqDiffEntriesToDictGP(dI, dGP, nGPFile, cDir)
    fGP.close()
    return dGP

def SUB_calcPropFit(dI, dGP, s1, s2):
    # proportion of genotypes carrying this allele
    propCGenTp = dGP[s1][2]*dGP[s2][2]
    # fitness of genotypes carrying this allele
    fitCGenTp = SUB_G2Fit(s1, s2, dGP[s1][1], dGP[s2][1], dI['minGTFit'],
                          dI['maxGTFit'], dI['domType'], dI['hetAdvM'],
                          dI['deltaIn'], dI['lambdaDeB'], dI['fixedHA'],
                          dI['limFit'])
    return propCGenTp, fitCGenTp

def SUB_calcAvNumDiff(dGP, dSDGP, s1, s2, avND):
    if s1 == s2:    # homozygote - no sequence differences
        return avND
    numSD = 0
    keyV1, keyV2 = (dGP[s1][0], dGP[s2][0]), (dGP[s2][0], dGP[s1][0])
    if keyV1 in dSDGP:
        numSD = dSDGP[keyV1]
    elif keyV2 in dSDGP:
        numSD = dSDGP[keyV2]
    else:
        print('ERROR:', keyV1, 'and', keyV2, 'not in seq. diff. dictionary.')
    avND += numSD*dGP[s1][2]*dGP[s2][2]
    return avND

def calcHetAdvFromRawData(dI, dGP, dSqDfGP):
    dHetAdv = {'propOfHom': 0.0, 'propOfHet': 0.0, 'avFitHom': 0.0,
               'avFitHet': 0.0, 'avNumDiff': 0.0, 'avOverdom': 0.0}
    propOfHom, propOfHet = 0.0, 0.0
    avFitHom, avFitHet, avOverdom = 0.0, 0.0, 0.0
    avNDiff = 0.0
    for curSeq in dGP:
        for othSeq in dGP:
            if othSeq == curSeq:
                propCHom, fitCHom = SUB_calcPropFit(dI, dGP, curSeq, curSeq)
                propOfHom += propCHom
                avFitHom += propCHom*fitCHom
            else:
                propCHet, fitCHet = SUB_calcPropFit(dI, dGP, curSeq, othSeq)
                propOfHet += propCHet
                avFitHet += propCHet*fitCHet
                avOverdom += propCHet*(fitCHet - max(dGP[curSeq][1],
                                                     dGP[othSeq][1]))
            avNDiff = SUB_calcAvNumDiff(dGP, dSqDfGP, curSeq, othSeq, avNDiff)
    dHetAdv['propOfHom'] = propOfHom
    dHetAdv['propOfHet'] = propOfHet
    if propOfHom > 0.0:
        dHetAdv['avFitHom'] = avFitHom/propOfHom
    if propOfHet > 0.0:
        dHetAdv['avFitHet'] = avFitHet/propOfHet
    dHetAdv['avOverdom'] = avOverdom
    dHetAdv['avNumDiff'] = avNDiff
    if round(propOfHom + propOfHet, ROUND_PREC - 4) != 1.0:
        print('Proportions of homozygotes and heterozygotes sum to',
              propOfHom + propOfHet)
    return dHetAdv

def calcAlleleFitVals(dGP, lSeqs, numAll = 0):
    lFits = [0.0]*len(lSeqs)
    lProps = [0.0]*len(lSeqs)
    minF, maxF, sdAFit, sdWtAFit, dBThr = 1.0, 0.0, 0.0, 0.0, 0.0
    ct, nUnfAll = 0, 0
    for curSeq in lSeqs:
        curFit = dGP[curSeq][1]
        if curFit < minF:
            minF = curFit
        if curFit > maxF:
            maxF = curFit
        lFits[ct] = curFit
        lProps[ct] = dGP[curSeq][2]
        ct += 1
    if len(lSeqs) > 0:
        dBThr = (len(lSeqs) - 1)/len(lSeqs)*hmean(lFits)
    for curSeq in lSeqs:
        if dGP[curSeq][1] <= dBThr:
            nUnfAll += 1
    if len(lFits) > 1:
        sdAFit = numpy.std(numpy.array(lFits), dtype = numpy.float64, ddof = 1)
    if numAll > 0:
        mnWtAFit, sdWtAFit = weightedMeanAndSD(numpy.array(lFits),
                                               numpy.array(lProps), numAll)
    return (maxF - minF, maxF, dBThr, nUnfAll, sdAFit, sdWtAFit)

def getAvSequenceDiff(dGP):
    avSqDAbsPc, avSqDWtPc = 0.0, 0.0
    for valGP in dGP.values():
        avSqDAbsPc += valGP[-2]/len(dGP)
        avSqDWtPc += valGP[-1]*valGP[2]
    return avSqDAbsPc, avSqDWtPc

def calcPopFit(dGP):
    popF = 0.0
    for curSeq in dGP:
        popF += dGP[curSeq][2]*dGP[curSeq][3]
    return popF

def getSpecInfo(dI, dGP, maxAllID, curTS):
    lenShare = len(dI['dSpecAna']['shareOfGP'])
    lenCopies = len(dI['dSpecAna']['copiesOfAllThr'])
    lenAge = len(dI['dSpecAna']['ageThrNumGen'])
    lShareGPN, lCpAllN, lAgeThrN = [0]*lenShare, [0]*lenCopies, [0]*lenAge
    lLAllSqs = [[] for i in range(lenShare + lenCopies + lenAge)]
    testShareIdx = [True]*lenShare
    shareCalcOpen, copiesCalcOpen = True, True
    curGPShare = 0.0
    lAll = sortAlleles(dGP, list(dGP), dI['dAllSortM']['GPcur'])
    for keyGP in lAll:
        valGP = dGP[keyGP]
        if shareCalcOpen:
            curGPShare += valGP[2]
            for idx, val in enumerate(dI['dSpecAna']['shareOfGP']):
                if testShareIdx[idx]:
                    lShareGPN[idx] += 1
                    lLAllSqs[idx].append(keyGP)
                    if round(curGPShare,
                             ROUND_PREC - 4) >= round(1.0 - val, ROUND_PREC):
                        testShareIdx[idx] = False   # this index is done
                        if idx == 0:    # now we're done (list ascending)
                            shareCalcOpen = False
        if copiesCalcOpen:
            curNumAllCopies = round(valGP[2]*2*dI['popSize'])
            for idx, val in enumerate(dI['dSpecAna']['copiesOfAllThr']):
                if curNumAllCopies >= val:
                    lCpAllN[idx] += 1
                    lLAllSqs[lenShare + idx].append(keyGP)
                else:
                    if idx == 0:    # now we're done (list ascending)
                        copiesCalcOpen = False
                    break
        for idx, val in enumerate(dI['dSpecAna']['ageThrNumGen']):
            if len(valGP) <= 6: # dGP does not contain "time step emerged" info
                allAgeThr = round(val*2*dI['popSize']*dI['pPtMutAll'])
                if valGP[0] <= maxAllID - allAgeThr:
                    lAgeThrN[idx] += 1
                    lLAllSqs[lenShare + lenCopies + idx].append(keyGP)
                else:
                    break
            else:               # dGP add. contains "time step emerged" info
                if valGP[4] <= curTS - val:
                    lAgeThrN[idx] += 1
                    lLAllSqs[lenShare + lenCopies + idx].append(keyGP)
                else:
                    break
    return (lShareGPN + lCpAllN + lAgeThrN, lLAllSqs)

def complementCheckDictAgeResult(dAR, nAgeCl, nAllFile):
    nAllCalc = [0]*3
    shareAllCalc = [0.0]*3
    for oneCl in range(nAgeCl):
        if oneCl + 1 not in dAR:
            dAR[oneCl + 1] = [0, 0.0]*3
        else:
            for tp in range(3):
                nAllCalc[tp] += dAR[oneCl + 1][tp*2]
                shareAllCalc[tp] += dAR[oneCl + 1][tp*2 + 1]
    for tp in range(3):
        if nAllCalc[tp] != nAllFile:
            print('Error: number of alleles of type', tp + 1, 'unequal (',
                  nAllCalc[tp], ',', nAllFile, ')')
        if round(shareAllCalc[tp], ROUND_PREC - 4) != 1.0:
            print('Error: total allele share of type', tp + 1, 'is',
                  round(shareAllCalc[tp], ROUND_PREC - 4), 'and not 1.0')

def SUB_doClassAnalysis(theDict, curUnit, curVal, curGPVal, numCl):
    for allCl in range(numCl):
        if curGPVal <= curUnit*(allCl + 1):
            addListValToDict(theDict, allCl + 1, curVal)
            break
    if curGPVal > curUnit*numCl:             # add these to last class
        addListValToDict(theDict, numCl, curVal)

def doClassAnalysis(dI, dGP, dAgeR, maxAllID):
    numAgeCl = dI['dOutput']['numAgeCLasses']
    for valGP in dGP.values():
        # age class analysis
        valToAdd = [1, valGP[2], 0, 0.0, 0, 0.0]
        if len(valGP) <= 6:     # dGP does not contain "time step emerged" info
            IDUnit = (maxAllID + 1)/numAgeCl
            SUB_doClassAnalysis(dAgeR, IDUnit, valToAdd, valGP[0], numAgeCl)
        else:                   # dGP add. contains "time step emerged" info
            timeUnit = dI['numIter']/numAgeCl
            SUB_doClassAnalysis(dAgeR, timeUnit, valToAdd, valGP[4], numAgeCl)
        # fitness class analysis
        valToAdd = [0, 0.0, 1, valGP[2], 0, 0.0]
        fitUnit = (dI['maxAFit'] - dI['minAFit'])/numAgeCl
        SUB_doClassAnalysis(dAgeR, fitUnit, valToAdd, valGP[1], numAgeCl)
        # frequency class analysis
        valToAdd = [0, 0.0, 0, 0.0, 1, valGP[2]]
        propUnit = 1.0/numAgeCl
        SUB_doClassAnalysis(dAgeR, propUnit, valToAdd, valGP[2], numAgeCl)
    complementCheckDictAgeResult(dAgeR, numAgeCl, len(dGP))

def analyseGPDetailed(dI, dGP, dAnaRes, dSpcRes, dAgeRes, nFile):
    curTS = getTimeStepFromFileName(nFile, -1, 0, 8)
    # calculate raw data analysis file data
    dSeqDiffGP = getSeqDiffMatrix(dGP)
    dHetAdv = calcHetAdvFromRawData(dI, dGP, dSeqDiffGP)
    (rangeFit, maxFit, dBThresh,
     nUnfitAll, sdAllFit, sdWtAllFit) = calcAlleleFitVals(dGP, list(dGP),
                                                          2*dI['popSize'])
    avHetAdvAbs = dHetAdv['avFitHet'] - dHetAdv['avFitHom']
    avPercPDHet, avHetAdvPerc, avOverdHet, avOverdHetPerc = 0.0, 0.0, 0.0, 0.0
    if dHetAdv['avFitHom'] > 0.0:
        avHetAdvPerc = avHetAdvAbs/dHetAdv['avFitHom']*100.0
        avOverdHet = dHetAdv['avOverdom']/dHetAdv['avFitHom']*100.0
        if dHetAdv['propOfHet'] > 0.0:
            avPercPDHet = round(dHetAdv['avNumDiff']/dHetAdv['propOfHet']/
                                dI['numAACpA']*100.0, 2)
            avOverdHet = dHetAdv['avOverdom']/dHetAdv['propOfHet']
            avOverdHetPerc = avOverdHet/dHetAdv['avFitHom']*100.0
    avSqDffAbsPc, avSqDffPrWtPc = getAvSequenceDiff(dGP)
    lDataA = [len(dGP), round(dHetAdv['avNumDiff'], 2), avPercPDHet,
              dHetAdv['propOfHom']*100.0, dHetAdv['propOfHet']*100.0,
              dHetAdv['avFitHom'], dHetAdv['avFitHet'], avHetAdvAbs,
              avHetAdvPerc, dHetAdv['avOverdom'], avOverdHet, avOverdHetPerc,
              rangeFit, rangeFit/maxFit*100.0, dBThresh, nUnfitAll,
              nUnfitAll/len(dGP)*100.0, sdAllFit, sdWtAllFit, avSqDffAbsPc,
              avSqDffPrWtPc, calcPopFit(dGP)]
    dAnaRes[curTS] = lDataA
    # calculate special analysis file data
    maxID = max(getListOfDictVals(dGP, 0))
    (lAllNums, lLAllSqs) = getSpecInfo(dI, dGP, maxID, curTS)
    lDataS = [len(dGP)] + lAllNums
    for curLAllSeqs in [list(dGP)] + lLAllSqs:
        rangeFitPerc = 0.0
        (rangeFit, maxFit, dBThresh,
         nUnfitAll, sdAllFit, sdWtAllFit) = calcAlleleFitVals(dGP, curLAllSeqs,
                                                              2*dI['popSize'])
        if maxFit > 0.0:
            rangeFitPerc = rangeFit/maxFit*100.0
        lDataS += [rangeFitPerc, sdAllFit, sdWtAllFit, dBThresh, nUnfitAll]
    dSpcRes[curTS] = lDataS
    # calculate age class analysis file data
    if curTS == dI['numIter']:
        doClassAnalysis(dI, dGP, dAgeRes, maxID)

def writeResultsDictToFile(dRes, cRp, cD, nDRes, nF, nFEnd, fKey = True):
    fRes = open(cD + '/' + nDRes + '/' + nF + '_Rep' + str(cRp) + nFEnd, 'a')
    for oneKey in sorted(dRes):
        if fKey:
            fRes.write(str(oneKey))
        for idx, val in enumerate(dRes[oneKey]):
            if not fKey and idx == 0:
                fRes.write(str(val))
            else:
                fRes.write(' ' + str(val))
        fRes.write('\n')
    fRes.close()

def getDictAlleleStatus(dGPP, dGPC):
    dLost, dPrev, dCur, dNew = {}, {}, {}, {}
    for oneSeq in dGPP:
        if oneSeq in dGPC:
            dPrev[oneSeq] = dGPP[oneSeq]
            dCur[oneSeq] = dGPC[oneSeq]
        else:
            dLost[oneSeq] = dGPP[oneSeq]
    for oneSeq in dGPC:
        if oneSeq not in dGPP:
            dNew[oneSeq] = dGPC[oneSeq]
    return dLost, dPrev, dCur, dNew

def SUB_writeAlleleCompLines(dI, theDict, theTSt, theStat, cRp, cDir):
    fRes = open(cDir + '/' + dI['nD_Results'] + '/' + dI['nF_AllRes'] +
                '_Rep' + str(cRp) + dI['nF_End'], 'a')
    for oneSeq in theDict:
        fRes.write(str(theTSt) + ' ' + theStat + ' ')
        for oneEl in theDict[oneSeq][:-1]:
            fRes.write(str(oneEl) + ' ')
        fRes.write(str(theDict[oneSeq][-1]) + '\n')
    fRes.close()

def writeToAlleleCompFile(dI, dGPPrev, dGPCur, nFResP, nFResC, cRp, cDir):
    prevTSt = getTimeStepFromFileName(nFResP, -1, 0, 8)
    curTSt = getTimeStepFromFileName(nFResC, -1, 0, 8)
    dLost, dPrev, dCur, dNew = getDictAlleleStatus(dGPPrev, dGPCur)
    SUB_writeAlleleCompLines(dI, dLost, prevTSt, 'L', cRp, cDir)
    SUB_writeAlleleCompLines(dI, dPrev, prevTSt, 'P', cRp, cDir)
    SUB_writeAlleleCompLines(dI, dCur, curTSt, 'C', cRp, cDir)
    SUB_writeAlleleCompLines(dI, dNew, curTSt, 'N', cRp, cDir)

def getSimuStatus(dictI, curDir):
    numStRp = 0
    # first count the number of rate result files
    pFRateRes = curDir + '/' + dictI['nD_Results'] + '/' + dictI['nF_RateRes']
    for curRep in range(1, dictI['numRep'] + 1):
        if os.path.isfile(pFRateRes + '_Rep' + str(curRep) + dictI['nF_End']):
            numStRp = max(numStRp, curRep)
    numFinRp = numStRp
    # then check whether the last rate result file is finished
    pFRp1 = pFRateRes + '_Rep' + str(1) + dictI['nF_End']
    pFRpMx = pFRateRes + '_Rep' + str(numStRp) + dictI['nF_End']
    if os.path.isfile(pFRp1) and os.path.isfile(pFRpMx):
        fRp1 = open(pFRp1, 'r')
        lFRp1 = fRp1.readlines()
        fRp1.close()
        fRpMx = open(pFRpMx, 'r')
        lFRpMx = fRpMx.readlines()
        fRpMx.close()
        if len(lFRpMx) < len(lFRp1):    # last rate result file is not finished
            numFinRp -= 1
    return numFinRp, numStRp

def analyseRawDataFiles(dictI, curRp, curDir):
    dictAnaRes, dictSpcRes, dictAgeRes = {}, {}, {}
    lNFRes = getListOfResultFileNames(curDir, curRp, dictI['nD_Results'],
                                      dictI['nF_GPRes'])
    lNFRes = sortListOfFileNames(lNFRes)
    writeHeaderOfRawDataAnalysisFile(dictI, curRp, curDir)
    writeHeaderOfSpecialAnalysisFile(dictI, curRp, curDir)
    writeHeaderOfAgeClassAnalysisFile(dictI, curRp, curDir)
    writeHeaderOfAlleleComparisonFile(dictI, curRp, curDir)
    for nFRes in lNFRes:
        if len(nFRes) >= len(dictI['nF_GPRes']):
            if nFRes[:len(dictI['nF_GPRes'])] == dictI['nF_GPRes']:
                dictGP = parseGPResFile(dictI, nFRes, curDir)
                analyseGPDetailed(dictI, dictGP, dictAnaRes, dictSpcRes,
                                  dictAgeRes, nFRes)
                cIdx = lNFRes.index(nFRes)
                if cIdx > 0:
                    nResFileP = lNFRes[cIdx - 1]
                    dictGPPrev = parseGPResFile(dictI, nResFileP, curDir)
                    writeToAlleleCompFile(dictI, dictGPPrev, dictGP, nResFileP,
                                          nFRes, curRp, curDir)
    writeResultsDictToFile(dictAnaRes, curRp, curDir, dictI['nD_Results'],
                           dictI['nF_AnaRes'], dictI['nF_End'])
    writeResultsDictToFile(dictSpcRes, curRp, curDir, dictI['nD_Results'],
                           dictI['nF_SpecRes'], dictI['nF_End'])
    writeResultsDictToFile(dictAgeRes, curRp, curDir, dictI['nD_Results'],
                           dictI['nF_AgeRes'], dictI['nF_End'])
    print('Printed raw data analysis files for repetition', curRp)

def addToValList(theDict, lSplits, splStr1, splStr2, str1, str2, nUsed, nF):
    if splStr1 == str1 and splStr2 == str2:
        if int(lSplits[1][3:]) <= nUsed:
            theDict[str1].append(nF)

def getRepResFiles(dI, pFR, nFs):
    dFinFs = {dI['nF_RateRes']: [], dI['nF_AnaRes']: [], dI['nF_SpecRes']: [],
              dI['nF_AgeRes']: [], dI['nF_HetAdv']: []}
    lNFiles = os.listdir(pFR)
    for nFile in lNFiles:
        nFileNoExt = nFile.split('.')[0]
        lULSplit = nFileNoExt.split('_')
        if len(lULSplit) > 1:           # only these files are candidates
            nF12 = [lULSplit[0], lULSplit[1]]
            addToValList(dFinFs, lULSplit, nF12[0], nF12[1][0],
                         dI['nF_RateRes'], 'R', nFs, nFile)
            addToValList(dFinFs, lULSplit, nF12[0], nF12[1][0],
                         dI['nF_AnaRes'], 'R', nFs, nFile)
            addToValList(dFinFs, lULSplit, nF12[0], nF12[1][0],
                         dI['nF_SpecRes'], 'R', nFs, nFile)
            addToValList(dFinFs, lULSplit, nF12[0], nF12[1][0],
                         dI['nF_AgeRes'], 'R', nFs, nFile)
            addToValList(dFinFs, lULSplit, nF12[0], nF12[1][0],
                         dI['nF_HetAdv'], 'R', nFs, nFile)
    return dFinFs

def calcMeans(pFR, numFs, stIdx, fTp, dFinFs):
    tbHdr = ''
    lfDtM = []
    for nF in dFinFs[fTp]:
        fCur = open(pFR + '/' + nF, 'r')
        lFCur = fCur.readlines()
        tbHdr = lFCur[0]
        numLines = len(lFCur[1:])
        if dFinFs[fTp].index(nF) == 0:  # first file
            lfDtM = [[]]*numLines
        for lineIdx in range(numLines):
            lNewVs = lFCur[1:][lineIdx].split(' ')
            lNewVs = [int(s) for s in lNewVs[:stIdx]] + \
                [float(s) for s in lNewVs[stIdx:]]
            if len(lfDtM[lineIdx]) > 0:        # not the first file
                for elIdx in range(stIdx, len(lNewVs)):
                    lfDtM[lineIdx][elIdx] += lNewVs[elIdx]/numFs
            else:                               # first file
                lfDtM[lineIdx] = lNewVs[:stIdx]
                lfDtM[lineIdx] += [f/numFs for f in lNewVs[stIdx:]]
        fCur.close()
    return tbHdr, lfDtM

def calcSDs(pFR, numFs, stIdx, fTp, dFinFs, lfDtM):
    lfDtS = []
    for nF in dFinFs[fTp]:
        fCur = open(pFR + '/' + nF, 'r')
        lFCur = fCur.readlines()
        numLines = len(lFCur[1:])
        if dFinFs[fTp].index(nF) == 0:  # first file
            lfDtS = [[]]*numLines
        for lineIdx in range(numLines):
            lNewVs = lFCur[1:][lineIdx].split(' ')
            lNewVs = [int(s) for s in lNewVs[:stIdx]] + \
                [float(s) for s in lNewVs[stIdx:]]
            if len(lfDtS[lineIdx]) > 0:        # not the first file
                for elIdx in range(stIdx, len(lNewVs)):
                    smd = (lNewVs[elIdx] - lfDtM[lineIdx][elIdx])* \
                        (lNewVs[elIdx] -lfDtM[lineIdx][elIdx])/(numFs - 1)
                    lfDtS[lineIdx][elIdx] += smd
            else:
                lfDtS[lineIdx] = lNewVs[:stIdx]
                lfDtS[lineIdx] += [(lNewVs[i] - lfDtM[lineIdx][i])*
                                   (lNewVs[i] - lfDtM[lineIdx][i])/
                                   (numFs - 1)
                                   for i in range(stIdx, len(lNewVs))]
        fCur.close()
    # calculate the standard deviations from the variances
    for lineIdx in range(len(lfDtS)):
        for elIdx in range(stIdx, len(lfDtS[lineIdx])):
            lfDtS[lineIdx][elIdx] = sqrt(lfDtS[lineIdx][elIdx])
    return lfDtS

def mergeRepRes(dictI, curDir, numFs):
    pFRes = curDir + '/' + dictI['nD_Results']
    dictFinFs = getRepResFiles(dictI, pFRes, numFs)
    for fType in dictFinFs:
        startIdx = 1
        if fType == dictI['nF_HetAdv']:
            startIdx = 0
        tabHdr, lfDatM = calcMeans(pFRes, numFs, startIdx, fType, dictFinFs)
        lfDatS = calcSDs(pFRes, numFs, startIdx, fType, dictFinFs, lfDatM)
        pFAvs = pFRes + '/' + fType + '_Av' + dictI['nF_End']
        pFSDs = pFRes + '/' + fType + '_SD' + dictI['nF_End']
        writeListToFile(pFAvs, tabHdr, lfDatM)
        writeListToFile(pFSDs, tabHdr, lfDatS)
        print('Printed average result file', fType + '_Av' + dictI['nF_End'])
        print('Printed st. dev. result file', fType + '_SD' + dictI['nF_End'])
    print('Finished calculating and printing means and standard deviations.')

# MAIN PROGRAM #
dictInp = genDictInp(iniAllFDist, iniAllProp, newAllFDist, limFitness,
                     dominType, deltaIn, gaussMuIn, gaussSigIn, minFitn,
                     maxFitn, minGTFitn, maxGTFitn, dictDefFit, lostThr,
                     disregThr, numAllThr, numRepetitions, numIterations,
                     numARIter, popSize, numAllI, numAmAcPerAll, hetAdvM,
                     lambdaDeB, fixedHA, modDisp, modRate, modFile, dAllSortM,
                     dFlags, dMutInfo, dSpecAna, dClOutc, dOutput, ND_RESULTS,
                     NF_RATERES, NF_ANARES, NF_SPECRES, NF_AGERES, NF_ALLRES,
                     NF_PHYL, NF_GPRES, NF_SGPRES, NF_TGPRES, NF_HA, NF_DIVRES,
                     NF_REPRES, NF_END)
numFinRep, numStRep = getSimuStatus(dictInp, os.getcwd())
for curRep in range(1, numStRep + 1):
    analyseRawDataFiles(dictInp, curRep, os.getcwd())
if numFinRep > 1:
    mergeRepRes(dictInp, os.getcwd(), numFinRep)
