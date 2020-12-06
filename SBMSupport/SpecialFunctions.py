# -*- coding: utf-8 -*-
################################################################################
# SpecialFunctions.py #
################################################################################

import random, time
from operator import itemgetter

import InputVariables
from Constants import *
from GenericFunctions import (drawFromPList, randPartitionRange,
                              calcCondBinProbs, SUB_G2Fit, getListOfDictVals,
                              getDictNumEvents, drawFromBinom, drawFromPoisson,
                              drawFromMultinom, safeRemoveEl, compare2Seq,
                              capValue, getStartKeyVal, getKeyOfExtVal)

# ### FUNCTIONS ################################################################
def genDictInp(iniAllFDist, iniAllProp, newAllFDist, limFit, domType, deltaIn,
               gaussMuIn, gaussSigIn, minAFit, maxAFit, minGTFit, maxGTFit,
               dictDefFit, lostThr, disregThr, numAllThr, numRep, numIter,
               numARIter, popSize, numAllIni, numAACpA, hetAdvM, lambdaDeB,
               fixedHA, modDisp, modRate, modFile, dAllSortM, dFlags, dMutInfo,
               dSpecAna, dClOutc, dOutput, nD_Results, nF_RateRes, nF_AnaRes,
               nF_SpecRes, nF_AgeRes, nF_AllRes, nF_Phyl, nF_GPRes, nF_SGPRes,
               nF_TGPRes, nF_HetAdv, nF_DivRes, nF_RepRes, nF_End):
    theDict = {'iniAllFDist': iniAllFDist,
               'iniAllProp': iniAllProp,
               'newAllFDist': newAllFDist,
               'limFit': limFit,
               'domType': domType,
               'deltaIn': deltaIn,
               'gaussMuIn': gaussMuIn,
               'gaussSigIn': gaussSigIn,
               'minAFit': minAFit,
               'maxAFit': maxAFit,
               'minGTFit': minGTFit,
               'maxGTFit': maxGTFit,
               'dictDefFit': dictDefFit,
               'lostThr': lostThr,
               'disregThr': min(disregThr, lostThr),
               'numAllThr': numAllThr,
               'numRep': numRep,
               'numIter': numIter,
               'numARIter': numARIter,
               'popSize': popSize,
               'numAllIni': numAllIni,
               'numAACpA': numAACpA,
               'hetAdvM': hetAdvM,
               'lambdaDeB': lambdaDeB,
               'fixedHA': fixedHA,
               'modDisp': modDisp,
               'modRate': modRate,
               'modFile': modFile,
               'dAllSortM': dAllSortM,
               'dFlags': dFlags,
               'dMutInfo': dMutInfo,
               'dSpecAna': dSpecAna,
               'dClOutc': dClOutc,
               'dOutput': dOutput,
               'pPtMutSglAA': dMutInfo['pPtMutSglAA'],
               'pPtMutAll': 1.0-(1.0-dMutInfo['pPtMutSglAA'])**numAACpA,
               'pMicConvMut': dMutInfo['pMicConvMut'],
               'pNewSeqMut': dMutInfo['pNewSeqMut'],
               'lCBinPrs': calcCondBinProbs(numAACpA, dMutInfo['pPtMutSglAA']),
               'nD_Results': nD_Results,
               'nF_RateRes': nF_RateRes,
               'nF_AnaRes': nF_AnaRes,
               'nF_SpecRes': nF_SpecRes,
               'nF_AgeRes': nF_AgeRes,
               'nF_AllRes': nF_AllRes,
               'nF_Phyl': nF_Phyl,
               'nF_GPRes': nF_GPRes,
               'nF_SGPRes': nF_SGPRes,
               'nF_TGPRes': nF_TGPRes,
               'nF_HetAdv': nF_HetAdv,
               'nF_DivRes': nF_DivRes,
               'nF_RepRes': nF_RepRes,
               'nF_End': nF_End}
    return theDict

def genDictResRepRes(dictI):
    dictRs = {'timeStep': 0, 'numSglPtMutTot': 0, 'numAllPtMutTot': 0,
              'numMicConvMutTot': 0, 'numNewSeqMutTot': 0,
              'numAllTot': dictI['numAllIni']}
    dClOutcInc = {}
    for oneCl in dictI['dClOutc']:
        dClOutcInc[oneCl] = 0
    dictRpRs = {'clOutcInc': dClOutcInc,
                'dAllPers': {},
                'dHugeNAllPers': {}}
    return dictRs, dictRpRs

def resetDictRes(dictI, dictRs):
    for oneKey in dictRs:
        if oneKey == 'numAllTot':
            dictRs[oneKey] = dictI['numAllIni']
        else:
            dictRs[oneKey] = 0

def assignAllFit(dI, allID):
    allFit = 0.0
    if dI['iniAllFDist'] == A_INI_FDIST_UNIDET:
        allFit = round(dI['minAFit'] + (dI['maxAFit'] - dI['minAFit'])/
                       (2*dI['numAllIni'])*(2*allID + 1), ROUND_PREC)
    elif dI['iniAllFDist'] == A_INI_FDIST_UNIRAND:
        allFit = random.uniform(dI['minAFit'], dI['maxAFit'])
    elif dI['iniAllFDist'] == A_INI_FDIST_GAUSS:
        allFit = random.gauss(dI['gaussMuIn'], dI['gaussSigIn'])
        allFit = capValue(allFit, dI['minAFit'], dI['maxAFit'],
                          dI['limFit'] in [L_FIT_1, L_FIT_MAX])
    elif dI['iniAllFDist'] == A_INI_FDIST_DEFINED:
        allFit = dI['dictDefFit'][allID%len(dI['dictDefFit'])]
    else: pass
    return allFit

def newAllFit(dI, fOrigA):
    allFit = 0.0
    if dI['newAllFDist'] == A_NEW_FDIST_UNIRAND:
        if dI['dFlags']['mutFitDepend']:
            allFit = random.uniform(fOrigA - dI['dMutInfo']['newFitVarUniR'],
                                    fOrigA + dI['dMutInfo']['newFitVarUniR'])
        else:
            allFit = random.uniform(dI['minAFit'], dI['maxAFit'])
    elif dI['newAllFDist'] == A_NEW_FDIST_GAUSS:
        if dI['dFlags']['mutFitDepend']:
            allFit = random.gauss((fOrigA + dI['gaussMuIn'])/2.0,
                                  dI['dMutInfo']['newFitVarGauss'])
        else:
            allFit = random.gauss(dI['gaussMuIn'], dI['gaussSigIn'])
    else: pass
    allFit = capValue(allFit, dI['minAFit'], dI['maxAFit'])
    return allFit

def genNewSequence(dI, nAAc):
    AAcS = ''
    for j in range(dI['numAACpA']):
        AAcS += DICT_AMINO_ACIDS[random.randint(1, nAAc)][2]
    return AAcS

def genDictGenepools(dictI):
    dictGP = {}
    dictGPTot = {}
    numAAc = len(DICT_AMINO_ACIDS)
    for i in range(dictI['numAllIni']):
        AAcSeq = genNewSequence(dictI, numAAc)
        dictGP[AAcSeq] = [i, assignAllFit(dictI, i), 0.0, 0.0, 0, -1]
        # ID, fitness, proportion, marginal fitness, emerged, mut. type
        dictGPTot[AAcSeq] = dictGP[AAcSeq][1]
        # fitness
    return dictGP, dictGPTot

def calcAllProp(dictI, dictGP):
    numAll = len(dictGP)
    lAllS = list(dictGP)
    if dictI['iniAllProp'] == A_PROP_UNIRAND:
        lProps = randPartitionRange(numAll)
        for allS in lAllS:
            dictGP[allS][2] = lProps[lAllS.index(allS)]
    elif dictI['iniAllProp'] == A_PROP_UNIDET:
        for allS in lAllS:
            dictGP[allS][2] = round(1.0/numAll, ROUND_PREC)
    else: pass

def convGene2Fit(dI, dGP, seqAll1, seqAll2):
    fAllC1 = dGP[seqAll1][1]    # fitness of allele on chromosome 1
    fAllC2 = dGP[seqAll2][1]    # fitness of allele on chromosome 2
    return SUB_G2Fit(seqAll1, seqAll2, fAllC1, fAllC2, dI['minGTFit'],
                     dI['maxGTFit'], dI['domType'], dI['hetAdvM'],
                     dI['deltaIn'], dI['lambdaDeB'], dI['fixedHA'],
                     dI['limFit'])

def genDictGenotypeFit(dictI, dictGP):
    dictGTF = {}
    lAllS = list(dictGP)
    for allS in lAllS:
        aIDC1 = dictGP[allS][0]
        for aIdxC2 in range(lAllS.index(allS) + 1):
            aIDC2 = dictGP[lAllS[aIdxC2]][0]
            keyGTF = (min(aIDC1, aIDC2), max(aIDC1, aIDC2)) # ID[0] <= ID[1]
            dictGTF[keyGTF] = convGene2Fit(dictI, dictGP, allS, lAllS[aIdxC2])
    return dictGTF

def calcMargFit(allS1, dGP, dGTF, roundPrec = 16):
    mF = 0.0
    valS1 = dGP[allS1]
    for allS2 in dGP:
        valS2 = dGP[allS2]
        # marginal fitness = sum_over_GTs_with_all_i(prop_all_j*fit_GT_ij)
        mF += valS2[2]*dGTF[(min(valS1[0], valS2[0]),
                             max(valS1[0], valS2[0]))]
    return round(mF, roundPrec)

def updateMargPopFit(dGP, dGTF, roundPrec = 16):
    popF = 0.0
    for allS in dGP:
        margF = calcMargFit(allS, dGP, dGTF, roundPrec)
        dGP[allS][3] = margF
        valS = dGP[allS]
        popF += valS[2]*margF
    return popF

def updateDictRes(dI, dGP, dR, popF = -1.0):
    dR['propZyg'] = {'propHomoZyg': 0.0, 'propHeteroZyg': 1.0}
    dR['lAllIDsPers'] = getListOfDictVals(dGP, 0, True)
    dR['lAllSeqsPers'] = list(dGP)
    for allS in dGP:
        # check if some alleles in dGP are below the lost threshold
        valS = dGP[allS]
        if valS[2] < dI['lostThr']:
            dR['lAllIDsPers'].remove(valS[0])
            dR['lAllSeqsPers'].remove(allS)
        dR['propZyg']['propHomoZyg'] += valS[2]*valS[2]
    dR['nAllPers'] = len(dR['lAllIDsPers'])
    dR['propZyg']['propHeteroZyg'] = 1.0 - dR['propZyg']['propHomoZyg']
    if popF >= 0:   # only set if a new population fitness value is given 
        dR['popFitness'] = round(popF, ROUND_PREC)

def updateDicts(dictI, dictGnP, dictGTF, dictR):
    # initiate result dictionary
    popFit = updateMargPopFit(dictGnP, dictGTF, ROUND_PREC)
    updateDictRes(dictI, dictGnP, dictR, popFit)

def checkIfAllLost(dI, dGP, cSeq, tS):
    if dGP[cSeq][2] < dI['disregThr']:    # this allele is lost
        safeRemoveEl(dGP, cSeq, 2, dGP[cSeq][2])

def adjProportions(dI, dGP, dGPT, dGTF, dR, dMutCt, curS, newS, mutT = 1):
    newProp = min(1/(2*dI['popSize']), dGP[curS][2])    # 'atomic': 1/num_all
    dGP[curS] = dGP[curS][:2] + [dGP[curS][2] - newProp, 0.0] + dGP[curS][4:]
    if newS in dGP:     # new sequence in current gene pool (back-mut.)
        dR['numAllTot'] -= 1
        dGP[newS][2] += newProp
    else:               # this sequence in not in the current gene pool
        newAllID = dR['numAllTot'] - 1
        if 'numAllPtMutTS' in dMutCt:
            newAllID += dMutCt['numAllPtMutTS']
        elif 'numMicConvMutTS' in dMutCt:
            newAllID += dMutCt['numMicConvMutTS']
        elif 'numNewSeqMutTS' in dMutCt:
            newAllID += dMutCt['numNewSeqMutTS']
        else: pass
        if newS not in dGPT:    # new allele has never existed before
            newFit = newAllFit(dI, dGPT[curS])
            dGPT[newS] = newFit
        else:
            newFit = dGPT[newS]
        dGP[newS] = [newAllID, newFit, newProp, 0.0, dR['timeStep'], mutT]
        dMutCt['lNewSeqs'].append(newS)
    # delete current allele from gene pool if proportion is too low
    checkIfAllLost(dI, dGP, curS, dR['timeStep'])

def updateDictGenotypeFit(dI, dGP, dGTF, dMutCt):
    for cSeq in dGP:
        cAllID = dGP[cSeq][0]
        for nSeq in dMutCt['lNewSeqs']:
            nAllID = dGP[nSeq][0]
            dGTF[(cAllID, nAllID)] = convGene2Fit(dI, dGP, cSeq, nSeq)

def doUpdates(dI, dGP, dGTF, dR, dMutCt, dNMut):
    popF = dR['popFitness']
    if 'numAllPtMutTS' in dMutCt:
        if dMutCt['numAllPtMutTS'] >= 1:
            dR['numSglPtMutTot'] += dMutCt['numSglPtMutTS']
            dR['numAllPtMutTot'] += dMutCt['numAllPtMutTS']
            dR['numAllTot'] += dMutCt['numAllPtMutTS']
            updateDictGenotypeFit(dI, dGP, dGTF, dMutCt)
    elif 'numMicConvMutTS' in dMutCt:
        if dMutCt['numMicConvMutTS'] >= 1:
            dR['numMicConvMutTot'] += dMutCt['numMicConvMutTS']
            dR['numAllTot'] += dMutCt['numMicConvMutTS']
            updateDictGenotypeFit(dI, dGP, dGTF, dMutCt)
    elif 'numNewSeqMutTS' in dMutCt:
        if dMutCt['numNewSeqMutTS'] >= 1:
            dR['numNewSeqMutTot'] += dMutCt['numNewSeqMutTS']
            dR['numAllTot'] += dMutCt['numNewSeqMutTS']
            updateDictGenotypeFit(dI, dGP, dGTF, dMutCt)
    else: pass
    dR['popFitness'] = popF

def newSeqOfPtMut(cSeq, lDrPos):
    for mutPos in lDrPos:
        # draw a new amino acid on that position
        newAAcIdx = random.randint(1, len(DICT_AMINO_ACIDS))
        newAAc = DICT_AMINO_ACIDS[newAAcIdx][2]
        while newAAc == cSeq[mutPos]:
            newAAcIdx = random.randint(1, len(DICT_AMINO_ACIDS))
            newAAc = DICT_AMINO_ACIDS[newAAcIdx][2]
        nSeq = (cSeq[:mutPos] + newAAc + cSeq[mutPos + 1:])
        cSeq = nSeq
    return nSeq

def ptMut(dI, dGP, dGPT, dGTF, dNMut, dR):
    dMutCt = {'numSglPtMutTS': 0, 'numAllPtMutTS': 0, 'lNewSeqs': []}
    for curSeq in dNMut:
        for mutIdx in range(dNMut[curSeq]):
            numPtMutCurAll = drawFromPList(dI['lCBinPrs'])
            dMutCt['numSglPtMutTS'] += numPtMutCurAll
            dMutCt['numAllPtMutTS'] += 1
            lPtMutPos = random.sample(range(dI['numAACpA']), numPtMutCurAll)
            newSeq = newSeqOfPtMut(curSeq, lPtMutPos)
            # update the allele proportion in the genepool dictionary
            adjProportions(dI, dGP, dGPT, dGTF, dR, dMutCt, curSeq, newSeq,
                           T_PT_MUT)
            if curSeq not in dGP:   # it could have been removed just before
                break
    # if >= 1 mutation, update marginal / population fitness, dictGenotypeFit
    doUpdates(dI, dGP, dGTF, dR, dMutCt, dNMut)

def chooseDonor(dGP):
    drIdx = drawFromPList(getListOfDictVals(dGP, 2))
    chosenSeq = ''
    seqIdx = 0
    for oneSeq in dGP:
        if seqIdx == drIdx:
            chosenSeq = oneSeq
            break
        seqIdx += 1
        if seqIdx >= len(dGP):
            print('Error: Index exceeds dictionary length.')
            return ERR_CODE
    return chosenSeq

def micConvMut(dI, dGP, dGPT, dGTF, dNMut, dR):
    dMutCt = {'numMicConvMutTS': 0, 'lNewSeqs': []}
    for curSeq in dNMut:
        for mutIdx in range(dNMut[curSeq]):
            dMutCt['numMicConvMutTS'] += 1
            # draw a random allele as donor
            donorAll = chooseDonor(dGP)
            # draw the donor segment length from a Gaussian distribution
            lenSeg = random.gauss(dI['dMutInfo']['muLMicConvSeq'],
                                  dI['dMutInfo']['sigLMicConvSeq'])
            lenSeg = round(capValue(lenSeg, 1, dI['numAACpA']))
            # draw start pos. of donor and acceptor allele, create new allele
            stPDonor = random.randint(0, dI['numAACpA'] - lenSeg)
            donorSeg = ''
            if random.randint(0, 1) == 0:
                donorSeg = donorAll[stPDonor:stPDonor + lenSeg]
            else:   # invert the sequence
                donorSeg = donorAll[stPDonor + lenSeg - 1:stPDonor - 1:-1]
            stPAcc = random.randint(0, dI['numAACpA'] - lenSeg)
            newSeq = curSeq[:stPAcc] + donorSeg + curSeq[stPAcc + lenSeg:]
            # update the allele proportion in the genepool dictionary
            adjProportions(dI, dGP, dGPT, dGTF, dR, dMutCt, curSeq, newSeq,
                           T_MIC_CONV_MUT)
            if curSeq not in dGP:   # it could have been removed just before
                break
    # if >= 1 mutation, update marginal / population fitness, dictGenotypeFit
    doUpdates(dI, dGP, dGTF, dR, dMutCt, dNMut)

def newSeqMut(dI, dGP, dGPT, dGTF, dNMut, dR):
    dMutCt = {'numNewSeqMutTS': 0, 'lNewSeqs': []}
    for curSeq in dNMut:
        for mutIdx in range(dNMut[curSeq]):
            dMutCt['numNewSeqMutTS'] += 1
            newSeq = genNewSequence(dI, len(DICT_AMINO_ACIDS))
            # update the allele proportion in the genepool dictionary
            adjProportions(dI, dGP, dGPT, dGTF, dR, dMutCt, curSeq, newSeq,
                           T_NEW_SEQ_MUT)
            if curSeq not in dGP:   # it could have been removed just before
                break
    # if >= 1 mutation, update marginal / population fitness, dictGenotypeFit
    doUpdates(dI, dGP, dGTF, dR, dMutCt, dNMut)

def doMutations(dictI, dictGP, dictGPTot, dictGTF, dictR):
    if dictI['pPtMutSglAA'] > 0:
        # do point mutations
        dictNMut = getDictNumEvents(dictGP, 2, 2*dictI['popSize'],
                                    dictI['pPtMutAll'])
        ptMut(dictI, dictGP, dictGPTot, dictGTF, dictNMut, dictR)
    if dictI['pMicConvMut'] > 0:
        # do micro-conversion mutations
        dictNMut = getDictNumEvents(dictGP, 2, 2*dictI['popSize'],
                                    dictI['pMicConvMut'])
        micConvMut(dictI, dictGP, dictGPTot, dictGTF, dictNMut, dictR)
    if dictI['pNewSeqMut'] > 0:
        # do mutations with total sequence exchange
        dictNMut = getDictNumEvents(dictGP, 2, 2*dictI['popSize'],
                                    dictI['pNewSeqMut'])
        newSeqMut(dictI, dictGP, dictGPTot, dictGTF, dictNMut, dictR)

def propUpdatePerGenMultin(dictI, dictGP, dictGTF, dictR):
    lAllCur, lNumAll = list(dictGP), [] # some keys might be deleted
    lPDet = [0.0]*len(lAllCur)  # list of new deterministic proportions
    curIdx = 0
    for curSeq in lAllCur:
        # update formula: prop(now) = prop(before)*margFit/popFit
        if dictGP[curSeq][2] > 0:
            newPDet = dictGP[curSeq][2]*dictGP[curSeq][3]/dictR['popFitness']
            lPDet[curIdx] = newPDet
            curIdx += 1
    if dictR['timeStep'] <= dictI['numIter']:
        lNumAll = drawFromMultinom(2*dictI['popSize'], lPDet)
    curIdx = 0
    for curSeq in lAllCur:
        if dictR['timeStep'] <= dictI['numIter']:
            dictGP[curSeq][2] = lNumAll[curIdx]/(2*dictI['popSize'])
        else:
            dictGP[curSeq][2] = lPDet[curIdx]
        curIdx += 1
        checkIfAllLost(dictI, dictGP, curSeq, dictR['timeStep'])

def cleanDictGTF(dictGP, dictGTF):
    lAllIDs = getListOfDictVals(dictGP, 0)
    lIDPairs = list(dictGTF)
    for curIDPair in lIDPairs:
        if curIDPair[0] not in lAllIDs or curIDPair[1] not in lAllIDs:
            del dictGTF[curIDPair]

def oneTimeStepDynamics(dictI, dictGnP, dictGnPTot, dictGTF, dictR, timeS):
    dictR['timeStep'] = timeS
    # update the allele proportions according to the per-generation update rule
    InputVariables.stT_PropUp = time.time()
    propUpdatePerGenMultin(dictI, dictGnP, dictGTF, dictR)
    InputVariables.elT_PropUp += time.time() - InputVariables.stT_PropUp
    if timeS <= dictI['numIter']:
        InputVariables.stT_CleanD = time.time()
        # clean the genotype fitness dictionary after fixed intervals
        if timeS%dictI['dMutInfo']['dGTFCleanTS'] == 0:
            cleanDictGTF(dictGnP, dictGTF)
        InputVariables.elT_CleanD += time.time() - InputVariables.stT_CleanD
        InputVariables.stT_Mut = time.time()
        # possible allele mutations, including marginal fitness update
        doMutations(dictI, dictGnP, dictGnPTot, dictGTF, dictR)
        InputVariables.elT_Mut += time.time() - InputVariables.stT_Mut
    InputVariables.stT_FitUp = time.time()
    # update the marginal fitnesses of all alleles and the population fitness
    popFitNew = updateMargPopFit(dictGnP, dictGTF, ROUND_PREC)
    dictR['popFitness'] = round(popFitNew, ROUND_PREC)
    InputVariables.elT_FitUp += time.time() - InputVariables.stT_FitUp

def getMostCommonAllele(dGP, lAll, lAllUs, aIdx, startK, startV):
    dLIdx = 2
    isMax = True
    outString = 'Error: No allele has a positive frequency...'
    if aIdx == 0:
        # get the lowest value first if sorting after highest - and vice versa
        startK, startV = getStartKeyVal(dGP, lAll, dLIdx, not isMax)
    kExtV = getKeyOfExtVal(dGP, lAll, lAllUs, startK, startV, dLIdx, isMax,
                           outString)
    return startK, startV, kExtV

def getFittestAllele(dGP, lAll, lAllUs, aIdx, startK, startV):
    dLIdx = 1
    isMax = True
    outString = 'Error: No allele has a positive fitness...'
    if aIdx == 0:
        # get the lowest value first if sorting after highest - and vice versa
        startK, startV = getStartKeyVal(dGP, lAll, dLIdx, not isMax)
    kExtV = getKeyOfExtVal(dGP, lAll, lAllUs, startK, startV, dLIdx, isMax,
                           outString)
    return startK, startV, kExtV

def getAlleleLowestID(dGP, lAll, lAllUs, aIdx, startK, startV):
    dLIdx = 0
    isMax = False
    outString = 'Error: No allele has an ID...'
    if aIdx == 0:
        # get the lowest value first if sorting after highest - and vice versa
        startK, startV = getStartKeyVal(dGP, lAll, dLIdx, not isMax)
    kExtV = getKeyOfExtVal(dGP, lAll, lAllUs, startK, startV, dLIdx, isMax,
                           outString)
    return startK, startV, kExtV

def getAlleleHighestID(dGP, lAll, lAllUs, aIdx, startK, startV):
    dLIdx = 0
    isMax = True
    outString = 'Error: No allele has an ID...'
    if aIdx == 0:
        # get the lowest value first if sorting after highest - and vice versa
        startK, startV = getStartKeyVal(dGP, lAll, dLIdx, not isMax)
    kExtV = getKeyOfExtVal(dGP, lAll, lAllUs, startK, startV, dLIdx, isMax,
                           outString)
    return startK, startV, kExtV

def sortAlleles(dGP, lAll, fSortMode):
    if fSortMode == M_NONE:   # do not sort
        return list(dGP)
    lAllSort = []
    lAllUsed = []
    stK, stV, kExtV = '', 0, ''
    for allIdx in range(len(lAll)):
        if fSortMode == M_MOST_COMMON_ALLELE:   # sort - most common allele
            stK, stV, kExtV = getMostCommonAllele(dGP, lAll, lAllUsed, allIdx,
                                                  stK, stV)
        elif fSortMode == M_FITTEST_ALLELE:     # sort - fittest allele
            stK, stV, kExtV = getFittestAllele(dGP, lAll, lAllUsed, allIdx,
                                               stK, stV)
        elif fSortMode == M_LOWEST_ID:          # sort - lowest allele ID
            stK, stV, kExtV = getAlleleLowestID(dGP, lAll, lAllUsed, allIdx,
                                                stK, stV)
        elif fSortMode == M_HIGHEST_ID:         # sort - highest allele ID
            stK, stV, kExtV = getAlleleHighestID(dGP, lAll, lAllUsed, allIdx,
                                                  stK, stV)
        else: pass
        lAllSort.append(kExtV)
    return lAllSort

def getSeqDiffMatrix(dGP):
    dSeqDiff = {}
    lSortedGPLines = sorted(dGP.items(), key = itemgetter(1))
    for curLine in lSortedGPLines:
        seqA1 = curLine[0]
        for idxAll2 in range(lSortedGPLines.index(curLine)):
            seqA2 = lSortedGPLines[idxAll2][0]
            dSeqDiff[(dGP[seqA1][0],
                      dGP[seqA2][0])] = compare2Seq(seqA1, seqA2)
    return dSeqDiff

def getNumberOfSeqDiffs(dI, dGP, dSDGP, lAll, cSeq):
    numSDPerc, wtSDPerc = 0.0, 0.0
    for oSeq in lAll:
        if cSeq != oSeq:
            keyV1 = (dGP[cSeq][0], dGP[oSeq][0])
            keyV2 = (dGP[oSeq][0], dGP[cSeq][0])
            if keyV1 in dSDGP:
                numSD = dSDGP[keyV1]
            elif keyV2 in dSDGP:
                numSD = dSDGP[keyV2]
            else:
                print('ERROR:', keyV1, 'and', keyV2, 'both nonexistent.')
            numSDPerc += numSD/(dI['numAACpA']*(len(lAll) - 1))*100.0
            wtSDPerc += numSD/dI['numAACpA']*dGP[oSeq][2]/(1.0 -
                                                           dGP[cSeq][2])*100.0
    return numSDPerc, wtSDPerc
