# -*- coding: utf-8 -*-
################################################################################
# WriteAndPrint.py #
################################################################################

import os, pprint, pickle
##from operator import itemgetter

from Constants import *
from SpecialFunctions import (updateDictRes, sortAlleles, getSeqDiffMatrix,
                              getNumberOfSeqDiffs)
from GenericFunctions import (addToDictOfLists, getDictKey, calcRunMean,
                              addEntryToDVList)

# ### FUNCTIONS ################################################################
def deleteAllFiles(dictI, cDir):
    resDir = os.path.join(cDir, dictI['nD_Results'])
    if not os.path.isdir(resDir):
        os.mkdir(resDir)
    for nOneF in os.listdir(resDir):
        os.remove(os.path.join(resDir, nOneF))
 
def printCurRepet(cRep, dictI, dictRR):
    print('-------- Repetition', cRep, 'of', dictI['numRep'],
          '--------', end = '')
    if len(dictRR['dHugeNAllPers']) > 0:
        print(' (up to', max(dictRR['dHugeNAllPers']),
              'alleles persisting) ', end = '')
    else:
        print(' (no great run so far) ', end = '')
    print('--------')

def openOutFiles(dictI, dictGP, cRp, cDir):
    # open rate result file
    fRes = open(cDir + '/' + dictI['nD_Results'] + '/' +
                dictI['nF_RateRes'] + '_Rep' + str(cRp) + dictI['nF_End'], 'w')
    fRes.write('TimeStep NumberAllelesCurrent PopulationFitness\n')
    fRes.close()
    
def writeToRateResult(dictI, dictR, cRp, tSt, cDir):
    fRes = open(cDir + '/' + dictI['nD_Results'] + '/' +
                dictI['nF_RateRes'] + '_Rep' + str(cRp) + dictI['nF_End'], 'a')
    fRes.write(str(tSt) + ' ' + str(len(dictR['lAllIDsPers'])) + ' ' +
               str(dictR['popFitness']) + '\n')
    fRes.close()

def printGPResult(dictGP, dictR, tSt):
    print('--------------------------------------------------')
    print('Time step', tSt, '- Population Fitness', dictR['popFitness'],
          '- Genotype dictionary:')
    pprint.pprint(dictGP)

def writeGPResult(dictI, dictGP, dictR, cRp, tSt, cDir):
    fRes = open(cDir + '/' + dictI['nD_Results'] + '/' +
                dictI['nF_GPRes'] + '_Rep' + str(cRp) + '_timeStep' +
                str(tSt) + dictI['nF_End'], 'w')
    fRes.write('Count Allele AlleleCharacteristic AlleleFitness' +
               ' AlleleProportion MarginalFitness Emerged MutationType' +
               ' SequenceDiffAbsPerc SequenceDiffPropWtPerc\n')
    lAll = sortAlleles(dictGP, dictR['lAllSeqsPers'],
                       dictI['dAllSortM']['GPcur'])
    dictSeqDiffGP = getSeqDiffMatrix(dictGP)
    for curSeq in lAll:
        numSeqDiffPerc, wtSeqDiffPerc = \
            getNumberOfSeqDiffs(dictI, dictGP, dictSeqDiffGP, lAll, curSeq)
        fRes.write(str(lAll.index(curSeq) + 1) + ' ')
        fRes.write(str(dictGP[curSeq][0]) + ' ' + curSeq + ' ')
        for oneEl in dictGP[curSeq][1:]:
            fRes.write(str(oneEl) + ' ')
        fRes.write(str(numSeqDiffPerc) + ' ' + str(wtSeqDiffPerc) + '\n')
    fRes.close()

def writeFinalData(dictI, dictGTF, dictR, cRp, cDir):
    dFin_S = {'numZygotes': {'numHom': 0, 'numHet': 0},
              'hetAdv': {'avFHom': 0.0, 'avFHet': 0.0, 'curHetAdv': 0.0}}
    dFin_E = {'numZygotes': {'numHom': 0, 'numHet': 0},
              'hetAdv': {'avFHom': 0.0, 'avFHet': 0.0, 'curHetAdv': 0.0}}
    for oneKey, oneVal in dictGTF.items():
        if oneKey[0] == oneKey[1]:  # homozygote
            dFin_S['numZygotes']['numHom'] += 1
            avFitHom = calcRunMean(oneVal, dFin_S['hetAdv']['avFHom'],
                                   dFin_S['numZygotes']['numHom'])
            dFin_S['hetAdv']['avFHom'] = avFitHom
            if oneKey[0] in dictR['lAllIDsPers']:
                dFin_E['numZygotes']['numHom'] += 1
                avFitHom = calcRunMean(oneVal, dFin_E['hetAdv']['avFHom'],
                                       dFin_E['numZygotes']['numHom'])
                dFin_E['hetAdv']['avFHom'] = avFitHom
        else:                       # heterozygote
            dFin_S['numZygotes']['numHet'] += 1
            avFitHet = calcRunMean(oneVal, dFin_S['hetAdv']['avFHet'],
                                   dFin_S['numZygotes']['numHet'])
            dFin_S['hetAdv']['avFHet']= avFitHet
            if (oneKey[0] in dictR['lAllIDsPers'] and
                oneKey[1] in dictR['lAllIDsPers']):
                dFin_E['numZygotes']['numHet'] += 1
                avFitHet = calcRunMean(oneVal, dFin_E['hetAdv']['avFHet'],
                                       dFin_E['numZygotes']['numHet'])
                dFin_E['hetAdv']['avFHet']= avFitHet
    startHetAdv = dFin_S['hetAdv']['avFHet'] - dFin_S['hetAdv']['avFHom']
    dFin_S['hetAdv']['curHetAdv'] = startHetAdv
    endHetAdv = dFin_E['hetAdv']['avFHet'] - dFin_E['hetAdv']['avFHom']
    dFin_E['hetAdv']['curHetAdv'] = endHetAdv
    if dictI['dOutput']['level'] == O_LEVEL_HIGH:
        pprint.pprint(dictGTF)
        print('Start # homozygotes:', dFin_S['numZygotes']['numHom'])
        print('Start # heterozygotes:', dFin_S['numZygotes']['numHet'])
        print('Start avFit homozygotes:', dFin_S['hetAdv']['avFHom'])
        print('Start avFit heterozygotes:', dFin_S['hetAdv']['avFHet'])
        print('Start heterozygote advantage:', startHetAdv)
        print('End # homozygotes:', dFin_E['numZygotes']['numHom'])
        print('End # heterozygotes:', dFin_E['numZygotes']['numHet'])
        print('End avFit homozygotes:', dFin_E['hetAdv']['avFHom'])
        print('End avFit heterozygotes:', dFin_E['hetAdv']['avFHet'])
        print('End heterozygote advantage:', endHetAdv)
    fRes = open(cDir + '/' + dictI['nD_Results'] + '/' + dictI['nF_HetAdv'] +
                '_Rep' + str(cRp) + dictI['nF_End'], 'w')
    fRes.write('StartNumHomozygotes StartNumHeterozygotes' +
               ' StartAvFitnessHomozygotes StartAvFitnessHeterozygotes' +
               ' StartHeterozygoteAdvantage' +
               ' EndNumHomozygotes EndNumHeterozygotes' +
               ' EndAvFitnessHomozygotes EndAvFitnessHeterozygotes' +
               ' EndHeterozygoteAdvantage\n')
    fRes.write(str(dFin_S['numZygotes']['numHom']) + ' ' +
               str(dFin_S['numZygotes']['numHet']) + ' ' +
               str(dFin_S['hetAdv']['avFHom']) + ' ' +
               str(dFin_S['hetAdv']['avFHet']) + ' ' +
               str(startHetAdv) + ' ' +
               str(dFin_E['numZygotes']['numHom']) + ' ' +
               str(dFin_E['numZygotes']['numHet']) + ' ' +
               str(dFin_E['hetAdv']['avFHom']) + ' ' +
               str(dFin_E['hetAdv']['avFHet']) + ' ' +
               str(endHetAdv) + '\n')
    fRes.close()

def AnalyseGP(dictI, dictGP, cRp, cDir):
    dictSeqDiffGP = getSeqDiffMatrix(dictGP)
    fRes = open(cDir + '/' + dictI['nD_Results'] + '/' + dictI['nF_DivRes'] +
                'GPFinal'+ '_Rep' + str(cRp) + dictI['nF_End'], 'w')
    fRes.write('Allele1 Allele2 NumberOfDifferences PercentDifferent\n')
    lDictItems = sorted(dictSeqDiffGP)
    for curAllPair in lDictItems:
        numDiff = dictSeqDiffGP[curAllPair]
        fRes.write(str(curAllPair[0]) + ' ' + str(curAllPair[1]) + ' ' +
                   str(numDiff) + ' ' +
                   str(round(numDiff/dictI['numAACpA']*100.0, 2)) + '\n')
    fRes.close()

def removeCurRepFiles(dictI, cRp, cDir):
    for nOneF in os.listdir(cDir + '/' + dictI['nD_Results']):
        lFBits = nOneF.split('.')[0].split('_')
        for oneBit in lFBits:
            if oneBit[:3] == 'Rep':
                if oneBit[3:] == str(cRp):
                    os.remove(cDir + '/' + dictI['nD_Results'] + '/' + nOneF)

def checkNumAllPers(dI, dR, dRR, cRp, fSucc):
    print('# alleles persisting:', dR['nAllPers'])
    if dR['nAllPers'] >= dI['numAllThr'] + dI['dClOutc']['huge']:
        dRR['clOutcInc']['huge'] += 1
        if fSucc:
            addEntryToDVList(dRR['dHugeNAllPers'], dR['nAllPers'], cRp)
    elif dR['nAllPers'] >= dI['numAllThr'] + dI['dClOutc']['large']:
        dRR['clOutcInc']['large'] += 1
    elif dR['nAllPers'] <= dI['numAllThr'] + dI['dClOutc']['wee']:
        dRR['clOutcInc']['wee'] += 1
    elif dR['nAllPers'] <= dI['numAllThr'] + dI['dClOutc']['small']:
        dRR['clOutcInc']['small'] += 1
    else:
        dRR['clOutcInc']['normal'] += 1
    dRR['dAllPers'][cRp] = dR['nAllPers']

def writeResults(dictI, dictGP, dictGPT, dictGTF, dictR, dictRR, curRp, tStep,
                 curDir):
    fGoOn = True
    if (dictR['nAllPers'] >= dictI['numAllThr'] + dictI['dClOutc']['huge'] or
        dictI['dFlags']['breakEarly'] == False):
        if tStep%dictI['modDisp'] == 0:
            numTSRep = dictI['numIter'] + dictI['numARIter']
            numTSTot = dictI['numRep']*numTSRep
            percRep = round(tStep/numTSRep*100.0, 2)
            percTot = round(((curRp - 1)*numTSRep + tStep)/numTSTot*100.0, 2)
            print('_Repetition ', curRp, ' of ', dictI['numRep'],
                  ' - time step ', tStep, ' of ', numTSRep, ' (', percRep,
                  '% rep / ', percTot, '% tot)_', sep = '')
            print('num. alleles (cur.): (', len(dictGP), ', ',
                  dictR['nAllPers'], ')',
                  ', num. alleles (total): ', dictR['numAllTot'], sep = '')
            if 'numSglPtMutTot' in dictR:
                print('single point mut.:', dictR['numSglPtMutTot'])
            if 'numAllPtMutTot' in dictR:
                print('allele point mut.:', dictR['numAllPtMutTot'])
            if 'numMicConvMutTot' in dictR:
                print('micro-conversion mut.:', dictR['numMicConvMutTot'])
            if 'numNewSeqMutTot' in dictR:
                print('new sequence mut.:', dictR['numNewSeqMutTot'])
            print('Length of dictGTF = ', len(dictGTF),
                  ', length of dictGPT = ', len(dictGPT), '.', sep ='')
        if (tStep%dictI['modRate'] == 0 or tStep%dictI['modFile'] == 0):
            updateDictRes(dictI, dictGP, dictR)
        if tStep%dictI['modRate'] == 0:
            writeToRateResult(dictI, dictR, curRp, tStep, curDir)
        if tStep%dictI['modFile'] == 0:
            if dictI['dOutput']['level'] == O_LEVEL_HIGH:
                printGPResult(dictGP, dictR, tStep)
            writeGPResult(dictI, dictGP, dictR, curRp, tStep, curDir)
        if tStep == dictI['numIter']:
            writeFinalData(dictI, dictGTF, dictR, curRp, curDir)
            if dictI['dOutput']['level'] >= O_LEVEL_STD:
                AnalyseGP(dictI, dictGP, curRp, curDir)
            checkNumAllPers(dictI, dictR, dictRR, curRp, fGoOn)
    else:
        fGoOn = False
        checkNumAllPers(dictI, dictR, dictRR, curRp, fGoOn)
        removeCurRepFiles(dictI, curRp, curDir)
    return fGoOn

def writeRepRes(dictI, dictRR, cDir):
    fCur = open(cDir + '/' + dictI['nD_Results'] + '/' +
                dictI['nF_RepRes'] + dictI['nF_End'], 'w')
    lClOutc = sorted(dictI['dClOutc'])
    for oneKey in lClOutc[:-1]:
        fCur.write(oneKey + ' ')
    fCur.write(lClOutc[-1] + '\n')
    for oneKey in lClOutc[:-1]:
        fCur.write(str(dictRR['clOutcInc'][oneKey]) + ' ')
    fCur.write(str(dictRR['clOutcInc'][lClOutc[-1]]) + '\n')
    pprint.pprint(dictRR['dHugeNAllPers'])
    fCur.close()
