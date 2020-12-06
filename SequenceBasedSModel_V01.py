# -*- coding: utf-8 -*-
################################################################################
# SequenceBasedSModel_V01.py #
################################################################################
#!/opt/local/Library/Frameworks/Python.framework/Versions/3.2

import os, sys, random, time

lSubDirs = ['SBMSupport']
for oneSubDir in lSubDirs:
    sys.path.append(sys.path[0] + '/' + oneSubDir)

import InputVariables
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
from SpecialFunctions import (genDictInp, genDictResRepRes, resetDictRes,
                              genDictGenepools, calcAllProp,
                              genDictGenotypeFit, updateDicts,
                              oneTimeStepDynamics)
from GenericFunctions import zipFolder
from WriteAndPrint import (deleteAllFiles, printCurRepet, openOutFiles,
                           writeResults, writeRepRes)
from AnalyseRawData import analyseRawDataFiles, mergeRepRes

# start timing
InputVariables.startTime = time.time()

print('SequencebasedQModel.py')
print('Simulates the evolution of a gene in a population over time',
      'maximising marginal allele fitness')

# MAIN PROGRAM #
actDir = os.getcwd()
print('Actual directory:', actDir)
random.seed()
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
dictRes, dictRepRes = genDictResRepRes(dictInp)
deleteAllFiles(dictInp, actDir)
InputVariables.elT_Init += time.time() - InputVariables.startTime
for curRep in range(1, dictInp['numRep'] + 1):
    InputVariables.stT_Init = time.time()
    resetDictRes(dictInp, dictRes)
    printCurRepet(curRep, dictInp, dictRepRes)
    dictGeneP, dictGenePTot = genDictGenepools(dictInp)
    calcAllProp(dictInp, dictGeneP)
    dictGenTypeFit = genDictGenotypeFit(dictInp, dictGeneP)
    updateDicts(dictInp, dictGeneP, dictGenTypeFit, dictRes)
    openOutFiles(dictInp, dictGeneP, curRep, actDir)
    InputVariables.elT_Init += time.time() - InputVariables.stT_Init
    for timeStep in range(numIterations + numARIter + 1):
        InputVariables.stT_TSDyn = time.time()
        if timeStep > 0:
            oneTimeStepDynamics(dictInp, dictGeneP, dictGenePTot,
                                dictGenTypeFit, dictRes, timeStep)
        InputVariables.elT_TSDyn += time.time() - InputVariables.stT_TSDyn
        InputVariables.stT_Write = time.time()
        goOn = writeResults(dictInp, dictGeneP, dictGenePTot, dictGenTypeFit,
                            dictRes, dictRepRes, curRep, timeStep, actDir)
        InputVariables.elT_Write += time.time() - InputVariables.stT_Write
        if not goOn:
            break
    InputVariables.stT_DatAn = time.time()
    analyseRawDataFiles(dictInp, curRep, actDir)
    InputVariables.elT_DatAn += time.time() - InputVariables.stT_DatAn
InputVariables.stT_DatAn = time.time()
if dictInp['numRep'] > 1:
    mergeRepRes(dictInp, actDir, dictInp['numRep'])
writeRepRes(dictInp, dictRepRes, actDir)
#zipFolder(dictInp['nD_Results'], dictInp['nD_Results'], True)
InputVariables.elT_DatAn += time.time() - InputVariables.stT_DatAn
print('-----------\nelT_Init =', InputVariables.elT_Init, '\nelT_TSDyn =',
      InputVariables.elT_TSDyn, '\nelT_Write =', InputVariables.elT_Write,
      '\nelT_DatAn =', InputVariables.elT_DatAn, '\n-----------')
print('elT_Mut =', InputVariables.elT_Mut, '\nelT_PropUp =',
      InputVariables.elT_PropUp, '\nelT_FitUp =', InputVariables.elT_FitUp,
      '\nelT_CleanD =', InputVariables.elT_CleanD, '\n-----------')
elapsedTime = round(time.time() - InputVariables.startTime, 2)
print('Total time elapsed: ', elapsedTime, ' seconds (',
      round(elapsedTime/numIterations*1000, 2),
      ' seconds per 1000 time steps).\nALL RUNS FINISHED.\n', sep = '')
