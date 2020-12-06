# -*- coding: utf-8 -*-
###############################################################################
# PandasRawDataAnalysis.py #
###############################################################################
import os

import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.graphics.regressionplots as regplt

# --- CONSTANTS ---------------------------------------------------------------
P_RAW_DB = os.path.join('..', '..', '62_AlteProj_Glasgow', '4_AllSimulations',
                        '5_SequenceBasedModel', 'Sims')
S_RES_INP = '_Results'

S_C_ALL_FIT = 'AlleleFitness'
S_C_ALL_PROP = 'AlleleProportion'
S_C_EMGD = 'Emerged'
S_C_SQ_D_ABS = 'SequenceDiffAbsPerc'
S_C_SQ_D_PWT = 'SequenceDiffPropWtPerc'

S_C_ALL_FIT_S = 'AllFit'
S_C_ALL_PROP_S = 'AllProp'
S_C_EMGD_S = 'Emgd'
S_C_SQ_D_ABS_S = 'SqDfAbs'
S_C_SQ_D_PWT_S = 'SqDfPWt'

S_FILE = 'File'
S_NM_PLT_F = 'sNmPltF'
S_GP_RES = 'GenepoolResult'
S_REP = 'Rep'
S_TIME_STEP = 'timeStep'
T_S_F_SPL = (S_GP_RES, S_REP, S_TIME_STEP)
S_MIN_ALL_AGE = 'minAge'

S_SLOPE = 'Slope'
S_P_VAL_SL = 'pValue'
S_R2 = 'R2'
S_AIC = 'AIC'

SP_CSV = ';'
SP_S = '_'
SP_DOT = '.'
SP_N = ' '

S_PDF = 'pdf'
S_CSV = 'csv'

# --- DERIVED CONSTANTS -------------------------------------------------------
D_S_COL = {S_C_ALL_FIT: S_C_ALL_FIT_S,
           S_C_ALL_PROP: S_C_ALL_PROP_S,
           S_C_EMGD: S_C_EMGD_S,
           S_C_SQ_D_ABS: S_C_SQ_D_ABS_S,
           S_C_SQ_D_PWT: S_C_SQ_D_PWT_S}

# --- DIRECTORIES -------------------------------------------------------------
# set 1 (new)
LD1 = ['831S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_10kInd_10kGen',
       '715S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_25kInd_10kGen',
       '716S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_100kInd_10kGen',
       '717S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_250kInd_10kGen',
       '718S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_1mInd_10kGen',
       '719S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_2p5mInd_10kGen',
       '922S_Sto_3rdMod_PtMutOnly_STD_Delt008_Sig001_10mInd_10kGen',
       '837S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_10kInd_10kGen',
       '720S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_25kInd_10kGen',
       '721S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_100kInd_10kGen',
       '722S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_250kInd_10kGen',
       '723S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_1mInd_10kGen',
       '724S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_2p5mInd_10kGen',
       '924S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0075_Sig001_10mInd_10kGen']

# set 2 (new)
LD2 = ['800S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_10kInd_10kGen',
       '801S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_25kInd_10kGen',
       '802S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_100kInd_10kGen',
       '803S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_250kInd_10kGen',
       '804S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_1mInd_10kGen',
       '805S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_2p5mInd_10kGen',
       '740S_Sto_3rdMod_PtMutOnly_STD_Delt024_Sig002_10mInd_10kGen',
       '806S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_10kInd_10kGen',
       '807S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_25kInd_10kGen',
       '808S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_100kInd_10kGen',
       '809S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_250kInd_10kGen',
       '810S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_1mInd_10kGen',
       '811S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_2p5mInd_10kGen',
       '743S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA024_Sig002_10mInd_10kGen']

# set 3 (new)
LD3 = ['830S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_10kInd_10kGen',
       '705S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_25kInd_10kGen',
       '706S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_100kInd_10kGen',
       '707S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_250kInd_10kGen',
       '708S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_1mInd_10kGen',
       '709S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_2p5mInd_10kGen',
       '921S_Sto_3rdMod_PtMutOnly_STD_Delt021_Sig0005_10mInd_10kGen',
       '836S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_10kInd_10kGen',
       '710S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_25kInd_10kGen',
       '711S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_100kInd_10kGen',
       '712S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_250kInd_10kGen',
       '713S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_1mInd_10kGen',
       '714S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_2p5mInd_10kGen',
       '923S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0204_Sig0005_10mInd_10kGen']

# set 4 (new)
LD4 = ['949S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_10kInd_10kGen',
       '950S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_25kInd_10kGen',
       '951S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_100kInd_10kGen',
       '952S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_250kInd_10kGen',
       '946S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_1mInd_10kGen',
       '953S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_2p5mInd_10kGen',
       '954S_Sto_3rdMod_PtMutOnly_STD_Delt00125_Sig00015_10mInd_10kGen',
       '955S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_10kInd_10kGen',
       '956S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_25kInd_10kGen',
       '957S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_100kInd_10kGen',
       '958S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_250kInd_10kGen',
       '947S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_1mInd_10kGen',
       '959S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_2p5mInd_10kGen',
       '960S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0011_Sig00015_10mInd_10kGen']

# set 5 (new)
LD5 = ['888S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_10kInd_10kGen',
       '889S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_25kInd_10kGen',
       '890S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_100kInd_10kGen',
       '891S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_250kInd_10kGen',
       '879S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_1mInd_10kGen',
       '880S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_2p5mInd_10kGen',
       '912S_Sto_3rdMod_PtMutOnly_STD_Delt00145_Sig0005_10mInd_10kGen',
       '869S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_10kInd_10kGen',
       '870S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_25kInd_10kGen',
       '871S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_100kInd_10kGen',
       '872S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_250kInd_10kGen',
       '864S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_1mInd_10kGen',
       '867S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_2p5mInd_10kGen',
       '913S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0013_Sig0005_10mInd_10kGen']

# set 6 (new)
LD6 = ['892S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_10kInd_10kGen',
       '893S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_25kInd_10kGen',
       '894S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_100kInd_10kGen',
       '895S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_250kInd_10kGen',
       '881S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_1mInd_10kGen',
       '882S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_2p5mInd_10kGen',
       '914S_Sto_3rdMod_PtMutOnly_STD_Delt0021_Sig002_10mInd_10kGen',
       '873S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_10kInd_10kGen',
       '874S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_25kInd_10kGen',
       '875S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_100kInd_10kGen',
       '876S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_250kInd_10kGen',
       '877S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_1mInd_10kGen',
       '878S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_2p5mInd_10kGen',
       '915S_Sto_FxHAAv3rdMod_PtMutOnly_STD_fHA0019_Sig002_10mInd_10kGen']

dSet = {1: LD1, 2: LD2, 3: LD3, 4: LD4, 5: LD5, 6: LD6}

# --- NAME COMPONENTS RELATED TO DIRECTORIES ----------------------------------
dSNmCmpDir = {# Set1 / DAA
              '831S': ('Set1', 'A', '01'),
              '715S': ('Set1', 'A', '02'),
              '716S': ('Set1', 'A', '03'),
              '717S': ('Set1', 'A', '04'),
              '718S': ('Set1', 'A', '05'),
              '719S': ('Set1', 'A', '06'),
              '922S': ('Set1', 'A', '07'),
              # Set1 / AO
              '837S': ('Set1', 'B', '01'),
              '720S': ('Set1', 'B', '02'),
              '721S': ('Set1', 'B', '03'),
              '722S': ('Set1', 'B', '04'),
              '723S': ('Set1', 'B', '05'),
              '724S': ('Set1', 'B', '06'),
              '924S': ('Set1', 'B', '07'),
              # Set2 / DAA
              '800S': ('Set2', 'A', '01'),
              '801S': ('Set2', 'A', '02'),
              '802S': ('Set2', 'A', '03'),
              '803S': ('Set2', 'A', '04'),
              '804S': ('Set2', 'A', '05'),
              '805S': ('Set2', 'A', '06'),
              '740S': ('Set2', 'A', '07'),
              # Set2 / AO
              '806S': ('Set2', 'B', '01'),
              '807S': ('Set2', 'B', '02'),
              '808S': ('Set2', 'B', '03'),
              '809S': ('Set2', 'B', '04'),
              '810S': ('Set2', 'B', '05'),
              '811S': ('Set2', 'B', '06'),
              '743S': ('Set2', 'B', '07'),
              # Set3 / DAA
              '830S': ('Set3', 'A', '01'),
              '705S': ('Set3', 'A', '02'),
              '706S': ('Set3', 'A', '03'),
              '707S': ('Set3', 'A', '04'),
              '708S': ('Set3', 'A', '05'),
              '709S': ('Set3', 'A', '06'),
              '921S': ('Set3', 'A', '07'),
              # Set3 / AO
              '836S': ('Set3', 'B', '01'),
              '710S': ('Set3', 'B', '02'),
              '711S': ('Set3', 'B', '03'),
              '712S': ('Set3', 'B', '04'),
              '713S': ('Set3', 'B', '05'),
              '714S': ('Set3', 'B', '06'),
              '923S': ('Set3', 'B', '07'),
              # Set4 / DAA
              '949S': ('Set4', 'A', '01'),
              '950S': ('Set4', 'A', '02'),
              '951S': ('Set4', 'A', '03'),
              '952S': ('Set4', 'A', '04'),
              '946S': ('Set4', 'A', '05'),
              '953S': ('Set4', 'A', '06'),
              '954S': ('Set4', 'A', '07'),
              # Set4 / AO
              '955S': ('Set4', 'B', '01'),
              '956S': ('Set4', 'B', '02'),
              '957S': ('Set4', 'B', '03'),
              '958S': ('Set4', 'B', '04'),
              '947S': ('Set4', 'B', '05'),
              '959S': ('Set4', 'B', '06'),
              '960S': ('Set4', 'B', '07'),
              # Set5 / DAA
              '888S': ('Set5', 'A', '01'),
              '889S': ('Set5', 'A', '02'),
              '890S': ('Set5', 'A', '03'),
              '891S': ('Set5', 'A', '04'),
              '879S': ('Set5', 'A', '05'),
              '880S': ('Set5', 'A', '06'),
              '912S': ('Set5', 'A', '07'),
              # Set5 / AO
              '869S': ('Set5', 'B', '01'),
              '870S': ('Set5', 'B', '02'),
              '871S': ('Set5', 'B', '03'),
              '872S': ('Set5', 'B', '04'),
              '864S': ('Set5', 'B', '05'),
              '867S': ('Set5', 'B', '06'),
              '913S': ('Set5', 'B', '07'),
              # Set6 / DAA
              '892S': ('Set6', 'A', '01'),
              '893S': ('Set6', 'A', '02'),
              '894S': ('Set6', 'A', '03'),
              '895S': ('Set6', 'A', '04'),
              '881S': ('Set6', 'A', '05'),
              '882S': ('Set6', 'A', '06'),
              '914S': ('Set6', 'A', '07'),
              # Set6 / AO
              '873S': ('Set6', 'B', '01'),
              '874S': ('Set6', 'B', '02'),
              '875S': ('Set6', 'B', '03'),
              '876S': ('Set6', 'B', '04'),
              '877S': ('Set6', 'B', '05'),
              '878S': ('Set6', 'B', '06'),
              '915S': ('Set6', 'B', '07')}

# --- INPUT -------------------------------------------------------------------
usedSet = 1

dXInfo = {S_MIN_ALL_AGE: 100000}   # 0 / 100 / 1000 / 100000

sResOLS = 'ResultOLS'

sCAAge = 'AlleleAge'
sNmPltF = sCAAge

pPlots = '_Visualisations'
pResDfr = '_ResultTable'

# lTSUsed = (list(range(100000, 1000000, 100000)) +
#            list(range(1000000, 11000001, 1000000)))
lTSUsed = list(range(1000000, 11000001, 1000000))

# lSXOLS = [S_C_ALL_FIT, S_C_ALL_PROP, S_C_SQ_D_ABS]
# lSXOLS = [S_C_ALL_FIT, S_C_SQ_D_ABS]
lSXOLS = [S_C_SQ_D_ABS]      # S_C_SQ_D_ABS / S_C_SQ_D_PWT / S_C_ALL_FIT

# --- PLOT INPUT --------------------------------------------------------------
title_PltXY = None                           # title of plot
xLbl_PltXY = None                            # x-label of plot
yLbl_PltXY = 'Allele age'                    # y-label of plot
xLimB = None            # bottom x-limit of plot (None / number)
xLimT = None            # top x-limit of plot (None / number)
yLimB = None            # bottom y-limit of plot (None / number (e.g. 16))
yLimT = None            # top y-limit of plot (None / number)
tpMark_PltXY = 'x'      # marker type of plot
szMark_PltXY = 2        # marker size of plot
ewMark_PltXY = 1        # marker edge width of plot
styLn_PltXY = ''        # line style of plot
wdthLn_PltXY = 1        # line width of plot
lClr_PltXY = [(0.9, 0., 0.), (0.8, 0.4, 0.), (0.6, 0.6, 0.),
              (0., 0.8, 0.), (0., 0.4, 0.4), (0., 0., 0.9),
              (0.5, 0., 0.5), (0.25, 0.25, 0.25), (0.5, 0.5, 0.5),
              (0.75, 0.75, 0.75), (1., 0.4, 0.4), (1., 0.7, 0.),
              (1., 1., 0.), (0.3, 1., 0.3), (0., 1., 1.),
              (0.5, 0.5, 1.), (1., 0., 1.)]

# --- DERIVED INPUT -----------------------------------------------------------
lSD = dSet[usedSet]
tSHdRes = (S_FILE, S_TIME_STEP)
sStNm = SP_S.join([D_S_COL[sX] for sX in lSXOLS])
if lSXOLS == [S_C_ALL_FIT]:
    xLbl_PltXY = 'Intrinsic merit of allele'
elif lSXOLS == [S_C_SQ_D_ABS] or lSXOLS == [S_C_SQ_D_PWT]:
    xLbl_PltXY = 'Average sequence difference of allele'

# --- DERIVED PLOT INPUT ------------------------------------------------------
dPltXY = {S_NM_PLT_F: sNmPltF,
          'title': title_PltXY,
          'xLbl': xLbl_PltXY,
          'yLbl': yLbl_PltXY,
          'xLimB': xLimB,
          'xLimT': xLimT,
          'yLimB': yLimB,
          'yLimT': yLimT,
          'tpMark': tpMark_PltXY,
          'szMark': szMark_PltXY,
          'ewMark': ewMark_PltXY,
          'styLn': styLn_PltXY,
          'wdthLn': wdthLn_PltXY,
          'lClr': lClr_PltXY}

# --- FUNCTIONS ---------------------------------------------------------------
def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def joinToPath(pF = '', nmF = 'Dummy.txt'):
    if len(pF) > 0:
        createDir(pF)
        return os.path.join(pF, nmF)
    else:
        return nmF

def getInfoLSNmFRDRes(pRF, tInfoFSpl, cSep = SP_S, dtSep = SP_DOT):
    dInf = {}    # (repetition, timeStep): file_name (info dict.)
    for nmF in os.listdir(pRF):
        if nmF.startswith(tInfoFSpl[0]):
            # dInf[nmF] = [None, None]
            lSpl = nmF.split(dtSep)[0].split(cSep)
            assert len(lSpl) == len(tInfoFSpl)
            keyD = tuple([int(cSpl[len(tInfoFSpl[k + 1]):]) for k, cSpl in
                              enumerate(lSpl[1:])])
            dInf[keyD] = nmF
    return dInf

def addToDictL(cD, cK, cE):
    if cK in cD:
        cD[cK].append(cE)
    else:
        cD[cK] = [cE]

def getPR(pR, sCD, s0, dXI, dSNmC, sFXt = S_CSV, cSep = SP_S, dtSep = SP_DOT):
    sSt, sEnd = s0, dtSep + sFXt
    if S_MIN_ALL_AGE in dXI:
        if dXI[S_MIN_ALL_AGE] >  0:
            sSt += cSep + S_MIN_ALL_AGE + str(dXI[S_MIN_ALL_AGE])
    if S_TIME_STEP in dXI:
        if dXI[S_TIME_STEP] >=  0:
            sEnd = cSep + S_TIME_STEP + str(dXI[S_TIME_STEP]) + sEnd
    assert len(sCD) >= 4
    sK = sCD[:4]
    assert sK in dSNmC
    sClf = cSep.join(dSNmC[sK])
    return joinToPath(pR, sSt + cSep + sClf + cSep + sCD + sEnd)

def concPdDfrS(lPdDfr, concAx = 0, verInt = True, srt = False, ignIdx = False,
               dropAx = None):
    d = pd.concat(lPdDfr, axis = concAx, verify_integrity = verInt, sort = srt,
                  ignore_index = ignIdx)
    if dropAx in [0, 1, 'index', 'columns']:
        d.dropna(axis = dropAx, inplace = True)
    return d

# --- PLOT FUNCTIONS ----------------------------------------------------------
def pltXYAxis(dfr, nmCX = None, nmCY = None, pltAxXY = None):
    minDfr, maxDfr = min(0, dfr.stack().min()), dfr.stack().max()
    if pltAxXY is None:
        return pltAxXY
    assert len(pltAxXY) >= 2
    if pltAxXY[0]:
        if nmCX is not None:
            minX, maxX = min(0, min(dfr.loc[:, nmCX])), max(dfr.loc[:, nmCX])
            plt.plot([minX, maxX], [0, 0], lw = 0.75, color = 'k')
        else:
            plt.plot([0, dfr.shape[0] - 1], [0, 0], lw = 0.75, color = 'k')
    if pltAxXY[1]:
        if nmCY is not None:
            minY, maxY = min(0, min(dfr.loc[:, nmCY])), max(dfr.loc[:, nmCY])
            plt.plot([0, 0], [minY, maxY], lw = 0.75, color = 'k')
        else:
            plt.plot([0, 0], [minDfr, maxDfr], lw = 0.75, color = 'k')

def setXLim(xLim = (None, None)):
    assert len(xLim) >= 2
    if xLim[0] is None:
        if xLim[1] is not None:
            plt.xlim(top = xLim[1])
    else:
        if xLim[1] is None:
            plt.xlim(bottom = xLim[0])
        else:
            plt.xlim(xLim)

def setYLim(yLim = (None, None)):
    assert len(yLim) >= 2
    if yLim[0] is None:
        if yLim[1] is not None:
            plt.ylim(top = yLim[1])
    else:
        if yLim[1] is None:
            plt.ylim(bottom = yLim[0])
        else:
            plt.ylim(yLim)

def decorateSaveFigLegOut(pF, cFig, dfr = None, sTtl = None, xLbl = None,
                          yLbl = None, xLim = None, yLim = None, nmCX = None,
                          nmCY = None, cLeg = None, pltAxXY = None):
    if pltAxXY is not None:
        assert len(pltAxXY) >= 2
        if dfr is not None:
            pltXYAxis(dfr, nmCX, nmCY, pltAxXY = pltAxXY)
    if sTtl is not None:
        plt.title(sTtl)
    if xLbl is not None:
        plt.xlabel(xLbl)
    if yLbl is not None:
        plt.ylabel(yLbl)
    setXLim(xLim)
    setYLim(yLim)
    if cLeg is not None:
        cFig.savefig(pF, bbox_extra_artists = (cLeg,), bbox_inches = 'tight')
    else:
        cFig.savefig(pF)

def fitModelOLS(cDfr, lSX, sY, sWave = ' ~ ', sPlus = ' + '):
    sMdl = sY
    if len(lSX) == 0:
        sMdl += sWave + '1'
    else:
        sMdl += sWave + lSX[0]
        for sX in lSX[1:]:
            sMdl += sPlus + sX
    return smf.ols(sMdl, cDfr).fit()

def pltDfrXY(cDfr, dIPlt, pF = 'Hugo.pdf', cOff = 0, pltAxXY = None,
             cModel = None):
    assert cDfr.shape[1] > 1
    sTtl, xLbl, yLbl = dIPlt['title'], dIPlt['xLbl'], dIPlt['yLbl']
    tpMark, szMark, ewMark = dIPlt['tpMark'], dIPlt['szMark'], dIPlt['ewMark']
    styLn, wdthLn, lClr = dIPlt['styLn'], dIPlt['wdthLn'], dIPlt['lClr']
    xLim = (dIPlt['xLimB'], dIPlt['xLimT'])
    yLim = (dIPlt['yLimB'], dIPlt['yLimT'])
    if cModel is not None:
        cFig = regplt.abline_plot(model_results = cModel)
        cAx = cFig.axes[0]
    else:
        cFig, cAx = plt.subplots()
    for k in range(1, cDfr.shape[1]):
        cClr = lClr[(cOff + k - 1) % len(lClr)]
        cAx.plot(cDfr.iloc[:, 0], cDfr.iloc[:, k], marker = tpMark,
                 ms = szMark, mew = ewMark, mec = cClr, mfc = cClr, ls = styLn,
                 lw = wdthLn, color = cClr)
    decorateSaveFigLegOut(pF, cFig, cDfr, sTtl, xLbl, yLbl, xLim,
                          yLim, nmCX = cDfr.columns[0], nmCY = cDfr.columns[k],
                          pltAxXY = pltAxXY)
    plt.close()

# --- MAIN --------------------------------------------------------------------
for sD in lSD:
    print('Processing directory', sD, '...')
    pRawDC = os.path.join(P_RAW_DB, sD, S_RES_INP)
    dInfoSF = getInfoLSNmFRDRes(pRawDC, T_S_F_SPL)
    dRes = {cSHdRes: [] for cSHdRes in tSHdRes}
    # for i, sCSel in enumerate(L_S_C_SEL[-2:]):
    #     print('Calculating', sCSel, '...')
    for cTS in lTSUsed:
        dXInfo[S_TIME_STEP] = cTS
        dfrM = pd.DataFrame(columns = lSXOLS + [sCAAge], dtype = 'float64')
        for tInfo, sNmF in sorted(dInfoSF.items()):
            if tInfo[1] == cTS:
                dfrRD = pd.read_table(joinToPath(pRawDC, sNmF), sep = SP_N)
                dfrSelC = dfrRD.loc[:, lSXOLS + [S_C_EMGD]]
                dfrSelC[sCAAge] = cTS - dfrSelC[S_C_EMGD]
                dfrM = concPdDfrS([dfrM, dfrSelC.loc[:, lSXOLS + [sCAAge]]],
                                  verInt = False, ignIdx = True)
        dfrM[sCAAge] = dfrM[sCAAge].astype('int64')
        dfrM = dfrM[dfrM[sCAAge] >= dXInfo[S_MIN_ALL_AGE]]
        cModOLS = fitModelOLS(dfrM, lSXOLS, sCAAge)
        lRes = [sD, cTS]
        assert len(lRes) == len(tSHdRes)
        for j, cRes in enumerate(lRes):
            dRes[tSHdRes[j]].append(cRes)
        for k in range(1, len(cModOLS.params)):
            sKey = S_SLOPE + SP_S + D_S_COL[cModOLS.params.index[k]]
            addToDictL(dRes, sKey, cModOLS.params.iloc[k])
        for k in range(1, len(cModOLS.pvalues)):
            sKey = S_P_VAL_SL + SP_S + D_S_COL[cModOLS.pvalues.index[k]]
            addToDictL(dRes, sKey, cModOLS.pvalues.iloc[k])
        addToDictL(dRes, S_R2, cModOLS.rsquared)
        addToDictL(dRes, S_AIC, cModOLS.aic)
        pPltF = getPR(pPlots, sD, sStNm, dXInfo, dSNmCmpDir, S_PDF)
        if not os.path.isfile(pPltF) and len(lSXOLS) == 1:
            pltDfrXY(dfrM.loc[:, [lSXOLS[0], sCAAge]], dPltXY, pF = pPltF,
                     cOff = 0, pltAxXY = (True, True), cModel = cModOLS)
    dXInfo[S_TIME_STEP] = -1
    pResF = getPR(pResDfr, sD, sResOLS + SP_S + sStNm, dXInfo, dSNmCmpDir)
    if not os.path.isfile(pResF):
        dfrRes = pd.DataFrame(dRes)
        dfrRes.to_csv(pResF, sep = SP_CSV)
        print('Saved results for directory', sD, '...')
print('*'*8, 'DONE! (type', lSXOLS, '/ min. allele age',
      dXInfo[S_MIN_ALL_AGE], '/ set', usedSet, ')', '*'*8)

###############################################################################
