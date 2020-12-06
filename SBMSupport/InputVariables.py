# -*- coding: utf-8 -*-
# ### INPUTS ###################################################################

from Constants import *

# other constants
NUM_REPETITIONS = 3
NUM_ITERATIONS = 10000000   # 10000000 (10 million generations / 40m years)
NUM_AFTERRUN_ITER = round(NUM_ITERATIONS/10) # NUM_ITERATIONS/10 (no mut.,...)
POP_SIZE = 100000          # 1000000
NUM_ALL_INI = 1             # 1
NUM_AMINOAC_P_ALL = 82      # 82
LOST_THRESHOLD = 1/(2*2*POP_SIZE)   # 1/(2*2*POP_SIZE)
DISREG_THRESHOLD = LOST_THRESHOLD   # LOST_THRESHOLD
NUM_ALL_THRESHOLD = 100     # 100
MIN_ALL_FIT = 0.0           # 0.0
MAX_ALL_FIT = 0.5           # 0.5 (0.47)
MIN_GT_FIT = 0.0            # 0.0
MAX_GT_FIT = 1.0            # 1.0 (0.47)
#DELTA_IN = (MAX_GT_FIT - MAX_ALL_FIT)/NUM_AMINOAC_P_ALL    # or a fixed value
DELTA_IN = 0.08/NUM_AMINOAC_P_ALL   # 0.24/NUM_AMINOAC_P_ALL
GAUSS_MU_IN = MIN_ALL_FIT + (MAX_ALL_FIT - MIN_ALL_FIT)/2   # /2
#GAUSS_SIG_IN = (MAX_ALL_FIT - MIN_ALL_FIT)/25              # /25
GAUSS_SIG_IN = 0.01         # 0.02
DICT_DEF_FIT = {0: GAUSS_MU_IN} # {0: GAUSS_MU_IN}
P_PT_MUT_SGLAA = 1.2195E-08 # 1.2195E-08 (-> pt. mutation on exon 2 = 1.0E-06)
P_MIC_CONV_MUT = 0.0        # 1.0E-08 (prob. for micro-conversion)
P_NEW_SEQ_MUT = 0.0         # 0.0 (prob. for new sequence mutation)
LEN_MIC_CONV_SEQ_MU = 15    # 15
LEN_MIC_CONV_SEQ_SIG = 5    # 5
HET_ADV_MODE = HAM_CONST_DELTA   # HAM_CONST_DELTA
LAMBDA_DEBOER = 1           # 1
FIXED_HET_ADV = 0.08        # 0.24
MOD_DISP = round(NUM_ITERATIONS/100)    # NUM_ITERATIONS/100
MOD_STEP = round(NUM_ITERATIONS/1000)   # NUM_ITERATIONS/1000
MOD_RATE = round(NUM_ITERATIONS/1000)   # NUM_ITERATIONS/1000
MOD_FILE = round(NUM_ITERATIONS/100)    # NUM_ITERATIONS/100
D_GTF_CLEAN_TS = 100        # 100
D_ALL_SORT_MODE = {'GPcur': M_MOST_COMMON_ALLELE}   # M_MOST_COMMON_ALLELE
D_FLAGS = {'breakEarly': False,
           'mutFitDepend': True}
D_MUT_INFO = {'pPtMutSglAA': P_PT_MUT_SGLAA,
              'pMicConvMut': P_MIC_CONV_MUT,
              'pNewSeqMut': P_NEW_SEQ_MUT,
              'newFitVarUniR': (MAX_ALL_FIT - MIN_ALL_FIT)/4.0,
              'newFitVarGauss': GAUSS_SIG_IN,
              'muLMicConvSeq': LEN_MIC_CONV_SEQ_MU,
              'sigLMicConvSeq': LEN_MIC_CONV_SEQ_SIG,
              'dGTFCleanTS': D_GTF_CLEAN_TS}
# values in lists of special analysis dictionary have to be ascending
D_SPEC_ANA = {'shareOfGP': [1.0E-06, 1.0E-05, 1.0E-04, 1.0E-03, 1.0E-02],
              'copiesOfAllThr': [10, 25, 100, 250, 1000],
              'ageThrNumGen': [10, 25, 100, 250, 1000]}
D_CLASS_OUTC = {'wee': -3,
                'small': -1,
                'normal': 0,
                'large': 1,
                'huge': 3}
D_OUTPUT = {'level': O_LEVEL_STD,
            'numAgeCLasses': N_AGE_CL_100}
ND_RESULTS = '_Results'
NF_RATERES = 'RateResult'
NF_ANARES = 'RawDataAnalysisResult'
NF_SPECRES = 'SpecialAnalysisResult'
NF_AGERES = 'AgeAnalysisResult'
NF_ALLRES = 'AlleleComparisonResult'
NF_PHYL = 'PhyloTreeCoordinates'
NF_GPRES = 'GenepoolResult'
NF_SGPRES = 'ShortSortedGenepoolResult'
NF_TGPRES = 'TotalGenepoolResult'
NF_HA = 'HeterozygoteAdvantageResult'
NF_DIVRES = 'DiversityResult'
NF_REPRES = '_RepetitionResult'
NF_END = '.dat'

# initiate variables
iniAllFDist = A_INI_FDIST_DEFINED   # A_INI_FDIST_DEFINED
iniAllProp = A_PROP_UNIDET          # A_PROP_UNIDET
newAllFDist = A_NEW_FDIST_GAUSS     # A_NEW_FDIST_GAUSS
limFitness = L_FIT_MAX              # L_FIT_MAX
dominType = D_SEQ_B_MAXPD_W_HA      # D_SEQ_B_MAXPD_W_HA

deltaIn = DELTA_IN                  # het. adv. per single sequence difference
gaussMuIn = GAUSS_MU_IN
gaussSigIn = GAUSS_SIG_IN
minFitn = MIN_ALL_FIT
maxFitn = MAX_ALL_FIT
minGTFitn = MIN_GT_FIT
maxGTFitn = MAX_GT_FIT
dictDefFit = DICT_DEF_FIT
lostThr = LOST_THRESHOLD
disregThr = DISREG_THRESHOLD
numAllThr = NUM_ALL_THRESHOLD
numRepetitions = NUM_REPETITIONS
numIterations = NUM_ITERATIONS
numARIter = NUM_AFTERRUN_ITER
popSize = POP_SIZE
numAllI = NUM_ALL_INI
numAmAcPerAll = NUM_AMINOAC_P_ALL
hetAdvM = HET_ADV_MODE
lambdaDeB = LAMBDA_DEBOER
fixedHA = FIXED_HET_ADV
modDisp = MOD_DISP
modRate = MOD_RATE
modFile = MOD_FILE
dAllSortM = D_ALL_SORT_MODE
dFlags = D_FLAGS
dMutInfo = D_MUT_INFO
dSpecAna = D_SPEC_ANA
dClOutc = D_CLASS_OUTC
dOutput = D_OUTPUT

# timings
startTime, elapsedTime = 0.0, 0.0
stT_Init, elT_Init = 0.0, 0.0
stT_TSDyn, elT_TSDyn = 0.0, 0.0
stT_Write, elT_Write = 0.0, 0.0
stT_DatAn, elT_DatAn = 0.0, 0.0
stT_Mut, elT_Mut = 0.0, 0.0
stT_PropUp, elT_PropUp = 0.0, 0.0
stT_FitUp, elT_FitUp = 0.0, 0.0
stT_CleanD, elT_CleanD = 0.0, 0.0
