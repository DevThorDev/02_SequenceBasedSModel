# -*- coding: utf-8 -*-
# ### CONSTANTS ################################################################
# ### flag constants ###
# different initial fitness distributions
A_INI_FDIST_UNIDET = 1
A_INI_FDIST_UNIRAND = 2
A_INI_FDIST_GAUSS = 3
A_INI_FDIST_DEFINED = 4

# different new (mutated) allele fitness distributions
A_NEW_FDIST_UNIRAND = 2
A_NEW_FDIST_GAUSS = 3

# different initial proportions of alleles
A_PROP_UNIDET = 1
A_PROP_UNIRAND = 2

# different fitness maxima
L_FIT_1 = 1
L_FIT_MAX = 2
L_FIT_UNLIM = 3

# different dominance types
D_HET_ADV_DEBOER = 5
D_HET_ADV_MAXPD = 6
D_HET_ADV_AVPD = 7
D_HET_ADV_TAKNEI = 8
D_SEQ_B_AVPD_NO_HA = 12
D_SEQ_B_MAXPD_NO_HA = 13
D_SEQ_B_AVPD_W_HA = 14
D_SEQ_B_MAXPD_W_HA = 15

# different mutation types
T_PT_MUT = 1
T_MIC_CONV_MUT = 2
T_NEW_SEQ_MUT = 3

# different allele sort modes
M_NONE = 0
M_MOST_COMMON_ALLELE = 1
M_FITTEST_ALLELE = 2
M_LOWEST_ID = 3
M_HIGHEST_ID = 4

# different output levels
O_LEVEL_LOW = 1
O_LEVEL_STD = 2
O_LEVEL_MED = 3
O_LEVEL_HIGH = 4    # for debugging only - don't use for simulation runs

# different numbers of age classes
N_AGE_CL_10 = 10
N_AGE_CL_25 = 25
N_AGE_CL_100 = 100
N_AGE_CL_250 = 250
N_AGE_CL_1000 = 1000

# different overdominance (het. adv.) modes
HAM_CONST_DELTA = 1
HAM_BAD_ALL_DELTA = 2
HAM_POSDEP_DELTA = 3

# ### other constants ###
ROUND_PREC = 15
ERR_CODE = -9999
DICT_AMINO_ACIDS = {1: ('Alanine', 'Ala', 'A'),
                    2: ('Arginine', 'Arg', 'R'),
                    3: ('Asparagine', 'Asn', 'N'),
                    4: ('Aspartic acid', 'Asp', 'D'),
                    5: ('Cysteine', 'Cys', 'C'),
                    6: ('Glutamic acid', 'Glu', 'E'),
                    7: ('Glutamine', 'Gln', 'Q'),
                    8: ('Glycine', 'Gly', 'G'),
                    9: ('Histidine', 'His', 'H'),
                    10: ('Isoleucine', 'Ile', 'I'),
                    11: ('Leucine', 'Leu', 'L'),
                    12: ('Lysine', 'Lys', 'K'),
                    13: ('Methionine', 'Met', 'M'),
                    14: ('Phenylalanine', 'Phe', 'F'),
                    15: ('Proline', 'Pro', 'P'),
                    16: ('Serine', 'Ser', 'S'),
                    17: ('Threonine', 'Thr', 'T'),
                    18: ('Tryptophan', 'Trp', 'W'),
                    19: ('Tyrosine', 'Tyr', 'Y'),
                    20: ('Valine', 'Val', 'V')}
