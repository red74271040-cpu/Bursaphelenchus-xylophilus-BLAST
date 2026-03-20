import math
import warnings
from Bio import SeqUtils, Seq
from Bio import BiopythonWarning

# RNA/RNA Thermodynamic lookup tables
RNA_NN3 = {
    'init': (6.40, 6.99), 'init_A/T': (3.85, 11.04), 'init_G/C': (0, 0),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT': (-7.09, -19.8), 'AT/TA': (-9.11, -25.8), 'TA/AT': (-8.50, -22.9),
    'CA/GT': (-11.03, -28.8), 'GT/CA': (-11.98, -31.3),
    'CT/GA': (-10.90, -28.5), 'GA/CT': (-13.21, -34.9),
    'CG/GC': (-10.88, -27.4), 'GC/CG': (-16.04, -40.6),
    'GG/CC': (-14.18, -35.0), 'GT/TG': (-13.83, -46.9),
    'GG/TT': (-17.82, -56.7), 'AG/TT': (-3.96, -11.6),
    'TG/AT': (-0.96, -1.8), 'TT/AG': (-10.38, -31.8), 'TG/GT': (-12.64, -38.9),
    'AT/TG': (-7.39, -21.0), 'CG/GT': (-5.56, -13.9), 'CT/GG': (-9.44, -24.7),
    'GG/CT': (-7.03, -16.8), 'GT/CG': (-11.09, -28.8)}

# Terminal mismatch table
DNA_TMM1 = {'AA/TA': (-3.1, -7.8), 'TA/AA': (-2.5, -6.3), 'CA/GA': (-4.3, -10.7), 'GA/CA': (-8.0, -22.5), 'AC/TC': (-0.1, 0.5), 'TC/AC': (-0.7, -1.3), 'CC/GC': (-2.1, -5.1), 'GC/CC': (-3.9, -10.6), 'AG/TG': (-1.1, -2.1), 'TG/AG': (-1.1, -2.7), 'CG/GG': (-3.8, -9.5), 'GG/CG': (-0.7, -19.2), 'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9), 'CT/GT': (-6.1, -16.9), 'GT/CT': (-7.4, -21.2), 'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 'CA/GC': (-2.6, -5.9), 'CC/GA': (-2.7, -6.0), 'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1), 'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0), 'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 'CC/GT': (-3.2, -8.0), 'CT/GC': (-3.9, -10.6), 'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8), 'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2), 'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 'CA/GG': (-3.9, -9.6), 'CG/GA': (-6.0, -15.5), 'GA/CG': (-4.3, -11.1), 'GG/CA': (-4.6, -11.4), 'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8), 'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 'CG/GT': (-3.8, -9.0), 'CT/GG': (-6.6, -18.7), 'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1), 'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)}

# Dangling ends table (RNA)
RNA_DE2 = {'.T/AA': (-4.9, -13.2), '.T/CA': (-0.9, -1.3), '.T/GA': (-5.5, -15.1), '.T/TA': (-2.3, -5.5), '.G/AC': (-9.0, -23.5), '.G/CC': (-4.1, -10.6), '.G/GC': (-8.6, -22.2), '.G/TC': (-7.5, -20.31), '.C/AG': (-7.4, -20.3), '.C/CG': (-2.8, -7.7), '.C/GG': (-6.4, -16.4), '.C/TG': (-3.6, -9.7), '.T/AG': (-4.9, -13.2), '.T/CG': (-0.9, -1.3), '.T/GG': (-5.5, -15.1), '.T/TG': (-2.3, -5.5), '.A/AT': (-5.7, -16.1), '.A/CT': (-0.7, -1.9), '.A/GT': (-5.8, -16.4), '.A/TT': (-2.2, -6.8), '.G/AT': (-5.7, -16.1), '.G/CT': (-0.7, -1.9), '.G/GT': (-5.8, -16.4), '.G/TT': (-2.2, -6.8), 'AT/.A': (-0.5, -0.6), 'CT/.A': (6.9, 22.6), 'GT/.A': (0.6, 2.6), 'TT/.A': (0.6, 2.6), 'AG/.C': (-1.6, -4.5), 'CG/.C': (0.7, 3.2), 'GG/.C': (-4.6, -14.8), 'TG/.C': (-0.4, -1.3), 'AC/.G': (-2.4, -6.1), 'CC/.G': (3.3, 11.6), 'GC/.G': (0.8, 3.2), 'TC/.G': (-1.4, -4.2), 'AT/.G': (-0.5, -0.6), 'CT/.G': (6.9, 22.6), 'GT/.G': (0.6, 2.6), 'TT/.G': (0.6, 2.6), 'AA/.T': (1.6, 6.1), 'CA/.T': (2.2, 8.1), 'GA/.T': (0.7, 3.5), 'TA/.T': (3.1, 10.6), 'AG/.T': (1.6, 6.1), 'CG/.T': (2.2, 8.1), 'GG/.T': (0.7, 3.5), 'TG/.T': (3.1, 10.6)}

def _check(seq, method):
    seq = ''.join(seq.split()).upper()
    seq = str(Seq.Seq(seq).back_transcribe())
    return seq

def salt_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None):
    Mon = Na + K + Tris / 2.0
    mon = Mon * 1e-3
    if method == 5:
        return 0.368 * (len(seq) - 1) * math.log(mon)
    return 16.6 * math.log10(mon)

def calculate_free_energy(seq, check=True, strict=True, c_seq=None, shift=0, nn_table=RNA_NN3,
          tmm_table=DNA_TMM1, de_table=RNA_DE2,
          dnac1=25, dnac2=0, selfcomp=False, Na=20, K=50, saltcorr=5):
    
    seq = str(seq)
    if not c_seq:
        c_seq = str(Seq.Seq(seq).complement())
    
    tmpseq, tmp_cseq = seq, c_seq
    deltaH, deltaS = 0, 0
    dH, dS = 0, 1

    if shift or len(seq) != len(c_seq):
        if shift > 0: tmpseq = '.' * shift + seq
        if shift < 0: tmp_cseq = '.' * abs(shift) + c_seq
        while tmpseq.startswith('..') or tmp_cseq.startswith('..'):
            tmpseq, tmp_cseq = tmpseq[1:], tmp_cseq[1:]
        
        if tmpseq.startswith('.') or tmp_cseq.startswith('.'):
            left_de = tmpseq[:2] + '/' + tmp_cseq[:2]
            if left_de in de_table:
                deltaH += de_table[left_de][dH]
                deltaS += de_table[left_de][dS]
            tmpseq, tmp_cseq = tmpseq[1:], tmp_cseq[1:]

    # Initiation
    deltaH += nn_table['init'][dH]
    deltaS += nn_table['init'][dS]

    # Nearest Neighbor Calculation
    for i in range(len(tmpseq) - 1):
        neighbors = tmpseq[i:i+2] + '/' + tmp_cseq[i:i+2]
        if neighbors in nn_table:
            deltaH += nn_table[neighbors][dH]
            deltaS += nn_table[neighbors][dS]
        elif neighbors[::-1] in nn_table:
            deltaH += nn_table[neighbors[::-1]][dH]
            deltaS += nn_table[neighbors[::-1]][dS]

    R = 1.987
    if saltcorr == 5:
        deltaS += salt_correction(Na=Na, K=K, method=5, seq=seq)

    tao = 273.15 + 22
    deltaG = (deltaH * 1000 - tao * deltaS) / 1000
    return deltaG