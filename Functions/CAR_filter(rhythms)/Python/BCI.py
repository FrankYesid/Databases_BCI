################### Functions  ###################
## 1. CAR
## 2. Filter
################ Database De Bci #################

###############   Function CAR   #################
# Inputs:
# SubjectsCAR = (subjects)
# subjects = X["X"]: Data to analysis
# --------------------------------------------------
# Outputs:
# Car = subjects[0,0][0,0] - mean(subjects[0,0][0,0],1)

# ############### FUNCTION rhythms #################
#  -------------------------> based on fncfiltband
#  Inputs:
# Xf = A_filterButter(Xcell,fs)
# X{trial}(samples x chan): Data to analysis

# fs: Sample frequency
# n : filter order
# freq : vector con las bandas
# --------------------------------------------------
# Outputs:
# Xf: Filtered data 
# Xf{trial}(samples x chan)

##################################################
# Copyright (C) 2018 Signal Processing and Recognition Group
# F.Y. Zapata Castano
##################################################

##################################################

#################### library #####################
import numpy as np
import numpy.matlib as npmat
import scipy.io as sio
from scipy import signal
import scipy.linalg
import argparse
import time 
import multiprocessing
import sys
###################################################
def CAR(subjects):
    SubjectsCAR = subjects
    fil = len(subjects)
    for l in range(0,fil): ## 9
        subject = subjects[l,0]
        fil = len(subject)
        SubjectCAR = subject         
        for k in range(0, fil): ## 273
            Sample = subject[k,0]
            fil,col = Sample.shape
            SampleCAR = Sample
            prom = np.mean(Sample,1)
            promedio = np.array(prom)
            prom_one = npmat.repmat(promedio,col,1)
            promedio = prom_one.T
            SampleCAR = Sample - promedio 
            SubjectCAR[k,0] = SampleCAR ## Subject
        SubjectsCAR[l,0] = SubjectCAR
    return SubjectsCAR
###################################################
################################################### 
def filter1(Freq,n,subjects,fs): ## ry
    ## fcn_filter, default=5 if n<5
    if n < 4:
        n = 5
    ## Filter_desing
    wn_freq = [(Freq[0]/(0.5*fs)), (Freq[1]/(0.5*fs))]
    ## [b_freq,a_freq] = butter(n,wn_freq)
    (b,a) = signal.butter(n,wn_freq,'bandpass')
    tr = len(subjects)
    Xfreq = subjects*0
    for k in range(0,tr):
        (s,c) = subjects[k,0].shape
        if c>s:
            subjects[k,0] = subjects[0,k].T
            [s,c] = subjects[k,0].shape
        for ch in range(0,c):
            x = subjects[k,0][:,ch]
            Xfreq[k,0][:,ch] = signal.filtfilt(b,a,x)
    return Xfreq
###################################################
###################### Main #######################
if __name__ == '__main__':
    ## Entradas de variables
    parser = argparse.ArgumentParser(description = "Program")
    parser.add_argument("Arch", help = "The potencials")
    parser.add_argument("car", help = "Common average reference (CAR)", nargs = "?")
    parser.add_argument("f", help = "Filter", nargs = "?")    
    parser.add_argument("f1", help = "frecuency 1 range", nargs = "?")
    parser.add_argument("f2",  help = "frecuency 2 range", nargs = "?")
    parser.add_argument("N",  help = "filter order", nargs = "?")
    args = parser.parse_args()
    Database = args.Arch              ## Carga de base de datos.
    database = sio.loadmat(Database)  ## 
    temporal = database
    subjects = temporal["X"]
    sub = len(subjects)
    N=len(subjects[0,0][0,0])
    fs = database["fs"][0,0]
    ## Manejo de la entrada de funcion solicitada.
    program_name = sys.argv[0]
    arguments = sys.argv[1:]
    count = len(arguments)
    name = 'dat'
    for parametros in sys.argv:
        if parametros == 'car':
            SubjectsCAR = CAR(subjects)
            temporal["X"] = SubjectsCAR
            print('Done CAR!')
            name='BCI_CAR.mat'
        elif parametros == 'f':
            F1 = arguments[arguments.index("f")+1]
            F2 = arguments[arguments.index("f")+2]
            Freq = []
            Freq.append(float(F1))
            Freq.append(float(F2))
            Num = arguments[arguments.index("f")+3] 
            subjects = temporal["X"]
            Xtemp = subjects*0
            for k in range(0,sub):
                Xtemp[k,0] = filter1(Freq,int(Num),subjects[k,0],fs)    
            temporal["X"] = Xtemp
            if name == 'car':
                name = 'BCI_CAR_filter.mat'
            else:
                name = 'BCI_filter.mat'
            print('Done FILTER!')
    print 'Nombre del archivo es: %s' % (name)
    sio.savemat(name,temporal)