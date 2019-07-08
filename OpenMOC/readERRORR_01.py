import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh, cholesky
from scipy.stats import norm

# --------------------------------------------------------------------------------------------------
# NJOY ERROR read demo
# --------------------------------------------------------------------------------------------------
# it is assumed that there is only a single temperature
# and a single background cross section
# and a single Legendre coefficient
# data generation from the feeding GROUPR module
# --------------------------------------------------------------------------------------------------
# according to the manual:
# -------------------------------------------------------------------------------------------------

myfile, NMAT ='u235_1100_groupr_error', 9228    # 9228
#myfile, NMAT='u238_1100_groupr_error', 9237   # 9237
#myfile, NMAT='o16_1100_groupr_errorr', 825   # 825
#myfile, NMAT='o16_600_groupr_error', 825     # 825
#myfile, NMAT='h01_600_groupr_error', 125     # 125


# read block size information
def readBlockSize(line):
    x0 = line[22:32+1]
    x1 = line[44:54+1]
    return int(x0), int(x1)

# read secondary block information
def readListSize(line):
    x0 = line[22:32+1]
    x1 = line[33:43+1]
    x2 = line[44:54+1]
    x3 = line[55:65+1]
    return int(x0), int(x1), int(x2), int(x3)

# read number of blocks
def readBlockNumbers(line):
    x0 = line[55:65+1]
    return int(x0)

# read number of errorr block entries
def readXSNumbers(line):
    x0 = line[45:55+1]
    return int(x0)

# regular data line
def readBlockData(line,n):
    x = []
    for i in range(n):
        r = line[0+i*11+1:10+i*11+1]
        r = r.replace('+', 'E+')
        r = r.replace('-', 'E-')
        x.append(float(r))
    return x

# read cross section from ERRORR format file
def ReadErrorR(myfile, mat, mf, mt):
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    Nentri = 0
    Nblocks = 0

    with open(myfile) as f:
        #x=f.readlines() # read all lines at once
        sum = 0
        cnt = 0
        block = []
        for line in f:
            #read_data = f.readline() # read a single line
            mat= line[66:69+1]
            mf = line[70:71 + 1]
            mt = line[72:76 + 1]
            #print(sum, '::', mat, ' ', mf, ' ', mt)
            if (int(mat)==mat0) and (int(mf)==mf0) and (int(mt)==mt0):
                #print(sum, '::', mat,' ',mf,' ',mt)
                #print(sum,' ** ',line)
                # mf 451 record -----------------------------------------------------
                if mt0==451:
                    if sum==1:
                        # read block size info
                        Nblocks, Nentri = readBlockSize(line)
                        #print(Nblocks, ' ', Nentri)
                        cnt=Nentri
                    else:
                        # read block data
                        if cnt>6:
                            x = readBlockData(line,6)
                            cnt = cnt - 6
                        else:
                            x = readBlockData(line,cnt)
                        block = block + x
                # mf 3 record -------------------------------------------------------
                elif mf0==3 and mt0>=1:
                    if sum==0:
                        # read block numbers
                        Nblocks = readXSNumbers(line)
                        block = [0] * Nblocks
                        #print(Nblocks)
                    else:
                        # pass
                        if cnt<Nblocks:
                            #print(' (2) ', line)
                            if Nblocks - cnt > 6:
                                x = readBlockData(line, 6)
                                block[(cnt):(cnt)+len(x)] = x
                                cnt = cnt + 6
                            else:
                                x = readBlockData(line, Nblocks - cnt)
                                block[(cnt):(cnt) + len(x)] = x
                                cnt = cnt + len(x)
                        else:
                            pass
                # mf (6,2) matrix record -------------------------------------------------------
                else:
                    pass
                sum +=1


    f.closed
    return Nblocks, block

# read covariance reaction pairs from ERRORR output file
def ReadCovPair(myfile, mat, mf):
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    Nentri = 0
    Nblocks = 0
    blocks = 0

    with open(myfile) as f:
        #x=f.readlines() # read all lines at once
        sum = 0
        cnt = 0
        block = []
        for line in f:
            #read_data = f.readline() # read a single line
            mat= line[66:69+1]
            mf = line[70:71 + 1]
            try:
                mx0 = int(line[23 - 1:33 - 0])
            except:
                mx0 = -1
            try:
                mz0 = int(line[34 - 1 :44 - 0])
            except:
                mz0 = -1
            try:
                mx1 = int(line[45 - 1 :55 - 0])
            except:
                mx1 = -1
            mt0 = line[72:76 + 0]
            #if (int(mat) == mat0) and (int(mf) == mf0):
            #    print(mx0," ", mz0, " ", mx1)
            if (int(mat)==mat0) and (int(mf)==mf0) and (mx0==0) and (mx1==0) and (mz0>0):
                #print(sum, '::', mat,' ',mf,' ',mt)
                #print(sum,' ** ',line)
                #print('reaction cov pair :',mt0,' -> ', mz0)
                block = block + [int(mt0), mz0]
                sum +=1

        Nblocks = sum
    f.closed
    return Nblocks, block

# read specific covariance matrix from ERRORR output file
# mt1: target reaction
# mt0: base reaction
def ReadCovMatrix(myfile, mat, mt1, mt0):
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    block = 0
    grp = 0

    with open(myfile) as f:
        #x=f.readlines() # read all lines at once
        sum = 0
        cnt = 0
        block = []
        ng2=0
        ng1=0
        for line in f:
            #read_data = f.readline() # read a single line
            mat= line[66:69+1]
            mf = line[70:71 + 1]
            try:
                mx0 = int(line[23 - 1:33 - 0])
            except:
                mx0 = -1
            try:
                mz0 = int(line[34 - 1 :44 - 0])
            except:
                mz0 = -1
            try:
                mx1 = int(line[45 - 1 :55 - 0])
            except:
                mx1 = -1
            m0  = line[72:76 + 0]

            #if (int(mat) == mat0) and (int(mf) == mf0):
            #    print(mx0," ", mz0, " ", mx1)
            if (int(mat)==mat0) and (int(mf)==mf0) and (mx0==0) and (mx1==0) and (mz0==mt1) and (int(m0)==mt0):
                #print(sum, '::', mat,' ',mf,' ',mt)
                #print(sum,' ** ',line)
                grp = int(line[56 - 1:66 + 0])
                print('reaction cov pair :',mt0,' -> ', mt1, ' with ', grp, ' groups.')
                #block = block + [int(mt0), mz0]
                block=np.zeros((grp,grp))
                sum +=1
            elif sum==1 and cnt<grp:
                #print(line)
                ng2, ng1, _, cnt = readListSize(line)
                sum +=1
                ng0=0
            elif sum==2:
                #print(line)
                if ng0 < ng2:
                    if ng2 - ng0 > 6:
                        x = np.asarray(readBlockData(line, 6))
                        block[cnt-1,(ng0+ng1-1):(ng0+ng1-1) + len(x)] = x
                        ng0 = ng0 + 6
                    else:
                        x = np.asarray(readBlockData(line, ng2 - ng0))
                        block[cnt-1,(ng0+ng1-1):(ng0+ng1-1) + len(x)] = x
                        ng0 = ng0 + len(x)
                if ng0 == ng2:
                    sum=1
            else:
                pass




        Nblocks = sum
    f.closed
    return block




# initialize energy group boundaries ############################
mat0, mf0, mt0 =NMAT, 1, 451
_ , block = ReadErrorR(myfile, mat0, mf0, mt0)

energyGroups=[]
if mt0==451:
    energyGroups = np.asarray(block)
print(energyGroups)
# ###############################################################


print('Total XS: ------------------')
# get total cross section #######################################
mat0, mf0, mt0 = NMAT, 3, 1
Nblocks , block = ReadErrorR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==1:
    totalXS = np.asarray(block)
print(totalXS)
# ###############################################################


print('Elastic XS: ------------------')
# get elastic cross section #####################################
mat0, mf0, mt0 = NMAT, 3, 2
Nblocks , block = ReadErrorR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==2:
    elasticXS = np.asarray(block[:])
print(elasticXS)
# ###############################################################

print('Fission XS: ------------------')
# get fission cross section #####################################
mat0, mf0, mt0 = NMAT, 3, 18
Nblocks , block = ReadErrorR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==18:
    fissionXS = np.asarray(block[:])
print(fissionXS)
# ###############################################################

print('Covariance reaction pairs: ------------------')
# get pairs of covariance data ##################################
mat0, mf0 = NMAT, 33
Nblocks , block = ReadCovPair(myfile, mat0, mf0)
pairs = np.asarray(block[:])
pairs = pairs.reshape(Nblocks, -1)
print(pairs)
# ###############################################################


print('Covariance matrix: ------------------')
# get specific cov matrix ##################################
mat0, mf0, mt1, mt0 = NMAT, 33, 1 , 1
# mt0: base reaction
# mt1: target reaction
covXS1_1 = ReadCovMatrix(myfile, mat0, mt1, mt0)
#print(covXS1_1)
# ###############################################################

# get specific cov matrix ##################################
mat0, mf0, mt1, mt0 = NMAT, 33, 2 , 2
covXS2_2 = ReadCovMatrix(myfile, mat0, mt1, mt0)
# ###############################################################

# get specific cov matrix ##################################
mat0, mf0, mt1, mt0 = NMAT, 33, 18 , 18
covXS18_18 = ReadCovMatrix(myfile, mat0, mt1, mt0)
# ###############################################################
# Colorplot of 2D array matplotlib
# https://stackoverflow.com/questions/16492830/colorplot-of-2d-array-matplotlib
#

x=plt.subplot(231)
x.plot(totalXS)
x.set_title('total XS')

x=plt.subplot(232)
x.plot(elasticXS)
x.set_title('elastic XS')

x=plt.subplot(233)
x.plot(fissionXS)
x.set_title('fission XS')

x=plt.subplot(234)
x.imshow(np.log(covXS1_1))
x.set_title('cov 1->1')

x=plt.subplot(235)
x.imshow(np.log(covXS2_2))
x.set_title('cov 2->2')

x=plt.subplot(236)
x.imshow(np.log(covXS18_18))
x.set_title('cov 18->18')

plt.tight_layout()
plt.show()

# https://scipy-cookbook.readthedocs.io/items/CorrelatedRandomSamples.html
# other resources:
# https://math.stackexchange.com/questions/268298/sampling-from-a-2d-normal-with-a-given-covariance-matrix
# https://stats.stackexchange.com/questions/32169/how-can-i-generate-data-with-a-prespecified-correlation-matrix
# https://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix
# https://en.wikipedia.org/wiki/Definiteness_of_a_matrix
#

# ########################################################
# we now generate sets of cross sections for total, elastic and fission
# ########################################################
num_samples = 100
num_variabl = 69  # number of energy groups

# Generate samples from three independent normally distributed random
# variables (with mean 0 and std. dev. 1).
x = norm.rvs(size=(num_variabl, num_samples))

# cholesky method
#c1_1 = cholesky(covXS1_1, lower=True)
# eigenvector method
evals, evecs = eigh(covXS1_1)
# clip small, negative eigenvalues
evals = evals.clip(0)
c1_1 = np.dot(evecs, np.diag(np.sqrt(evals)))

y1_1 = np.dot(c1_1, x) # we now have 'num_samples' column vectors of length 'num_variabl'
# we can calculate the covariance matrix from y and check against requirement
#re = np.cov(y)
#print(np.cov(y))

evals, evecs = eigh(covXS2_2)
evals = evals.clip(0)
c2_2 = np.dot(evecs, np.diag(np.sqrt(evals)))
y2_2 = np.dot(c2_2, x)

evals, evecs = eigh(covXS18_18)
evals = evals.clip(0)
c18_18 = np.dot(evecs, np.diag(np.sqrt(evals)))
y18_18 = np.dot(c18_18, x)

# ########################################################
# generate take the relative cov data and generate absolute XS values
# ########################################################
#yr1_1 = y1_1 * totalXS.reshape(num_variabl,-1)+ totalXS.reshape(num_variabl,-1)
#yr2_2 = y1_1 * elasticXS.reshape(num_variabl,-1)+ elasticXS.reshape(num_variabl,-1)
#yr18_18 = y1_1 * fissionXS.reshape(num_variabl,-1)+ fissionXS.reshape(num_variabl,-1)

yr1_1 = y1_1 + totalXS.reshape(num_variabl,-1)
yr2_2 = y2_2 + elasticXS.reshape(num_variabl,-1)
yr18_18 = y18_18 + fissionXS.reshape(num_variabl,-1)


x=plt.subplot(131)
#x.plot(np.log(yr1_1))
x.plot(yr1_1)
x.set_title('total XS')

x=plt.subplot(132)
x.plot(yr2_2)
x.set_title('elastic XS')

x=plt.subplot(133)
x.plot(yr18_18)
x.set_title('fission XS')

plt.tight_layout()
plt.show()

# ########################################################
# save total, elastic and fission cross section ensemble
# ########################################################
#np.savez('o16_1100_numpy_ensemble', yr1_1, yr2_2, yr18_18)
#np.savez('h01_600_numpy_ensemble', yr1_1, yr2_2)

#daten=np.load('h01_600_numpy_ensemble.npz')
#print(daten.files)