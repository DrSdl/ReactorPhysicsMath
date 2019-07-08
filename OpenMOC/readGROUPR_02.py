import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------
# NJOY GROUPR read demo
# --------------------------------------------------------------------------------------------------
# it is assumed that there is only a single temperature
# and a single background cross section for the group
# and a single Legendre coefficient
# data generation GROUPR
# --------------------------------------------------------------------------------------------------
# according to the manual:
# The NJOY Nuclear Data Processing System, Version 2016 (LA-UR-17-20093, page 229)
#
# For simple cross section vectors (MF=3), NG2 is 2, and A contains the two Fortran arrays
# FLUX(NL,NZ), SIGMA(NL,NZ) in that order.
#
# For ratio quantities like fission nubar, NG2 is 3, and A contains
# FLUX(NL,NZ), RATIO(NL,NZ), SIGMA(NL,NZ).
#
# For transfer matrices (MF=6, 16, 21, etc.), A contains
# FLUX(NL,NZ), MATRIX(NL,NZ,NG2-1).
#
# For delayed neutron spectra (MF=5), NL is used to index the time groups, NZ is 1, and
# there is only one incident energy record (IG=IGN). The array A contains
# LAMBDA(NL), CHID(NL,NG2-1),
# where LAMBDA contains the delayed-neutron time constants and CHID contains the spectra.
# -------------------------------------------------------------------------------------------------

myfile, NMAT ='u235_1100_groupr', 9228
#myfile, NMAT ='u238_1100_groupr', 9237
#myfile, NMAT ='o16_1100_groupr', 825
#myfile, NMAT ='o16_600_groupr',  825
#myfile, NMAT ='h01_600_groupr',  125

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

# read MF5 block header
def readBlockMF5(line):
    x0 = line[23-1:33]
    x1 = line[56-1:66]
    return int(x0), int(x1)

# regular data line
def readBlockData(line,n):
    x = []
    for i in range(n):
        r = line[0+i*11+1:10+i*11+1]
        r = r.replace('+', 'E+')
        r = r.replace('-', 'E-')
        x.append(float(r))
    return x

# read cross section from GROUPR format file
def ReadGroupR(myfile, mat, mf, mt):
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    Nentri = 0
    Nblocks = 0

    with open(myfile) as f:
        #x=f.readlines() # read all lines at once
        sum = 0
        cnt = 0
        block = []
        blocc = []
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
                    #print(mf0, ' ', sum)
                    if sum==0:
                        # read block numbers
                        Nblocks = readBlockNumbers(line)
                        block = [0] * Nblocks
                        #cnt=Nblocks
                        #print(line, ' ', Nblocks)
                    else:
                        # pass
                        if cnt==0:
                            #print('(1) ',line, " ", cnt)
                            ng2, ig2lo, nw, ig = readListSize(line)
                            cnt = nw
                        elif cnt>0:
                            #print(' (2) ', line)
                            x = readBlockData(line, cnt)
                            cnt = cnt - 6
                            if block[ig-1]==0:
                                block[ig - 1] = x
                            else:
                                block[ig - 1] = block[ig - 1] + x
                        else:
                            pass
                        if(cnt<0 and ig<Nblocks):
                            cnt = 0
                # mf (6,2) matrix record -------------------------------------------------------
                elif (mf0 == 6 and mt0 == 2) or (mf0 == 6 and mt0 == 221):
                    # print(mf0, ' ', sum)
                    if sum == 0:
                        # read block numbers
                        Nblocks = readBlockNumbers(line)
                        block = [[0] * Nblocks for i in range(Nblocks)]
                        #block = [[0] * Nblocks ] * Nblocks
                    else:
                        # pass
                        if cnt == 0:
                            # print('(1) ',line, " ", cnt)
                            ng2, ig2lo, nw, ig = readListSize(line)
                            cnt = nw
                        elif cnt > 6:
                            x = readBlockData(line, 6)
                            if (cnt==nw):
                                for j in range(2):
                                    block[ig - 1][ig2lo+j-1] = x[2+j*2]
                                ig2lo = ig2lo + 2
                            if (cnt<nw):
                                for j in range(3):
                                    block[ig - 1][ig2lo+j-1] = x[j*2]
                                ig2lo = ig2lo + 3
                            cnt = cnt - 6
                        elif cnt > 0 and cnt <=6:
                            x = readBlockData(line, cnt)
                            if (cnt == nw):
                                for j in range(int((nw-2)/2)):
                                    block[ig - 1][ig2lo + j-1] = x[2 + j * 2]
                                ig2lo = ig2lo + int((nw-2)/2)
                            if (cnt<nw):
                                for j in range(int((cnt-0)/2)):
                                    block[ig - 1][ig2lo + j-1] = x[j * 2]
                                ig2lo = ig2lo + int((cnt-0)/2)
                            cnt = cnt - 6
                            # print(" *** ", np.asarray(block))
                        else:
                            pass
                        if (cnt < 0 and ig < Nblocks):
                            cnt = 0
                # mf 5 record -------------------------------------------------------
                elif mf0==5 and mt0==18:
                    #print(mf0, ' ', sum)
                    if sum==0:
                        # read block numbers
                        groups, Nblocks = readBlockMF5(line)
                        block = np.zeros(Nblocks)
                        #print(groups, Nblocks)
                    else:
                        # pass
                        if cnt==0:
                            #print('(1) ',line, " ", cnt)
                            ng2, ig2lo, nw, ig = readListSize(line)
                            cnt = nw
                        elif cnt>0:
                            #print(' (2) ', line)
                            if cnt >= 6:
                                x = readBlockData(line, 6)
                                block[(Nblocks-cnt):(Nblocks-cnt) + len(x)]=x
                                cnt = cnt - 6
                            else:
                                x = readBlockData(line, cnt)
                                block[(Nblocks - cnt):(Nblocks - cnt) + len(x)] = x
                                cnt = cnt - 6
                        else:
                            pass
                        if(cnt<0 and ig<Nblocks):
                            cnt = 0
                else:
                    pass
                sum +=1


    f.closed
    return Nblocks, block

print('thermal scattering matrix: ------------------')
# get thermal scattering matrix ################################
mat0, mf0, mt0 =NMAT, 6, 221
#mat0, mf0, mt0 =825, 6, 221
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==6 and mt0==221:
    thermalXSM = np.asarray(block[::])
print(thermalXSM)
print(thermalXSM[0:5,0:5])
# ###############################################################


print('elastic scattering matrix: ------------------')
# get elastic scattering matrix #################################
mat0, mf0, mt0 =NMAT, 6, 2
#mat0, mf0, mt0 =825, 6, 2
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==6 and mt0==2:
    elasticXSM = np.asarray(block[::])
print(elasticXSM)
print(elasticXSM[0:5,0:5])
# ###############################################################
# remember:
# np.sum(thermalXSM[0])  == np.sum(elasticXSM[0])
# np.sum(thermalXSM[K])  == np.sum(elasticXSM[K])
# thermalXSM represents the broadened elasticXSM

# initialize energy group boundaries ############################
mat0, mf0, mt0 =NMAT, 1, 451
_ , block = ReadGroupR(myfile, mat0, mf0, mt0)

energyGroups=[]
if mt0==451:
    energy0 = block[0]
    energy1 = block[len(block)-2]
    energyGroups=np.asarray( [energy0] + block[2:len(block)-1] )
print(energyGroups)
# ###############################################################

print('Total XS: ------------------')
# get total cross section #######################################
mat0, mf0, mt0 =NMAT, 3, 1
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==1:
    totalXS = np.asarray(block[:])
    totalXS = totalXS.reshape(Nblocks,-1)
print(totalXS)
# ###############################################################

print('Elastic XS: ------------------')
# get elastic cross section #####################################
mat0, mf0, mt0 =NMAT, 3, 2
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==2:
    elasticXS = np.asarray(block[:])
    elasticXS = elasticXS.reshape(Nblocks,-1)
print(elasticXS)
# ###############################################################

print('Fission XS: ------------------')
# get fission cross section #####################################
mat0, mf0, mt0 =NMAT, 3, 18
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==18:
    fissionXS = np.asarray(block[:])
    fissionXS = fissionXS.reshape(Nblocks,-1)
print(fissionXS)
# ###############################################################


## get cross section ############################################
#mat0, mf0, mt0 =9228, 3, 51
#Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)
#
#if mf0==3 and mt0==51:
#    block = [ [0,0] if x==0 else x for x in block]
#    prodXS = np.asarray(block[:])
#    prodXS = prodXS.reshape(Nblocks,-1)
#print(prodXS)
# ###############################################################

print('nubar: ------------------')
# get nnubar total ##############################################
mat0, mf0, mt0 =NMAT, 3, 452
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==3 and mt0==452:
    block = [[0, 0] if x == 0 else x for x in block]
    nubarXS = np.asarray(block[:])
    nubarXS = nubarXS.reshape(Nblocks,-1)
print(nubarXS)
# ###############################################################

print('chi prompt: ------------------')
# get chi prompt ##############################################
mat0, mf0, mt0 =NMAT, 5, 18
Nblocks , block = ReadGroupR(myfile, mat0, mf0, mt0)

if mf0==5 and mt0==18:
    block = [[0, 0] if x == 0 else x for x in block]
    chiXS = np.asarray(block[:])
    chiXS = chiXS.reshape(Nblocks,-1)
print(chiXS)
# ###############################################################

np.savez('u235_1100_numpy_tables', energyGroups, totalXS, elasticXS, fissionXS, elasticXSM, nubarXS, chiXS)
#np.savez('u238_1100_numpy_tables', energyGroups, totalXS, elasticXS, fissionXS, elasticXSM, nubarXS, chiXS)
#np.savez('o16_1100_numpy_tables', energyGroups, totalXS, elasticXS, elasticXSM)
#np.savez('o16_600_numpy_tables', energyGroups, totalXS, elasticXS, elasticXSM)
#np.savez('h01_600_numpy_tables', energyGroups, totalXS, elasticXS, elasticXSM)
