import h5py
import numpy
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt

Nmax = 67

#

fission = []
total   = []
phi     = []


for n in range(Nmax):
    fname = './210819/result_210819_'+str(n+1)+'.npz'
    fname = './result_210819_'+str(n+1)+'.npz'
    print(fname)
    #
    daten_01 = numpy.load(fname)
    d01_fi   = daten_01['arr_0'] # fission
    d01_to   = daten_01['arr_1'] # total XS
    d01_pi   = daten_01['arr_2'] # phi
    #
    #print(d01_fi.shape)
    #print(d01_to.shape)
    #print(d01_pi.shape)
    # ----------------------------------
    x = numpy.sum(d01_fi, axis = 1)
    x = numpy.sum(x     , axis = 1)
    if len(fission)==0:
        fission=numpy.copy(x)
    else:
        fission=numpy.vstack((fission,x))
    # ----------------------------------
    x = numpy.sum(d01_to, axis = 1)
    x = numpy.sum(x     , axis = 1)
    if len(total)==0:
        total=numpy.copy(x)
    else:
        total=numpy.vstack((total,x))
    # ----------------------------------
    x = numpy.sum(d01_pi, axis = 1)
    x = numpy.sum(x     , axis = 1)
    if len(phi)==0:
        phi=numpy.copy(x)
    else:
        phi=numpy.vstack((phi,x))

#print(fission)
#print(total)

#plt.hist(total, bins='auto')
#plt.show()

Ng=69
boxMe = []
# phi
# ----------------------------------------------
for n in range(Ng):
    boxMe.append(phi[:,n])

fig, ax = plt.subplots(figsize=(20,10))
ax.boxplot(boxMe)
#ax.set_ylim([0.0,0.05])
##ax.set_aspect(1.5)
plt.xlabel('group')
plt.ylabel('flux (rel. units)')
##ax.xaxis.label.set_size(0.5)
plt.xticks(fontsize=5)
plt.show()
# total
# ----------------------------------------------
boxMe = []

for n in range(Ng):
    boxMe.append(total[:,n])

fig, ax = plt.subplots(figsize=(20,10))
ax.boxplot(boxMe)
ax.set_ylim([0.0,1.5])
##ax.set_aspect(1.5)
plt.xlabel('group')
plt.ylabel('total XS (rel. units)')
##ax.xaxis.label.set_size(0.5)
plt.xticks(fontsize=5)
plt.show()
# fission
# ----------------------------------------------
boxMe = []

for n in range(Ng):
    boxMe.append(fission[:,n])

fig, ax = plt.subplots(figsize=(20,10))
ax.boxplot(boxMe)
ax.set_ylim([0.0,0.03])
##ax.set_aspect(1.5)
plt.xlabel('group')
plt.ylabel('fission XS (rel. units)')
##ax.xaxis.label.set_size(0.5)
plt.xticks(fontsize=5)
plt.show()


    


