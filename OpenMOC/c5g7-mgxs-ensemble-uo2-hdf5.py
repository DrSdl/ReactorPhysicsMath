import h5py
import numpy
from sklearn.preprocessing import normalize

###############################################################################
# This file writes all of the materials data (multi-group nuclear
# cross-sections) for the OECD's C5G7 deterministic neutron transport
# benchmark problem to an HDF5 file. The script uses the h5py Python package
# to interact with the HDF5 file format. This may be a good example for those
# wishing ot write their nuclear data to an HDF5 file to import using the
# OpenMOC 'materialize' Python module.
###############################################################################


# Create the file to store C5G7 multi-groups cross-sections
f = h5py.File('c5g7-mgxs-ensemble-uo2.h5')
f.attrs["# groups"] = 69

# Create a group to specify that MGXS are split by material (vs. cell)
material_group = f.create_group('material')

############################################################################
# create UO2 data which mimicks convariance data from an ENDF file
# Later we use real data from ENDF, here we make just random variations
# to test the prcodure of random sampling

N = 100  # number of samples; MUST be equal with size of daten_h01_600 et al.

############################################################################
# load microscopic cross sections generated from NJOY and GROUPR and ERRORR
#
daten_h01_600  =numpy.load('h01_600_numpy_ensemble.npz')   # total, elastic XS ensemble yr1_1, yr2_2
daten_o16_600  =numpy.load('o16_600_numpy_ensemble.npz')   # total, elastic XS ensemble yr1_1, yr2_2
daten_o16_1100 =numpy.load('o16_1100_numpy_ensemble.npz')  # total, elastic XS ensemble yr1_1, yr2_2
daten_u235_1100=numpy.load('u235_1100_numpy_ensemble.npz') # total, elastic and fission XS ensemble yr1_1, yr2_2, yr18_18
daten_u238_1100=numpy.load('u238_1100_numpy_ensemble.npz') # total, elastic and fission XS ensemble yr1_1, yr2_2, yr18_18

daten_u235_1100_err=numpy.load('u235_1100_numpy_tables.npz') # energyGroups, totalXS, elasticXS, fissionXS, elasticXSM, nubarXS, chiXS
daten_u238_1100_err=numpy.load('u238_1100_numpy_tables.npz') # energyGroups, totalXS, elasticXS, fissionXS, elasticXSM, nubarXS, chiXS
daten_o16_1100_err =numpy.load('o16_1100_numpy_tables.npz')  # energyGroups, totalXS, elasticXS, elasticXSM
daten_o16_600_err  =numpy.load('o16_600_numpy_tables.npz')   # energyGroups, totalXS, elasticXS, elasticXSM
daten_h01_600_err  =numpy.load('h01_600_numpy_tables.npz')   # energyGroups, totalXS, elasticXS, elasticXSM


h01_600_t_xs = daten_h01_600['arr_0']     # total cross section
h01_600_e_xs = daten_h01_600['arr_1']     # elastic cross section
h01_600_e_xsm= daten_h01_600_err['arr_3'] # elastic matrix

o16_600_t_xs = daten_o16_600['arr_0']     # total cross section
o16_600_e_xs = daten_o16_600['arr_1']     # elastic cross section
o16_600_e_xsm= daten_o16_600_err['arr_3'] # elastic matrix

o16_1100_t_xs = daten_o16_1100['arr_0']     # total cross section
o16_1100_e_xs = daten_o16_1100['arr_1']     # elastic cross section
o16_1100_e_xsm= daten_o16_1100_err['arr_3'] # elastic matrix

u235_1100_t_xs = daten_u235_1100['arr_0']  # total cross section
u235_1100_e_xs = daten_u235_1100['arr_1']  # elastic cross section
u235_1100_f_xs = daten_u235_1100['arr_2']  # fission cross section
u235_1100_e_xsm   = daten_u235_1100_err['arr_4']  # elastic matrix
u235_1100_f_nubar = daten_u235_1100_err['arr_5']  # nubar
u235_1100_f_chi   = daten_u235_1100_err['arr_6']  # chi

u238_1100_t_xs = daten_u238_1100['arr_0']  # total cross section
u238_1100_e_xs = daten_u238_1100['arr_1']  # elastic cross section
u238_1100_f_xs = daten_u238_1100['arr_2']  # fission cross section
u238_1100_e_xsm   = daten_u238_1100_err['arr_4']  # elastic matrix
u238_1100_f_nubar = daten_u238_1100_err['arr_5']  # nubar
u238_1100_f_chi   = daten_u238_1100_err['arr_6']  # chi

density_h20 = 0.662 # g/cm3
density_uo2 = 10.15 # g/cm3
enr = 0.04 # 4% enrichment U235

uo2_m =  270.072  # g/mol molecular weight
h20_m =  18.015   # g/mol molecular weight
mol   = 6.022E+23  # 1 mol of particles
barn  = 1.0E-24    # cm2

nd_h20 = density_h20/h20_m  # mol/cm3
nd_uo2 = density_uo2/uo2_m  # mol/cm3

n_u235_fuel = nd_uo2*mol*enr         # atoms/cm3
n_u238_fuel = nd_uo2*mol*(1.0-enr)   # atoms/cm3
n_o_fuel    = nd_uo2*mol*2.0         # atoms/cm3
n_o_mod     = nd_h20*mol             # atoms/cm3
n_h_mod     = nd_h20*mol*2.0         # atoms/cm3

# microscopic cross sections are in barn
# macroscopic cross sections are in 1/cm

xs_t_moderator = (n_h_mod * h01_600_t_xs + n_o_mod * o16_600_t_xs)*barn
xs_e_moderator = (n_h_mod * h01_600_e_xs + n_o_mod * o16_600_e_xs)*barn

# bring energy groups into right order, i.e. from high to low energy
xs_t_moderator = numpy.flip(xs_t_moderator, axis=0)
xs_e_moderator = numpy.flip(xs_e_moderator, axis=0)

xs_t_fuel = (n_o_fuel * o16_1100_t_xs + n_u238_fuel * u238_1100_t_xs + n_u235_fuel * u235_1100_t_xs)*barn
xs_e_fuel = (n_o_fuel * o16_1100_e_xs + n_u238_fuel * u238_1100_e_xs + n_u235_fuel * u235_1100_e_xs)*barn
xs_f_fuel = (                           n_u238_fuel * u238_1100_f_xs + n_u235_fuel * u235_1100_f_xs)*barn

xs_t_fuel = numpy.flip(xs_t_fuel, axis=0)
xs_e_fuel = numpy.flip(xs_e_fuel, axis=0)
xs_f_fuel = numpy.flip(xs_f_fuel, axis=0)

chi_fuel   = (n_u235_fuel * u235_1100_f_chi + n_u238_fuel * u238_1100_f_chi)/(n_u235_fuel + n_u238_fuel)
nubar_fuel = (n_u235_fuel * u235_1100_f_nubar + n_u238_fuel * u238_1100_f_nubar)/(n_u235_fuel + n_u238_fuel)

chi_fuel   = numpy.flip(chi_fuel, axis=0)
nubar_fuel = numpy.flip(nubar_fuel, axis=0)
nubar_fuel = nubar_fuel[:,1]

chi_fuel   = chi_fuel.flatten()
nubar_fuel = nubar_fuel.flatten()

xsm_e_fuel = (n_o_fuel * o16_1100_e_xsm + n_u238_fuel * u238_1100_e_xsm + n_u235_fuel * u235_1100_e_xsm)*barn
xsm_e_moderator = (n_h_mod * h01_600_e_xsm + n_o_mod * o16_600_e_xsm)*barn

xsm_e_fuel = numpy.flip(xsm_e_fuel,axis=0)
xsm_e_fuel = numpy.flip(xsm_e_fuel,axis=1)

xsm_e_moderator = numpy.flip(xsm_e_moderator,axis=0)
xsm_e_moderator = numpy.flip(xsm_e_moderator,axis=1)

#print(xs_t_fuel[:,0])
#print(xs_f_fuel[:,0])
#print(xsm_e_fuel)
#print(xsm_e_moderator)
#stop()
###############################################################################
################################      UO2      ################################
###############################################################################

for n in range(100):
    #
    ###############################################################################
    ##############################      UO2             ###########################
    ###############################################################################
    # Create a subgroup for UO2 materials data
    uo2_name = 'UO2-'+str(n)
    uo2 = material_group.create_group(uo2_name)

    # Total cross section -----------------------------------------------------------
    # -------------------------------------------------------------------------------
    sigma_t = xs_t_fuel[:,n]

    # Scattering cross section -------------------------------------------------------
    # -------------------------------------------------------------------------------
    # https://stackoverflow.com/questions/8904694/how-to-normalize-a-2-dimensional-numpy-array-in-python-less-verbose
    sigma_s = normalize(xsm_e_fuel, axis=1, norm='l1')
    sigma_s = sigma_s * xs_e_fuel[:,n]

    # Fission cross section ---------------------------------------------------------
    # -------------------------------------------------------------------------------
    sigma_f = xs_f_fuel[:,n]

    # Nu fission cross section ------------------------------------------------------
    # -------------------------------------------------------------------------------
    nu_sigma_f = sigma_f * nubar_fuel
    chi = chi_fuel

    # Create datasets for each cross-section type
    uo2.create_dataset('total', data=sigma_t)
    uo2.create_dataset('scatter matrix', data=sigma_s)
    uo2.create_dataset('fission', data=sigma_f)
    uo2.create_dataset('nu-fission', data=nu_sigma_f)
    uo2.create_dataset('chi', data=chi)

    ###############################################################################
    ##############################      Guide Tube      ###########################
    ###############################################################################

    # Create a subgroup for guide tube materials data
    gd_name = 'Guide_Tube-' + str(n)
    guide_tube = material_group.create_group(gd_name)

    sigma_t = xs_t_moderator[:,n]

    sigma_s = normalize(xsm_e_moderator, axis=1, norm='l1')
    sigma_s = sigma_s * xs_e_moderator[:,n]

    sigma_f = numpy.zeros(69)
    nu_sigma_f = numpy.zeros(69)
    chi = numpy.zeros(69)

    # Create datasets for each cross-section type
    guide_tube.create_dataset('total', data=sigma_t)
    guide_tube.create_dataset('scatter matrix', data=sigma_s)
    guide_tube.create_dataset('fission', data=sigma_f)
    guide_tube.create_dataset('nu-fission', data=nu_sigma_f)
    guide_tube.create_dataset('chi', data=chi)


    ###############################################################################
    ################################      Water      ##############################
    ###############################################################################

    # Create a subgroup for water materials data
    h20_name = 'Water-' + str(n)
    water = material_group.create_group(h20_name)

    sigma_t = xs_t_moderator[:, n]

    sigma_s = normalize(xsm_e_moderator, axis=1, norm='l1')
    sigma_s = sigma_s * xs_e_moderator[:, n]

    sigma_f = numpy.zeros(69)
    nu_sigma_f = numpy.zeros(69)
    chi = numpy.zeros(69)

    # Create datasets for each cross-section type
    water.create_dataset('total', data=sigma_t)
    water.create_dataset('scatter matrix', data=sigma_s)
    water.create_dataset('fission', data=sigma_f)
    water.create_dataset('nu-fission', data=nu_sigma_f)
    water.create_dataset('chi', data=chi)



# Close the hdf5 data file
f.close()
