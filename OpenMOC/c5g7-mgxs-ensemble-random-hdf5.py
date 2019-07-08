import h5py
import numpy


###############################################################################
# This file writes all of the materials data (multi-group nuclear
# cross-sections) for the OECD's C5G7 deterministic neutron transport
# benchmark problem to an HDF5 file. The script uses the h5py Python package
# to interact with the HDF5 file format. This may be a good example for those
# wishing ot write their nuclear data to an HDF5 file to import using the
# OpenMOC 'materialize' Python module.
###############################################################################


# Create the file to store C5G7 multi-groups cross-sections
f = h5py.File('c5g7-mgxs-ensemble.h5')
f.attrs["# groups"] = 7

# Create a group to specify that MGXS are split by material (vs. cell)
material_group = f.create_group('material')

############################################################################
# create UO2 data which mimicks convariance data from an ENDF file
# Later we use real data from ENDF, here we make just random variations
# to test the prcodure of random sampling

N = 100  # number of samples

###############################################################################
################################      UO2      ################################
###############################################################################
# cross sections are from:
# https://www.oecd-nea.org/science/docs/2005/nsc-doc2005-16.pdf
# https://www.oecd-nea.org/science/docs/2003/nsc-doc2003-16.pdf
# https://www.oecd-nea.org/science/docs/1994/nsc-doc94-28.pdf
# EPRI-CPM group structure
# http://serpent.vtt.fi/mediawiki/index.php/EPRI-CPM_69_group_structure
#
for n in range(100):
    # 
    # Create a subgroup for UO2 materials data
    uo2_name = 'UO2-'+str(n)
    uo2 = material_group.create_group(uo2_name)

    # Total cross section -----------------------------------------------------------
    # -------------------------------------------------------------------------------
    sigma_t = numpy.array([1.779490E-01, 3.298050E-01, 4.803880E-01,
                           5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01])
    # create random matrix according to above shape
    shp = sigma_t.shape
    dis = 0.1*numpy.random.rand(*shp)
    sigma_t = sigma_t + dis

    # Scattering cross section -------------------------------------------------------
    # -------------------------------------------------------------------------------
    sigma_s = numpy.array([1.275370E-01, 4.237800E-02, 9.437400E-06,
                           5.516300E-09, 0., 0., 0., 0., 3.244560E-01,
                           1.631400E-03, 3.142700E-09, 0., 0., 0., 0.,
                           0., 4.509400E-01, 2.679200E-03, 0., 0., 0.,
                           0., 0., 0., 4.525650E-01, 5.566400E-03, 0.,
                           0., 0., 0., 0., 1.252500E-04, 2.714010E-01,
                           1.025500E-02, 1.002100E-08, 0., 0., 0., 0.,
                           1.296800E-03, 2.658020E-01, 1.680900E-02,
                           0., 0., 0., 0., 0., 8.545800E-03, 2.730800E-01])
    # create random matrix according to above shape
    shp = sigma_s.shape
    dis = 0.01*numpy.random.rand(*shp)
    sigma_s = sigma_s + dis

    # Fission cross section ---------------------------------------------------------
    # -------------------------------------------------------------------------------
    sigma_f = numpy.array([7.212060E-03, 8.193010E-04, 6.453200E-03,
                           1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01])
    # create random matrix according to above shape
    shp = sigma_f.shape
    dis = 0.01*numpy.random.rand(*shp)
    sigma_f = sigma_f + dis

    # Nu fission cross section ------------------------------------------------------
    # -------------------------------------------------------------------------------
    nu_sigma_f = numpy.array([2.005998E-02, 2.027303E-03, 1.570599E-02,
                              4.518301E-02, 4.334208E-02, 2.020901E-01,
                              5.257105E-01])
    # create random matrix according to above shape
    shp = nu_sigma_f.shape
    dis = 0.01*numpy.random.rand(*shp)
    nu_sigma_f = nu_sigma_f + dis



    chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04,
                       1.17610E-07, 0., 0., 0.])

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
guide_tube = material_group.create_group('Guide Tube')

sigma_t = numpy.array([1.260320E-01, 2.931600E-01, 2.842400E-01,
                       2.809600E-01, 3.344400E-01, 5.656400E-01,
                       1.172150E+00])
sigma_s = numpy.array([6.616590E-02, 5.907000E-02, 2.833400E-04,
                       1.462200E-06, 2.064200E-08, 0., 0., 0., 2.403770E-01,
                       5.243500E-02, 2.499000E-04, 1.923900E-05,
                       2.987500E-06, 4.214000E-07, 0., 0., 1.832970E-01,
                       9.239700E-02, 6.944600E-03, 1.0803000E-03,
                       2.056700E-04, 0., 0., 0., 7.885110E-02, 1.701400E-01,
                       2.588100E-02, 4.929700E-03, 0., 0., 0., 3.733300E-05,
                       9.973720E-02, 2.067900E-01, 2.447800E-02, 0., 0., 0.,
                       0., 9.172600E-04, 3.167650E-01, 2.387700E-01, 0., 0.,
                       0., 0., 0., 4.979200E-02, 1.099120E+00])

sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

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
water = material_group.create_group('Water')

sigma_t = numpy.array([1.592060E-01, 4.129700E-01, 5.903100E-01,
                       5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00])
sigma_s = numpy.array([4.447770E-02, 1.134000E-01, 7.234700E-04,
                       3.749900E-06, 5.318400E-08, 0., 0., 0., 2.823340E-01,
                       1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06,
                       1.045500E-06, 0., 0., 3.452560E-01, 2.245700E-01,
                       1.699900E-02, 2.644300E-03, 5.034400E-04, 0., 0., 0.,
                       9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02,
                       0., 0., 0., 7.143700E-05, 1.391380E-01, 5.118200E-01,
                       6.122900E-02, 0., 0., 0., 0., 2.215700E-03, 6.999130E-01,
                       5.373200E-01, 0., 0., 0., 0., 0., 1.324400E-01,
                       2.480700E+00])
sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

# Create datasets for each cross-section type
water.create_dataset('total', data=sigma_t)
water.create_dataset('scatter matrix', data=sigma_s)
water.create_dataset('fission', data=sigma_f)
water.create_dataset('nu-fission', data=nu_sigma_f)
water.create_dataset('chi', data=chi)



# Close the hdf5 data file
f.close()
