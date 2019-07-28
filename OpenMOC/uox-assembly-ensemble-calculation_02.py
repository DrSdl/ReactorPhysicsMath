
# ******************************************************************************************

import matplotlib.pyplot as plt

import openmoc
import openmoc.plotter as plotter
import openmoc.process as process
from openmoc.materialize import load_from_hdf5

RunID = drsdlID


# ## Simulation Runtime Parameters

# In[2]:

num_threads = 4
azim_spacing = 0.05
num_azim = 16
tolerance = 1E-5
max_iters = 50


# ## Initialize Materials

# In[3]:

materials = load_from_hdf5(filename='c5g7-mgxs-ensemble-uo2.h5', directory='..')
print(materials.keys())


# ## Create Bounding Surfaces

# In[4]:

# Create ZCylinder for the fuel
fuel_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.54)



# ## Create Fuel Pins

# In[6]:

# Create ensemble of fuel pins
# 4.0%  UOX pin cell base case
Nmat = 100 # number of UO2 pin materials in the ensemble 
uox40_array = []
gd_array = []

for n in range(Nmat):
    uo2_name = 'UO2-' + str(n)
    uox40_cell = openmoc.Cell()
    uox40_cell.setFill(materials[uo2_name])
    uox40_cell.setNumRings(3)
    uox40_cell.setNumSectors(8)
    uox40_cell.addSurface(-1, fuel_radius)

    uox40 = openmoc.Universe(name=uo2_name)
    uox40.addCell(uox40_cell)
    uox40_array.append(uox40)

print(len(uox40_array))


# In[7]:

# Guide tube pin cell
for n in range(Nmat):
    gd_name = 'Guide_Tube-' + str(n)
    guide_tube_cell = openmoc.Cell()
    guide_tube_cell.setFill(materials[gd_name])
    guide_tube_cell.setNumRings(3)
    guide_tube_cell.setNumSectors(8)
    guide_tube_cell.addSurface(-1, fuel_radius)

    guide_tube = openmoc.Universe(name=gd_name)
    guide_tube.addCell(guide_tube_cell)
    gd_array.append(guide_tube)


# In[8]:

# Add moderator rings to each pin cell
sum=0
for uox40 in uox40_array:
    # Moderator definition
    moderator = openmoc.Cell()
    h20_name = 'Water-' + str(sum)
    moderator.setFill(materials[h20_name])
    moderator.addSurface(+1, fuel_radius)
    moderator.setNumRings(3)
    moderator.setNumSectors(8)
    uox40.addCell(moderator)
    sum +=1

sum=0    
for gd in gd_array:
    # Moderator definition
    moderator = openmoc.Cell()
    h20_name = 'Water-' + str(sum)
    moderator.setFill(materials[h20_name])
    moderator.addSurface(+1, fuel_radius)
    moderator.setNumRings(3)
    moderator.setNumSectors(8)
    gd.addCell(moderator)    
    sum +=1


# ## Create Fuel Assembly Class

# In[9]:

class FuelAssembly:
    def __init__(self):
        # define fuel assembly basic parameters
        # CellFills for the assembly
        self.assembly1_cell = openmoc.Cell(name='Fuel Assembly')
        self.assembly1 = openmoc.Universe(name='Fuel Assembly')
        self.assembly1.addCell(self.assembly1_cell)
        
        # A mixed enrichment PWR MOX fuel assembly
        self.assembly = openmoc.Lattice(name='UOX Assembly')
        self.assembly.setWidth(width_x=1.26, width_y=1.26)
        
        # Create a template to map to pin cell types
        self.template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 1, 1, 1],
                         [1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 4, 1, 1, 4, 1, 1, 1, 1, 1, 4, 1, 1, 4, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1],
                         [1, 1, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        self.root_cell=[]
        self.root_universe=[]
        
    def setUniverse(self, uox40, guide_tube):
        universes = {1 : uox40, 4 : guide_tube}
        # loop over all pins
        for i in range(17):
            for j in range(17):
                self.template[i][j] = universes[self.template[i][j]]
        # write pin layout into assembly structure        
        self.assembly.setUniverses([self.template])
    
    def setRootUniverse(self, boundary):
        # Root Cell
        self.root_cell = openmoc.Cell(name='Full Geometry')
        self.root_cell.setFill(self.assembly)
        self.root_cell.setRegion(boundary)
        # Root Universe
        self.root_universe = openmoc.Universe(name='Root Universe')
        self.root_universe.addCell(self.root_cell)
        
    def getRootUniverse(self):
        return self.root_universe


# Initialize FuelAssembly classes

fuelStack=[]

for n in range(Nmat):
    x = FuelAssembly()
    x.setUniverse(uox40_array[n],gd_array[n])

    # Create planes to bound the entire geometry
    boundary = openmoc.RectangularPrism(21.42, 21.42)
    boundary.setBoundaryType(openmoc.REFLECTIVE)

    x.setRootUniverse(boundary)
    fuelStack.append(x)


# ## Run the big XS generation loop

# In[ ]:

import numpy as np

# Aggregate the total cross sections for each Material
# into a dictionary to pass to the mesh tally
domains_to_coeffsF = {}
domains_to_coeffsT = {}
domains_to_coeffsSx7 = {}
domains_to_coeffsPhi = {}

# Create OpenMOC Mesh on which to tally fission rates
mesh = openmoc.process.Mesh()
mesh.dimension = [17*2, 17*2]
mesh.lower_left = [-21.42/2, -21.42/2]
mesh.upper_right= [ 21.42/2,  21.42/2]
mesh.width = [21.42/(17*2), 21.42/(17*2)]

# Array to store macroscopic XS results
Fission_results=[]
TotalXS_results=[]
Phi_results=[]

for n in range(RunID,RunID+1):
    # Initialize CMFD
    cmfd = openmoc.Cmfd()
    cmfd.setSORRelaxationFactor(1.5)
    cmfd.setLatticeStructure(17,17)
    #cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
    cmfd.setKNearest(3)

    # Initialize Geometry¶
    geometry = openmoc.Geometry()
    geometry.setRootUniverse(fuelStack[n].getRootUniverse())
    geometry.setCmfd(cmfd)

    # Initialize TrackGenerator
    track_generator = openmoc.TrackGenerator(geometry, num_azim, azim_spacing)
    track_generator.setNumThreads(num_threads)
    track_generator.generateTracks()

    # Run simulation
    solver = openmoc.CPUSolver(track_generator)
    solver.setConvergenceThreshold(tolerance)
    solver.setNumThreads(num_threads)
    solver.computeEigenvalue(max_iters)
    
    # Retrieve the Materials and number of groups from the geometry
    materials = geometry.getAllMaterials()
    num_groups = geometry.getNumEnergyGroups()
    
    # Get microscopic cross sections
    for material_id in materials:
        domains_to_coeffsT[material_id] = np.zeros(num_groups)
        domains_to_coeffsF[material_id] = np.zeros(num_groups)
        domains_to_coeffsPhi[material_id] = np.zeros(num_groups)
        for group in range(num_groups):
            domains_to_coeffsT[material_id][group] = materials[material_id].getSigmaTByGroup(group+1)
            domains_to_coeffsF[material_id][group] = materials[material_id].getSigmaFByGroup(group+1)
            domains_to_coeffsPhi[material_id][group] = 1.0 # gives the flux
        
    # Retrieve macroscopic cross sections
    T_rates = mesh.tally_on_mesh(solver, domains_to_coeffsT, domain_type='material', energy='by_group', volume='integrated')
    F_rates = mesh.tally_on_mesh(solver, domains_to_coeffsF, domain_type='material', energy='by_group', volume='integrated')
    Phi     = mesh.tally_on_mesh(solver, domains_to_coeffsPhi, domain_type='material', energy='by_group', volume='integrated')

    Fission_results.append(F_rates)
    TotalXS_results.append(T_rates)
    Phi_results.append(Phi)
    print('Case :', n, 'finished. -------------------------------------------------')


# In[ ]:

#    Fission_results.append(F_rates)
#    TotalXS_results.append(T_rates)
#    Phi_results.append(Phi)
print(len(TotalXS_results))


# In[ ]:

np.savez('result_210819_'+str(RunID), Fission_results, TotalXS_results, Phi_results)

