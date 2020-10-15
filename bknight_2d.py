import os
os.environ["UW_ENABLE_TIMING"] = "1"
import underworld as uw
import underworld.function as fn
import UWGeodynamics as GEO

order = int(os.getenv("UW_ORDER",1))
res   = int(os.getenv("UW_RESOLUTION",32))
itol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-6))
otol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-6))*10
penalty  = float(os.getenv("UW_PENALTY",1e-3))
do_IO  = bool(int(os.getenv("UW_ENABLE_IO","0")))
jobid = str(os.getenv("PBS_JOBID",os.getenv("SLURM_JOB_ID","0000000")))

import os
from datetime import datetime
import pytz
import numpy as np
import math
import matplotlib.pyplot as plt

import pandas as pd
# import glucifer
import underworld.visualisation as visualisation

uw.timing.start()

# import pickle
# import badlands
# GEO.__version__
# uw.__version__


# In[2]:


### surface processes function
# from UW2DsurfaceProcesses import diffusiveSurfaceErosion, velocitySurfaceErosion


# In[3]:

u = GEO.UnitRegistry

GEO.rcParams["initial.nonlinear.tolerance"] = 1e-2
GEO.rcParams['initial.nonlinear.max.iterations'] = 50

GEO.rcParams["nonlinear.tolerance"] = 1e-2
GEO.rcParams['nonlinear.max.iterations'] = 50

GEO.rcParams["popcontrol.particles.per.cell.2D"] = 30
GEO.rcParams["swarm.particles.per.cell.2D"] = 30


# In[4]:


restart = False


# ### Values to change

# In[5]:


conv_vel = 2.
crustalthickness = 25.0

crust_transition = 700.0
Update_material_LHS_Length = 200.0


# ### Setup of box

# In[6]:


u = GEO.UnitRegistry

Sticky_air = 30.0
x_box = 1200.0
y_box = 300.0 - Sticky_air


Temp_LAB = 1573.15 ### in K, 1300 in C
Bottom_temp = Temp_LAB + (380 * 0.5) ### Base of model temp


Depth_of_box = y_box * u.kilometer
model_length = x_box * u.kilometer

x_res = res*4 #2**4 #int(x / 10.)
y_res = res   #2**4 #int(y / 10.)


resolution = (x_box/x_res)
print(resolution)


# In[7]:


Total_Convergence = 600.0
Total_Time = Total_Convergence / (10*conv_vel)
The_Checkpoint_interval = Total_Time / 40.


# In[8]:


if uw.mpi.rank == 0:
    print(x_res,y_res)
    print(Total_Time)
    print(The_Checkpoint_interval)

uw.mpi.barrier()


# ## Scaling

# In[9]:


half_rate = 1. * u.centimeter / u.year
length_scale = 300.0 * u.kilometer
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1573.15 * u.degK
bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = length_scale
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT



### Scaling
if uw.mpi.rank == 0:
    print('Length, km = ', GEO.dimensionalise(1., u.kilometer))
    print('Time, Myr = ',GEO.dimensionalise(1., u.megayear))
    print('Pressure, MPa = ',GEO.dimensionalise(1., u.megapascal))
    print('Temperature, K = ',GEO.dimensionalise(1., u.degK))
    print('Mass, kg = ',GEO.dimensionalise(1.,u.kilogram))

    print('Velocity, cm/yr = ',GEO.dimensionalise(1., u.centimeter / u.year))
    print('Diffusivity, m^2/s = ',GEO.dimensionalise(1.,u.metre**2 / u.second))
    print('Density, kg/m^3 = ',GEO.dimensionalise(1.,u.kilogram / u.metre**3))
    print('Viscosity, Pa S = ',GEO.dimensionalise(1.,u.pascal * u.second))
    print('gravity, m/s^2 = ',GEO.dimensionalise(1.,u.meter / u.second**2))


# # Define the external geometry
#
# The first step is to define the geometry of our problem, essentially a box on which we will apply some physical constraints and that will contain a set of materials. We can think of it as an "universe".
# The "laws" and their effects are calculated on a mesh, that mesh discretized our universe into finite elements.
#
# The geodynamics module allows you to quickly set up a model by creating a *Model* object.
# A series of arguments are required to define a *Model*:
#
#     - The number of elements in each direction elementRes=(nx, ny);
#     - The minimum coordinates for each axis minCoord=(minX, minY)
#     - The maximum coordinates for each axis maxCoord=(maxX, maxY)
#     - A vector that defines the magnitude and direction of the gravity components gravity=(gx, gy)
#

# In[10]:


Model = GEO.Model(elementRes=(x_res, y_res),
                  minCoord=(0. * u.kilometer, -1.*Depth_of_box),
                  maxCoord=(model_length, Sticky_air * u.kilometer),
                  gravity=(0.0, -9.81 * u.meter / u.second**2))


# In[11]:


today = datetime.now(pytz.timezone('Australia/Melbourne'))
if restart == False:
    Model.outputDir = os.path.join("out_"+jobid)

    directory = os.path.join(os.getcwd(),Model.outputDir) 

    if uw.mpi.rank == 0:
        if not os.path.exists(directory):
            os.makedirs(directory)
    uw.mpi.barrier()

if restart == True:
    RestartDirectory = os.getcwd()
    directory = RestartDirectory


# # Add some additional swarm/mesh variables

# In[12]:


swarmStrainRateInv = Model.add_swarm_variable(name='SR_swarm', count=1, restart_variable=True)
swarmDisplacement = Model.add_swarm_variable(name='swarmDisplacement', count=1, restart_variable=True)


meshDisplacement = Model.add_mesh_variable(name='meshDisplacement', nodeDofCount=1)





meshViscousDissipation = Model.add_mesh_variable(name='VD_mesh', nodeDofCount=1 )


# In[13]:


Model.diffusivity = 1e-6 * u.metre**2 / u.second
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# # Add some Materials

# In[14]:


air = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0. * u.kilometer))


# In[15]:


## Upper crust
crust1 = Model.add_material(name="Crust1")
crust2 = Model.add_material(name="Crust2")
crust3 = Model.add_material(name="Crust3")
crust4 = Model.add_material(name="Crust4")


# In[16]:


#x_threshold = (GEO.nd(300.0*u.kilometer))

#Model.plasticStrain.data[x<x_threshold] = 0.



Fault_positionX = crust_transition * u.kilometer


sin_function = np.sign(np.sin(GEO.dimensionalise(Model.swarm.data[:,1], u.kilometer)/(1.6 * u.kilometer)))

Model.materialField.data[(sin_function>0) & (Model.swarm.data[:,0] < GEO.nd(Fault_positionX+15.*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer))] = crust1.index
Model.materialField.data[(sin_function<0) & (Model.swarm.data[:,0] < GEO.nd(Fault_positionX+15.*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer))] = crust2.index

Model.materialField.data[(sin_function>0) & (Model.swarm.data[:,0] >= GEO.nd(Fault_positionX+15.*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer))] = crust3.index
Model.materialField.data[(sin_function<0) & (Model.swarm.data[:,0] >= GEO.nd(Fault_positionX+15.*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer))] = crust4.index


# In[17]:


mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer(top=-crustalthickness* u.kilometer, bottom=-100.0 * u.kilometer))
mantle = Model.add_material(name="Mantle", shape=GEO.shapes.Layer(top=mantleLithosphere.bottom, bottom=Model.bottom))


# In[18]:


Sediment = Model.add_material(name="Sediment")


# In[19]:


def Update_Material_LHS():
    sin_function = np.sign(np.sin(GEO.dimensionalise(Model.swarm.data[:,1], u.kilometer)/(1.6 * u.kilometer)))
    Model.materialField.data[(Model.swarm.data[:,0] < GEO.nd(Update_material_LHS_Length*u.kilometer)) &  (Model.swarm.data[:,1] > GEO.nd(0.*u.kilometer)) ] = air.index
    Model.materialField.data[(sin_function>0) & (Model.swarm.data[:,0] < GEO.nd(Update_material_LHS_Length*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer)) &  (Model.swarm.data[:,1] >= GEO.nd(-1.*crustalthickness*u.kilometer))] = crust1.index
    Model.materialField.data[(sin_function<0) & (Model.swarm.data[:,0] < GEO.nd(Update_material_LHS_Length*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer)) &  (Model.swarm.data[:,1] >= GEO.nd(-1.*crustalthickness*u.kilometer))] = crust2.index
    Model.materialField.data[(Model.swarm.data[:,0] < GEO.nd(Update_material_LHS_Length*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(-1.*crustalthickness*u.kilometer)) &  (Model.swarm.data[:,1] >= GEO.nd(mantleLithosphere.bottom))] = mantleLithosphere.index
    Model.materialField.data[(Model.swarm.data[:,0] < GEO.nd(Update_material_LHS_Length*u.kilometer)) &  (Model.swarm.data[:,1] < GEO.nd(mantleLithosphere.bottom))] = mantle.index


# ### Material physical properties

# In[20]:


air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
# stickyAir.capacity = 100. * u.joule / (u.kelvin * u.kilogram)


# In[21]:


air.density = 1. * u.kilogram / u.metre**3
# stickyAir.density = 1. * u.kilogram / u.metre**3

mantleLithosphere.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantle.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)


mantleLithosphere.radiogenicHeatProd = 0.00 * u.microwatt / u.meter**3
mantle.radiogenicHeatProd = 0.00 * u.microwatt / u.meter**3


# In[22]:


crust1.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
crust2.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
crust3.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
crust4.radiogenicHeatProd = 0. * u.microwatt / u.meter**3

# crust5.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
# crust6.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
# crust7.radiogenicHeatProd = 0. * u.microwatt / u.meter**3
# crust8.radiogenicHeatProd = 0. * u.microwatt / u.meter**3

crust1.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
crust2.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
crust3.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
crust4.density = GEO.LinearDensity(2700. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

# crust5.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
# crust6.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
# crust7.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
# crust8.density = GEO.LinearDensity(2900. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

Sediment.density = GEO.LinearDensity(2600. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
Sediment.radiogenicHeatProd = 0. * u.microwatt / u.meter**3


# In[23]:


angle_D = 45 # in degrees
Fault_PositionX_LAB = Fault_positionX + ((mantleLithosphere.top - mantleLithosphere.bottom) * math.tan(math.radians(angle_D)))

FaultShape = GEO.shapes.Polygon(vertices=[(Fault_positionX, mantleLithosphere.top),
                                        (Fault_positionX + 35.* u.kilometer, mantleLithosphere.top),
                                        (Fault_PositionX_LAB+3.5* u.kilometer, mantleLithosphere.bottom),
                                         (Fault_PositionX_LAB, mantleLithosphere.bottom)])

Fault = Model.add_material(name="Fault", shape=FaultShape)
Fault.radiogenicHeatProd = 0.00 * u.microwatt / u.meter**3
Fault.density  = GEO.LinearDensity(reference_density=3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)


# In[24]:


Fig = visualisation.Figure(figsize=(1200,400), title="Material Field", quality=2)
# #Fig.Points(surface_tracers, pointSize=2.0)
# # Fig.Points(moho_tracers, pointSize=2.0, colour="red")
# # Fig.Points(FSE_Crust, pointSize=2.0)
# # Fig.Points(FSE_Mantle, pointSize=2.0)
Fig.Points(Model.swarm, Model.materialField, fn_size=2.0)
Fig.show()


# # Define Viscosities
#
# The rheology library contains some commonly used rheologies stored in a python dictionary structure. We can list the keys defining the rheologies as follows:

# In[25]:


rh = GEO.ViscousCreepRegistry()


# ### Values from Gerya (2009), page 271

# In[26]:


WQ_Ranalli_1995 = GEO.ViscousCreep(preExponentialFactor=3.2e-4/u.megapascal**2.3/u.second,
                                                      stressExponent=2.3,
                                                      activationVolume=0.,
                                                      activationEnergy=154 * u.kilojoules/u.mole,
                                                      waterFugacity=0.0,
                                                      grainSize=0.0,
                                                      meltFraction=0.,
                                                      grainSizeExponent=0.,
                                                      waterFugacityExponent=0.,
                                                      meltFractionFactor=0.0,
                                                      f=1.0,
                                                   name='Wet quartzite')


# In[27]:


Diabase_Dislocation_Ranalli_1995 = GEO.ViscousCreep(preExponentialFactor=2.0e-4/u.megapascal**3.4/u.second,
                                                      stressExponent=3.4,
                                                      activationVolume=0.,
                                                      activationEnergy=260 * u.kilojoules/u.mole,
                                                      waterFugacity=0.0,
                                                      grainSize=0.0,
                                                      meltFraction=0.,
                                                      grainSizeExponent=0.,
                                                      waterFugacityExponent=0.,
                                                      meltFractionFactor=0.0,
                                                      f=1.0,
                                                   name='Diabase',
                                                    mineral='Diabase')


# In[28]:



# In[29]:


# ### Weak Zone Rheology
#rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
#rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003
# rh.Wet_Olivine_Dislocation_Karato_and_Wu_1993

# ### Mantle Rheology
# rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
# rh.Dry_Olivine_Diffusion_Hirth_and_Kohlstedt_2003
# rh.Dry_Olivine_Dislocation_Hirth_and_Kohlstedt_2003

# ### Strong Crust
# rh.Dry_Quartz_Dislocation_Koch_et_al_1983
# rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998

# ### Weak Crust
# rh.Wet_Quartz_Dislocation_Tullis_et_al_2002
#rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990

# rh.Wet_Anorthite_Dislocation_Ribacki_et_al_2000

combined_viscosity_mantle = GEO.CompositeViscosity([rh.Dry_Olivine_Diffusion_Hirth_and_Kohlstedt_2003, rh.Dry_Olivine_Dislocation_Hirth_and_Kohlstedt_2003])
#rh.Dry_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
combined_viscosity_fault = GEO.CompositeViscosity([rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003, rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003])


#rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
#rh.Wet_Olivine_Dislocation_Karato_and_Wu_1993


# In[30]:


rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998


# In[31]:


Model.minViscosity = 1e19 * u.pascal * u.second
Model.maxViscosity = 1e24 * u.pascal * u.second

air.viscosity                = 1e19 * u.pascal * u.second


mantleLithosphere.viscosity  = combined_viscosity_mantle
mantle.viscosity             = combined_viscosity_mantle

Fault.viscosity              = combined_viscosity_fault


### Crust
crust1.viscosity              = rh.Wet_Quartz_Dislocation_Tullis_et_al_2002
crust2.viscosity              = rh.Wet_Quartz_Dislocation_Tullis_et_al_2002
crust3.viscosity              = rh.Wet_Quartz_Dislocation_Tullis_et_al_2002
crust4.viscosity              = rh.Wet_Quartz_Dislocation_Tullis_et_al_2002

Sediment.viscosity            = 0.5 * WQ_Ranalli_1995


# # Define Plasticity
#
# Plastic behavior is assigned using the same approach as for viscosities.

# In[32]:


pl = GEO.PlasticityRegistry()


# In[33]:


# Weak_crust_plasticity = DruckerPrager_SRweakening(
#     SR=Model.strainRate_2ndInvariant,
#     cohesion=10.0 * u.megapascal,
#     frictionCoefficient=0.3,
#     frictionAfterSoftening=0.15,
#     cohesionAfterSoftening=1.0 * u.megapascal,
#     StrainRateThreshold=1.0e-15*1./u.second,
# )


# Sediment_plasticity = DruckerPrager_SRweakening(
#     SR=Model.strainRate_2ndInvariant,
#     cohesion=10.0 * u.megapascal,
#     cohesionAfterSoftening=1.0 * u.megapascal,
#     frictionCoefficient=0.3,
#     frictionAfterSoftening=0.15,
#     StrainRateThreshold=1.0e-15*1./u.second)


# In[34]:

veryweak_crust_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.1,
                                     frictionAfterSoftening=0.05,
                                     epsilon1=0.5, epsilon2=1.5)

crust_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.3,
                                     frictionAfterSoftening=0.15,
                                     epsilon1=0.5, epsilon2=1.5)


Weak_crust_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.2,
                                     frictionAfterSoftening=0.1,
                                     epsilon1=0.5, epsilon2=1.5)

Strong_crust_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.4,
                                     frictionAfterSoftening=0.2,
                                     epsilon1=0.5, epsilon2=1.5)

Sediment_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.2,
                                     frictionAfterSoftening=0.1,
                                     epsilon1=0.5, epsilon2=1.5)



Mantle_plasticity =  GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=10.*u.megapascal,
                                     frictionCoefficient=0.6,
                                     frictionAfterSoftening=0.6,
                                     epsilon1=0.5, epsilon2=1.5)

Fault_plasticity = GEO.DruckerPrager(cohesion=10.* u.megapascal,
                                     cohesionAfterSoftening=1.*u.megapascal,
                                     frictionCoefficient=0.1,
                                     frictionAfterSoftening=0.05,
                                     epsilon1=0.5, epsilon2=1.5)


# In[35]:


mantleLithosphere.plasticity  = Mantle_plasticity
mantle.plasticity             = Mantle_plasticity
Fault.plasticity              = Fault_plasticity


# In[36]:


# Ucrust.plasticity              = Weak_crust_plasticity
# Lcrust.plasticity              = Weak_crust_plasticity

crust1.plasticity              = veryweak_crust_plasticity
crust2.plasticity              = veryweak_crust_plasticity
crust3.plasticity              = veryweak_crust_plasticity
crust4.plasticity              = veryweak_crust_plasticity
### Lower Crust
# crust5.plasticity              = Weak_crust_plasticity
# crust6.plasticity              = Weak_crust_plasticity
# crust7.plasticity              = Weak_crust_plasticity
# crust8.plasticity              = Weak_crust_plasticity

Sediment.plasticity            = Sediment_plasticity


# ## Temperature Boundary Conditions

# In[37]:


Model.set_temperatureBCs(top=273.15 * u.degK,
                         bottom=1573.15 * u.degK,
                         materials=[(air, 273.15 * u.degK)])


# ## Velocity Boundary Conditions

# In[38]:


Model.velocityField.data[:] = 0.


# In[39]:


# FigVel = visualisation.Figure(figsize=(1200,400), title="Velocity")
# # Fig.Surface(Model.mesh, GEO.dim(Model.temperature, u.degK))
# FigVel.VectorArrows(Model.mesh, Model.velocityField)
# # FigVel.show()


# In[40]:



def UpdateVelocity():
    global conv_vel

    conv_vel = conv_vel * u.centimeter/u.year

    conditionsA = [(Model.y < GEO.nd(0. * u.kilometre), GEO.nd(conv_vel)),
                   (True, GEO.nd(conv_vel) + Model.y * (GEO.nd((-2. * conv_vel) / GEO.nd(Sticky_air * u.kilometer))))]


    Left_wall_vel_top_changed = fn.branching.conditional(conditionsA)

## Test air boundary condition on side wall


    conditionsB = [(Model.y > GEO.nd(mantleLithosphere.bottom), Left_wall_vel_top_changed),
                   (True, (GEO.nd(conv_vel) + (Model.y-GEO.nd(mantleLithosphere.bottom)) * (GEO.nd(conv_vel) / GEO.nd(Depth_of_box+mantleLithosphere.bottom))))]


    Left_wall_vel = fn.branching.conditional(conditionsB)


    Model.set_velocityBCs(left = [Left_wall_vel, None],
                      right=[0., None],
                      top = [None, 0.],
                        bottom = [None, None])


UpdateVelocity()



# ## Initialize plastic strain

# ### Interface tracers

# In[41]:


# ### topography tracer every 1. km

# npoints = int(Model.maxCoord[0].to(u.kilometer).magnitude*(1/1))

# ### creates surface
# # x_surface = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints)
# # y_surface = 0. * u.kilometer

# coords = np.ndarray((npoints, 2))
# coords[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints)
# coords[:, 1] = GEO.nd(0. * u.kilometer)

# surface_tracers = Model.add_passive_tracers(name="SurfaceTracers", vertices=[coords[:, 0], coords[:, 1].min()])


# ### Grid Tracers

# In[42]:


coords_circ = GEO.circles_grid(radius=5.0*u.kilometer,
                    minCoord=[Model.minCoord[0], mantleLithosphere.top],
                    maxCoord=[Model.maxCoord[0], air.top])


FSE_Crust = Model.add_passive_tracers(name="FSE_Crust", vertices=coords_circ)

FSE_Crust.add_tracked_field(Model.pressureField,
                              name="tracers_press",
                              units=u.megapascal,
                              dataType="double")

FSE_Crust.add_tracked_field(Model.temperature,
                              name="tracers_temp",
                              units=u.degK,
                              dataType="double")

FSE_Crust.add_tracked_field(Model.strainRate_2ndInvariant,
                              name="tracers_SR",
                                units=1./u.second,
                              dataType="double")


# In[43]:


Fig = visualisation.Figure(figsize=(1200,400), title="Material Field", quality=0)
# Fig.Points(surface_tracers, pointSize=2.0, colour="red")
# Fig.Points(FSE_Crust, pointSize=2.0)
# Fig.Points(FSE_Mantle, pointSize=2.0)
Fig.Points(Model.swarm, Model.materialField, fn_size=2.0)
Fig.show()


# In[44]:


#   # # ### reasigns the swarm variables based on the diffusive surface
# #    Model.materialField.data[(Model.swarm.data[:,1] > GEO.nd(f(GEO.dimensionalise(swarm1[:,0], u.kilometer) / u.kilometer)*u.kilometer)) & (Model.materialField.data[:,0] > air.index)] = air.index
# #    Model.materialField.data[(Model.swarm.data[:,1] < GEO.nd(f(GEO.dimensionalise(swarm1[:,0], u.kilometer) / u.kilometer)*u.kilometer)) & (Model.materialField.data[:,0] == air.index)] = Sediment.index
# print((Model.swarm.data[:,1] > GEO.nd(f(GEO.dimensionalise(Model.swarm.data[:,0], u.kilometer) / u.kilometer)*u.kilometer)))


# Viscous dissipation calc on the swarm and mesh, for integration over the crust area

# In[45]:


def update_custom_fields():
    ## For projecting swarm variables onto mesh
    # Model.swarm.allow_parallel_nn = True

    # swarmViscousDissipation.data[:] = GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate(Model.swarm), 1./u.second) * GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate(Model.swarm), 1./u.second) * GEO.dimensionalise(Model.viscosityField.data[:], u.pascal * u.second)
    ### Calculation on mesh
    #Model.meshViscousDissipation.data[:] = 2. * Model.strainRate_2ndInvariant.evaluate( Model.mesh) * Model.strainRate_2ndInvariant.evaluate( Model.mesh) * Model.viscosityField.evaluate( Model.mesh)
    meshViscousDissipation.data[:] = GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate(Model.mesh), 1./u.second) * GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate(Model.mesh), 1./u.second) * GEO.dimensionalise(Model.projViscosityField.data[:], u.pascal * u.second)

    swarmStrainRateInv.data[:] = GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate( Model.swarm), 1./u.second)
    # viscositySwarm.data[:] = GEO.dimensionalise(Model.viscosityField.evaluate( Model.swarm), u.pascal * u.second)

    # print('mesh VD data: ', meshViscousDissipation.data[1:10])

update_custom_fields()


# In[46]:


meshDisplacement.data[:] = 0.
swarmDisplacement.data[:] = 0.

def Displacement_calc():
    vel_mesh_data = GEO.dimensionalise(Model.velocityField.evaluate(Model.mesh), u.kilometer/u.year).magnitude

    vel_mag = np.sqrt(vel_mesh_data[:,0]**2 + vel_mesh_data[:,1]**2)

    SR_swarm_data = GEO.dimensionalise(Model.strainRate_2ndInvariant.evaluate(Model.mesh), 1./u.year).magnitude

    with np.errstate(divide='ignore'):
        displacement_measurement = np.divide(SR_swarm_data[:,0], vel_mag[:])
        displacement_measurement[displacement_measurement == np.inf] = 0

    meshDisplacement.data[:,0] += displacement_measurement

    swarmDisplacement.data[:] = meshDisplacement.evaluate(Model.swarm)








# In[47]:


# Fig = visualisation.Figure(figsize=(1200,400), title="Displacement", quality=0)
# # Fig.Points(moho_tracers, pointSize=2.0, colour="red")
# # Fig.Points(FSE_Crust, pointSize=2.0)
# # Fig.Points(FSE_Mantle, pointSize=2.0)
# Fig.Surface(Model.mesh, meshDisplacement)
# # Fig.Points(Model.mesh, meshDisplacement)
# Fig.show()


# In[48]:


VD_Ucrust = []
VD_Lcrust = []
VD_Sed    = []
VD_model  = []

def viscous_dissipation_calc():
    VD_model_df = pd.DataFrame()
    c1 = Model.materialField > (crust1.index - 0.5)
    c2 = Model.materialField < (crust2.index + 0.5)

    lower_plate = c1 & c2

    c3 = Model.materialField > (crust3.index - 0.5)
    c4 = Model.materialField < (crust4.index + 0.5)

    upper_plate = c3 & c4

    model_material = Model.materialField > (crust1.index - 0.5)

    s1 = Model.materialField > (Sediment.index - 0.5)
    s2 = Model.materialField < (Sediment.index + 0.5)

    sed = s1 & s2


    # how to calculate VD
    vd = Model._viscosityFn * Model.strainRate_2ndInvariant**2
    # vd_dimensionalised = GEO.dimensionalise(vd, u.pascal / u.second)

    """lower plate VD"""
    clause = [ (lower_plate, vd),
                ( True   , 0.) ]

    fn_crust_vd = fn.branching.conditional( clause )

    VD_Lcrust.append(2. * Model.mesh.integrate(fn_crust_vd)[0])
    """upper plate VD"""
    clause = [ (upper_plate, vd),
                ( True   , 0.) ]

    fn_crust_vd = fn.branching.conditional( clause )

    VD_Ucrust.append(2. * Model.mesh.integrate(fn_crust_vd)[0])


    """Sediment VD"""
    clause = [ (sed, vd),
                ( True   , 0.) ]

    fn_crust_vd = fn.branching.conditional( clause )

    VD_Sed.append(2. * Model.mesh.integrate(fn_crust_vd)[0])


    '''Model VD'''

    clause = [ (model_material, vd),
                ( True   , 0.) ]

    fn_model_vd = fn.branching.conditional( clause )

    VD_model.append(2. * Model.mesh.integrate(fn_model_vd)[0])

    '''Save VD to DF'''
    VD_model_df['Lower Crust'] = VD_Lcrust
    VD_model_df['Upper Crust'] = VD_Ucrust
    VD_model_df['Sediment']    = VD_Sed

    VD_model_df['Total']       = VD_model

    VD_model_df['Lower Crust %'] = (VD_model_df['Lower Crust'] / VD_model_df['Total'] * 100)
    VD_model_df['Upper Crust %'] = (VD_model_df['Upper Crust'] / VD_model_df['Total'] * 100)
    VD_model_df['Sediment %']    = (VD_model_df['Sediment']    / VD_model_df['Total'] * 100)

    VD_model_df['Total %']       = (VD_model_df['Total'] / VD_model_df['Total'] *100)



    VD_model_df.to_csv(os.path.join(directory, 'VD_data.csv'))


# # Compute initial condition

# In[49]:


Model.init_model()


# In[50]:


### Custom temp gradient

for index, coord in enumerate(Model.mesh.data):
### Temperature in air
    if coord[1] > 0.:
        T = (273.15 * u.kelvin)
    #### Temperature across top 10 km of crust
    elif coord[1] < 0. and coord[1] >= GEO.nd(-10*u.kilometer):
            T = (273.15 * u.kelvin + (-1*GEO.dimensionalise(coord[1], u.kilometer) * 25. * u.kelvin/u.kilometer))
    #### Temperature for the lower crust and lithosphere
    elif coord[1] < GEO.nd(-10*u.kilometer) and coord[1] >= GEO.nd(mantleLithosphere.bottom):
            T = ((273.15+130.0) * u.kelvin + (-1*GEO.dimensionalise(coord[1], u.kilometer) * 12. * u.kelvin/u.kilometer))
#### Temperature for the Upper Mantle
    elif coord[1] < GEO.nd(mantleLithosphere.bottom):
        T = 1573.15* u.degK #(1573.15 * u.kelvin + (-1*GEO.dimensionalise(coord[1], u.kilometer) * 0.5 * u.kelvin/u.kilometer))

    Model.temperature.data[index] = GEO.nd(T)



# In[51]:


# # Only run this when in serial. Will fail in parallel

# if GEO.nProcs == 1:

#     distances, viscosities = GEO.extract_profile(Model.projViscosityField,
#                                                  line = [(180.* u.kilometer, 0.), (180.* u.kilometer, Model.bottom)])
#     distances, stresses = GEO.extract_profile(Model.projStressField,
#                                              line = [(180.* u.kilometer, 0.), (180.* u.kilometer, Model.bottom)])


#     Fig, (ax1) = plt.subplots(1,1,figsize=(7,7))
#     ax1.plot([0,1e23], [20,20])
#     ax1.plot(GEO.dimensionalise(viscosities, u.pascal * u.second), GEO.dimensionalise(distances, u.kilometer))
#     ax1.set_xlabel("Viscosity in Pa.s")
#     ax1.set_ylabel("Depth in kms")
#     ax1.set_ylim(100,-10)
#     ax1.set_title("Viscosity profile")


# ### Additional stuff at checkpointing

# In[ ]:


if restart == True:
    TotalConvergence = 15.312948351157392 * u.kilometer
else:
    TotalConvergence = 0.

def Additional_files():

    file_prefix = os.path.join(directory, 'viscosityField-%s' % Model.checkpointID)
    handle = Model.viscosityField.save('%s.h5' % file_prefix)


def CumulativeStrainCheck():
    Model.plasticStrain.data[Model.swarm.data[:,0]<(GEO.nd(100.0*u.kilometer))] = 0.

def Checkpoint_additional_stuff():
#     global TotalConvergence
#     TotalConvergence += GEO.dimensionalise((GEO.nd((conv_vel* u.centimeter / u.year) * GEO.dimensionalise(Model._dt, u.megayear))), u.kilometer)


    ### Stops strain on new incoming materials

    # CumulativeStrainCheck()

    if GEO.nd(round(Model.time,0) % round(The_Checkpoint_interval * 1e6 * u.years, 0)) == 0.:

        Additional_files()
        viscous_dissipation_calc()









# In[ ]:


# ### if superlu (local) /superludist (HPC)/mumps (local/HPC), high penalty
Model.solver.set_inner_method("mg")
Model.solver.set_inner_rtol(itol)
Model.solver.set_outer_rtol(otol)
Model.solver.set_penalty(penalty)


# In[ ]:


### additional functions for the model
Model.pre_solve_functions["A-pre"] = CumulativeStrainCheck
Model.pre_solve_functions["B-pre"] = Update_Material_LHS
Model.pre_solve_functions["C-pre"] = Checkpoint_additional_stuff
Model.pre_solve_functions["D-pre"] = update_custom_fields
### Erosion functions

# Model.post_solve_functions["A-post"] = Erosion_multiple_cpu
Model.post_solve_functions["B-post"] = Displacement_calc
Model.post_solve_functions["C-post"] = Update_Material_LHS


# In[ ]:


# Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[air.index], sedimentIndex=[Sediment.index],
#                                           XML="badlands-erosion.xml", resolution=10. * u.kilometer,
#                                           checkpoint_interval=The_Checkpoint_interval * u.megayears)


# In[ ]:

#
# Model.surfaceProcesses = diffusiveSurfaceErosion(
#     airIndex=air.index,
#     sedimentIndex=Sediment.index,
#     D= 100.0*u.meter**2/u.year,
#     dx=2 * u.kilometer,
#     surfaceElevation=air.bottom
# )


# In[ ]:


GEO.rcParams["default.outputs"].append("SR_swarm")

GEO.rcParams["default.outputs"].append("VD_mesh")

GEO.rcParams["default.outputs"].append("viscosityField")

GEO.rcParams["default.outputs"].append("meshDisplacement")

GEO.rcParams["default.outputs"].append("swarmDisplacement")


# In[ ]:


GEO.rcParams["default.outputs"].remove("projMeltField")

GEO.rcParams["default.outputs"].remove("projDensityField")

GEO.rcParams["default.outputs"].remove("projPlasticStrain")


# # Run the Model

# In[ ]:


if restart == True:
    Model.run_for(Total_Time * u.megayears, checkpoint_interval=The_Checkpoint_interval*u.megayears, restartStep=-1, restartDir=RestartDirectory)
else:
    Model.run_for(nstep=1)





# In[ ]:


uw.timing.stop()
module_timing_data_orig = uw.timing.get_data(group_by="routine")

# write out data
filename = "Res_{}_Nproc_{}_JobID_{}".format(res,uw.mpi.size,jobid)
import json
if module_timing_data_orig:
    module_timing_data = {}
    for key,val in module_timing_data_orig.items():
        module_timing_data[key[0]] = val
    module_timing_data["Other_data"]   = { "res":res, "nproc":uw.mpi.size }
    with open(filename+".json", 'w') as fp:
        json.dump(module_timing_data, fp,sort_keys=True, indent=4)

uw.timing.print_table(group_by="routine", output_file=filename+".txt", display_fraction=0.99)
