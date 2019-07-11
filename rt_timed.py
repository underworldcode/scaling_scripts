import os
os.environ["UW_ENABLE_TIMING"] = "1"
import underworld as uw
from underworld import function as fn
import glucifer
import math
import numpy as np
import time
time_post_import   = time.time()
time_launch_srun   = float(os.environ["TIME_LAUNCH_SRUN"])/1000.
time_launch_python = float(os.environ["TIME_LAUNCH_PYTHON"])/1000.
uw.timing.start()

if os.environ["UW_ENABLE_IO"] == "1":
    do_IO=True
else:
    do_IO=False

other_timing = {}
other_timing["Python_Import_Time"] = time_post_import - time_launch_python
other_timing["Container_Launch_Time"] = time_launch_python - time_launch_srun

res = 64
RESKEY = "UW_RESOLUTION"
if RESKEY in os.environ:
    res = int(os.environ[RESKEY])

PREFIX = os.environ["PREFIXSTRING"]

mesh = uw.mesh.FeMesh_Cartesian(elementRes  = (res, res, res),
                                minCoord    = ( 0., 0., 0., ),
                                maxCoord    = ( 1., 1., 1., ))

velocityField         = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=3 )
pressureField         = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 )
temperatureField      = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )
temperatureFieldDeriv = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=1 )

# initialise 
velocityField.data[:] = [0.,0.,0.]
pressureField.data[:] = 0.

for index, coord in enumerate(mesh.data):
    temperatureField.data[index] = coord[2]

temperatureFieldDeriv.data[:] = 0.


# Create a swarm.
swarm = uw.swarm.Swarm( mesh=mesh )

# Create a data variable. It will be used to store the material index of each particle.
materialIndex = swarm.add_variable( dataType="int", count=1 )

# Create a layout object, populate the swarm with particles.
swarmLayout = uw.swarm.layouts.PerCellSpaceFillerLayout( swarm=swarm, particlesPerCell=40 )
swarm.populate_using_layout( layout=swarmLayout )



# define these for convience. 
denseIndex = 0
lightIndex = 1

# material perturbation from van Keken et al. 1997
wavelength = 2.0
amplitude  = 0.02
offset     = 0.2
k = 2. * math.pi / wavelength

# Create function to return particle's coordinate
coord = fn.coord()

# Define the material perturbation, a function of the x coordinate (accessed by `coord[0]`).
perturbationFn = offset + amplitude*fn.math.cos( k*coord[0] )

# Setup the conditions list. 
# If z is less than the perturbation, set to lightIndex.
conditions = [ ( perturbationFn > coord[1] , lightIndex ),
               (                      True , denseIndex ) ]

# The swarm is passed as an argument to the evaluation, providing evaluation on each particle.
# Results are written to the materialIndex swarm variable.
fnc = fn.branching.conditional( conditions )
matdat = fnc.evaluate(swarm)
materialIndex.data[:] = matdat

store = glucifer.Store('{}_RT'.format(PREFIX),compress=False)

fig = glucifer.Figure( store, name="firstFig" )
fig.append( glucifer.objects.Points(swarm, materialIndex, pointSize=2, colourBar=False) )
fig.append( glucifer.objects.Surface(mesh, pressureField))
fig.append( glucifer.objects.VectorArrows( mesh, velocityField, scaling=1.0e2))


# Set a density of '0.' for light material, '1.' for dense material.
densityMap   = { lightIndex:0., denseIndex:1. }
densityFn    = fn.branching.map( fn_key = materialIndex, mapping = densityMap )

# Set a viscosity value of '1.' for both materials.
viscosityMap = { lightIndex:1., denseIndex:1. }
fn_viscosity  = fn.branching.map( fn_key = materialIndex, mapping = viscosityMap )

# Define a vertical unit vector using a python tuple.
z_hat = ( 0., 0., 1. )

# Create buoyancy force vector
buoyancyFn = -densityFn*z_hat


# Construct node sets using the mesh specialSets
iWalls = mesh.specialSets["MinI_VertexSet"] + mesh.specialSets["MaxI_VertexSet"]
jWalls = mesh.specialSets["MinJ_VertexSet"] + mesh.specialSets["MaxJ_VertexSet"]
kWalls = mesh.specialSets["MinK_VertexSet"] + mesh.specialSets["MaxK_VertexSet"]

allWalls = iWalls + jWalls + kWalls

# Prescribe degrees of freedom on each node to be considered Dirichlet conditions.
# In the x direction on allWalls flag as Dirichlet
# In the y direction on jWalls (horizontal) flag as Dirichlet
stokesBC = uw.conditions.DirichletCondition( variable        = velocityField,
                                             indexSetsPerDof = (allWalls, allWalls, kWalls))
advdiffBc = uw.conditions.DirichletCondition( variable        = temperatureField,
                                              indexSetsPerDof = kWalls )

stokes = uw.systems.Stokes( velocityField = velocityField,
                            pressureField = pressureField,
#                            voronoi_swarm = swarm,
                            conditions    = stokesBC,
                            fn_viscosity  = fn_viscosity, 
                            fn_bodyforce  = buoyancyFn )

solver = uw.systems.Solver( stokes )

# Create a system to advect the swarm
advector = uw.systems.SwarmAdvector( swarm=swarm, velocityField=velocityField, order=2 )

# Create a dummy temperature field.
advdiff = uw.systems.AdvectionDiffusion(velocityField=velocityField, phiField=temperatureField, phiDotField=temperatureFieldDeriv, 
                                        fn_diffusivity=1.,conditions=advdiffBc)


# functions for calculating RMS velocity
vdotv = fn.math.dot(velocityField,velocityField)
v2sum_integral  = uw.utils.Integral( mesh=mesh, fn=vdotv )
volume_integral = uw.utils.Integral( mesh=mesh, fn=1. )


# Get instantaneous Stokes solution
solver.solve()
# Calculate the RMS velocity.
vrms = math.sqrt( v2sum_integral.evaluate()[0] )


# update 
dt1 = advector.get_max_dt()
dt2 = advdiff.get_max_dt()
dt = min(dt1,dt2)
# Advect using this timestep size.
advector.integrate(dt)

advdiff.integrate(dt)

# Save things

if do_IO:
	meshFileHandle = mesh.save("{}_Mesh.h5".format(PREFIX))

	vFH = velocityField.save("{}_velocityField.h5".format(PREFIX))
	velocityField.xdmf( "{}_velocityField".format(PREFIX), vFH, "velocity", meshFileHandle, "Mesh" )

	swarmFileHandle = swarm.save("{}_Swarm.h5".format(PREFIX))
	mH = materialIndex.save("{}_materialIndex.h5".format(PREFIX))
	materialIndex.xdmf("{}_materialIndex".format(PREFIX), mH, "material", swarmFileHandle, "Swarm" )

	fig.save()

	# load things
	# first	 create analogues
	mesh_copy = uw.mesh.FeMesh_Cartesian(
                                 elementRes  = (res, res, res),
                                 minCoord    = (20., 20., 20.),
                                 maxCoord    = (33., 33., 33.))

	velocityField_copy = uw.mesh.MeshVariable( mesh=mesh_copy,         nodeDofCount=3 )

	swarm_copy = uw.swarm.Swarm(mesh = mesh_copy)
	materialIndex_copy = swarm_copy.add_variable( dataType="int", count=1 )

	# now load data and check loaded versions are identical to originals
	mesh_copy.load("{}_Mesh.h5".format(PREFIX))

	# test
	if not np.allclose(mesh_copy.data, mesh.data):
	    raise RuntimeError("Loaded mesh data does not appear to be identical to previous data.")
	velocityField_copy.load("{}_velocityField.h5".format(PREFIX))
	if not np.allclose(velocityField_copy.data, velocityField.data):
	    raise RuntimeError("Loaded velocity data does not appear to be identical to previous data.")


	swarm_copy.load("{}_Swarm.h5".format(PREFIX))

	if not np.allclose(swarm_copy.particleCoordinates.data, swarm.particleCoordinates.data):
		raise RuntimeError("Loaded swarm data does not appear to be identical to previous data.")
	materialIndex_copy.load("{}_materialIndex.h5".format(PREFIX))
	if not np.allclose(materialIndex_copy.data, materialIndex.data):
	    raise RuntimeError("Loaded material data does not appear to be identical to previous data.")

uw.timing.stop()
module_timing_data_orig = uw.timing.get_data(group_by="routine")

# write out data
filename = "{}_Res_{}_Nproc_{}_SlurmID_{}".format(os.environ["SLURM_JOB_NAME"],res,uw.mpi.size,os.environ["SLURM_JOB_ID"])
import json
if module_timing_data_orig:
    module_timing_data = {}
    for key,val in module_timing_data_orig.items():
        module_timing_data[key[0]] = val
    other_timing["Total_Runtime"] = uw.timing._endtime-uw.timing._starttime
    module_timing_data["Other_timing"] = other_timing
    module_timing_data["Other_data"]   = {"vrms":vrms, "res":res, "nproc":uw.mpi.size}
    with open(filename+".json", 'w') as fp:
        json.dump(module_timing_data, fp,sort_keys=True, indent=4)

uw.timing.print_table(group_by="routine", output_file=filename+".txt", display_fraction=0.99)
