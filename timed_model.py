import os
os.environ["UW_ENABLE_TIMING"] = "1"
import underworld as uw
from underworld import function as fn
import underworld.visualisation as vis
import numpy as np
import time
import math

order = int(os.getenv("UW_ORDER","2"))
res   = int(os.getenv("UW_RESOLUTION",16))
dim   = int(os.getenv("UW_DIM",3))

itol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-10))
otol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-10))*10
penalty  = float(os.getenv("UW_PENALTY",-1.))
max_its  = int(os.getenv("UW_MAX_ITS",-1))

soln_name = str(os.getenv("UW_MODEL","SolDB3d"))

do_IO  = bool(int(os.getenv("UW_ENABLE_IO","0")))

jobid = str(os.getenv("PBS_JOBID",os.getenv("SLURM_JOB_ID","0000000")))

picklename = str(os.getenv("PICKLENAME","None"))

# Find all available solutions. 
# Use ordered dict to preserve alphabetical ordering
import collections
solns_avail = collections.OrderedDict()
for _soln in dir(fn.analytic):
    if _soln[0] == "_": continue  # if private member, ignore
    # get soln class
    soln = getattr(fn.analytic,_soln)
    # check if actually soln
    if issubclass(soln, fn.analytic._SolBase):
        solns_avail[_soln] = soln
soln = solns_avail[soln_name]()


time_post_import   = time.time()
time_launch_mpi    = float(os.getenv("TIME_LAUNCH_MPI"   ,time_post_import))/1000.
time_launch_python = float(os.getenv("TIME_LAUNCH_PYTHON",time_post_import))/1000.

uw.timing.start()

other_timing = {}
other_timing["Python_Import_Time"] = time_post_import - time_launch_python
other_timing[   "MPI_Launch_Time"] = time_launch_python - time_launch_mpi

def rms_error(numeric, analytic, mesh, nodal_errors=False, normalise=False):
    '''
    Calculates the rms error.
    
    Returns
    -------
    abs, abs_scaled: float
        The absolute and scaled absolute errors.
    '''
    partCountMap = { "DQ0"  : 1,
                     "Q1"   : 2,
                     "DQ1"  : 2,
                     "DPC1" : 2,
                     "Q2"   : 3  }
    particleCount = partCountMap[ numeric.mesh.elementType.upper() ]
    intSwarm = uw.swarm.GaussIntegrationSwarm(mesh,particleCount) 
    if nodal_errors:
        analyticsoln = analytic # grab handle before replacing
        analytic = numeric.mesh.add_variable(numeric.nodeDofCount)
        analytic.data[:] = analyticsoln.evaluate(numeric.mesh) # careful here to eval on corresponding mesh
    
    if normalise:
        numeric  -= uw.utils.Integral( numeric,  mesh, integrationSwarm=intSwarm, integrationType=None).evaluate()[0]
        analytic -= uw.utils.Integral( analytic, mesh, integrationSwarm=intSwarm, integrationType=None).evaluate()[0]
 
    delta     = analytic - numeric
    delta_dot = fn.math.dot(delta,delta)
    analytic_dot = fn.math.dot(analytic,analytic)

    # l2 norms
    rms_err_abs = np.sqrt(uw.utils.Integral(    delta_dot, mesh, integrationSwarm=intSwarm, integrationType=None ).evaluate()[0])
    rms_sol_ana = np.sqrt(uw.utils.Integral( analytic_dot, mesh, integrationSwarm=intSwarm, integrationType=None ).evaluate()[0])
    rms_err_sca = rms_err_abs / rms_sol_ana
        
    return rms_err_abs, rms_err_sca


if order == 1:
    els = "Q1/dQ0"
elif order == 2:
    els = "Q2/dPc1"
else:
    raise ValueError("Provided system order should be 1 or 2.")

mesh          = uw.mesh.FeMesh_Cartesian(elementType=els, elementRes=(res,)*dim,minCoord=(0.,)*dim,maxCoord=(1.,)*dim)
velocityField = uw.mesh.MeshVariable(mesh, dim)
pressureField = uw.mesh.MeshVariable(mesh.subMesh, 1)

velocityField.data[:] = (0.,)*dim
pressureField.data[:] = 0.

bcs = soln.get_bcs(velocityField)
visc = soln.fn_viscosity
if soln.nonlinear==True:
    visc = soln.get_viscosity_nl(vel,press)
stokes = uw.systems.Stokes(velocityField, pressureField, fn_viscosity=visc, fn_bodyforce=soln.fn_bodyforce, conditions=[bcs,])
solver = uw.systems.Solver(stokes)

if max_its < 0:
    solver.set_inner_rtol(itol)
    solver.set_outer_rtol(otol)
else:
    solver.set_inner_method("nomg")
    solver.options.A11.ksp_rtol=1e-99; solver.options.A11.ksp_max_it = max_its
    solver.options.scr.ksp_rtol=1e-99; solver.options.scr.ksp_max_it = max_its

if soln.nonlinear==True:
    stokes.fn_viscosity = 1.
    solver.solve()
    stokes.fn_viscosity = visc
if penalty>=0.:
    solver.set_penalty(penalty)

# functions for calculating RMS velocity
vdotv = fn.math.dot(velocityField,velocityField)
v2sum_integral  = uw.utils.Integral( mesh=mesh, fn=vdotv )
volume_integral = uw.utils.Integral( mesh=mesh, fn=1. )

pdotp = fn.math.dot(pressureField,pressureField)
pressure_int = uw.utils.Integral( mesh=mesh, fn=pdotp )

# Get instantaneous Stokes solution
solver.solve()

stats=solver.get_stats()
solver.print_stats()

# If using `max_its`, velocity solution should not be 
# trust, so switch to fixed velocity field.
if max_its >= 0:
    for index,coord in enumerate(mesh.data):
        transcoord = coord - (0.5,0.5,0.)
        if   (coord[0] >= 0.) and (coord[1] >= 0.):
            velocityField.data[index] = (-1.,  0.,  0.)
        elif (coord[0] <= 0.) and (coord[1] >= 0.):
            velocityField.data[index] = ( 0., -1.,  0.)
        elif (coord[0] <= 0.) and (coord[1] <= 0.):
            velocityField.data[index] = ( 1.,  0.,  0.)
        elif (coord[0] >= 0.) and (coord[1] <= 0.):
            velocityField.data[index] = ( 0.,  1.,  0.)
        else:
            velocityField.data[index] = ( 0.,  0.,  0.)

# Calculate the RMS velocity.
vrms = math.sqrt( v2sum_integral.evaluate()[0] )
prms = math.sqrt( pressure_int.evaluate()[0] )

temperatureField      = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
temperatureFieldDeriv = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
temperatureField2     = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )

for index, coord in enumerate(mesh.data):
    temperatureField.data[index] = coord[2]
temperatureField2.data[:] = temperatureField.data[:]
temperatureFieldDeriv.data[:] = 0.

kWalls = mesh.specialSets["MinK_VertexSet"] + mesh.specialSets["MaxK_VertexSet"]
advdiffBc = uw.conditions.DirichletCondition( variable        = temperatureField,
                                              indexSetsPerDof = kWalls )

advdiff = uw.systems.AdvectionDiffusion(velocityField=velocityField, phiField=temperatureField, phiDotField=temperatureFieldDeriv, 
                                        fn_diffusivity=1.,conditions=advdiffBc, allow_non_q1=True)
# advdiff2 = uw.systems.AdvectionDiffusion(velocityField=velocityField, phiField=temperatureField, 
#                                        fn_diffusivity=1.,conditions=advdiffBc, allow_non_q1=True, method="SLCN")

# Create a swarm.
if uw.mpi.rank==0: print("Creating Swarm")
swarm = uw.swarm.Swarm( mesh=mesh, particleEscape=True)
# Create a layout object, populate the swarm with particles.
materialIndex = swarm.add_variable('int',1)
swarmLayout = uw.swarm.layouts.PerCellSpaceFillerLayout( swarm=swarm, particlesPerCell=20 )
if uw.mpi.rank==0: print("Populating swarm")
swarm.populate_using_layout( layout=swarmLayout )
# Create a system to advect the swarm
advector = uw.systems.SwarmAdvector( swarm=swarm, velocityField=velocityField, order=2 )

store = vis.Store('{}_RT'.format(jobid),compress=False)
fig = vis.Figure( store, name="firstFig" )
fig.append( vis.objects.Points(swarm, pointSize=2, colourBar=False) )
fig.append( vis.objects.Surface(mesh, pressureField))
fig.append( vis.objects.VectorArrows( mesh, velocityField, scaling=1.0e2))

# update 
# note that the following doesn't make much sense as we're 
# potentially integrating different objects with different time step 
# sizes. however, in particular for the swarm advection, 
# we'd like to use the max timestep as this will stress
# the communication overhead the most.
if uw.mpi.rank==0: print("Doing swarm advection")
dt = 0.1*advector.get_max_dt()
advector.integrate(dt)
if uw.mpi.rank==0: print("Doing advection diffusion")
dt = 0.1*advdiff.get_max_dt()
advdiff.integrate(dt)
# dt = advdiff2.get_max_dt()
# advdiff2.integrate(dt)


# Save things
if do_IO:
    meshFileHandle = mesh.save("{}_Mesh.h5".format(jobid))

    vFH = velocityField.save("{}_velocityField.h5".format(jobid))
    velocityField.xdmf( "{}_velocityField".format(jobid), vFH, "velocity", meshFileHandle, "Mesh" )

    swarmFileHandle = swarm.save("{}_Swarm.h5".format(jobid))
    mH = materialIndex.save("{}_materialIndex.h5".format(jobid))
    materialIndex.xdmf("{}_materialIndex".format(jobid), mH, "material", swarmFileHandle, "Swarm" )

    fig.save()

    # load things
    # first	 create analogues
    mesh_copy = uw.mesh.FeMesh_Cartesian(elementType=els, elementRes=(res,)*dim,minCoord=(20.,)*dim,maxCoord=(30.,)*dim)

    velocityField_copy = uw.mesh.MeshVariable( mesh=mesh_copy, nodeDofCount=dim )

    swarm_copy = uw.swarm.Swarm(mesh = mesh_copy)
    materialIndex_copy = swarm_copy.add_variable( dataType="int", count=1 )

    # now load data and check loaded versions are identical to originals
    mesh_copy.load("{}_Mesh.h5".format(jobid))

    # test
    if not np.allclose(mesh_copy.data, mesh.data):
        raise RuntimeError("Loaded mesh data does not appear to be identical to previous data.")
    velocityField_copy.load("{}_velocityField.h5".format(jobid))
    if not np.allclose(velocityField_copy.data, velocityField.data):
        raise RuntimeError("Loaded velocity data does not appear to be identical to previous data.")


    swarm_copy.load("{}_Swarm.h5".format(jobid))

    if not np.allclose(swarm_copy.particleCoordinates.data, swarm.particleCoordinates.data):
        raise RuntimeError("Loaded swarm data does not appear to be identical to previous data.")
    materialIndex_copy.load("{}_materialIndex.h5".format(jobid))
    if not np.allclose(materialIndex_copy.data, materialIndex.data):
        raise RuntimeError("Loaded material data does not appear to be identical to previous data.")
    # tidy up
    uw.mpi.barrier()
    if uw.mpi.rank==0:
        for filename_endpart in ("RT.gldb","Mesh.h5","velocityField.h5","velocityField.xdmf","Swarm.h5","materialIndex.h5","materialIndex.xdmf"):
            os.remove('{}_{}'.format(jobid,filename_endpart))

if picklename != "None":
    errv = rms_error( velocityField, soln.fn_velocity, mesh, nodal_errors=False)
    errp = rms_error( pressureField, soln.fn_pressure, mesh, nodal_errors=False, normalise=True )

    if uw.mpi.rank==0:
        velocity_key = "Velocity"
        pressure_key = "Pressure"
        try:
            # try and load existing results
            with open(picklename,'rb') as f:
                import pickle
                soln_results = pickle.load(f)
        except:
            # if failed, it's most prob because the file doesn't 
            # exist. in this case, create empty dict.
            soln_results = collections.OrderedDict()
        
        if (soln_name,order,velocity_key) in soln_results:
            err_pre = soln_results[ (soln_name,order,pressure_key) ]
            err_vel = soln_results[ (soln_name,order,velocity_key) ]
        else:
            err_pre = collections.OrderedDict()
            err_vel = collections.OrderedDict()

        err_vel[res] = errv
        err_pre[res] = errp
        soln_results[(soln_name,order,velocity_key)] = err_vel
        soln_results[(soln_name,order,pressure_key)] = err_pre
        with open(picklename,'wb') as f:
            import pickle
            f.write(pickle.dumps(soln_results))



uw.timing.stop()
module_timing_data_orig = uw.timing.get_data(group_by="routine")

# write out data
filename = "Res_{}_Nproc_{}_JobID_{}".format(res,uw.mpi.size,jobid)
import json
if module_timing_data_orig:
    module_timing_data = {}
    for key,val in module_timing_data_orig.items():
        module_timing_data[key[0]] = val
    other_timing["Total_Runtime"] = uw.timing._endtime-uw.timing._starttime
    module_timing_data["Other_timing"] = other_timing
    module_timing_data["Other_data"]   = { "res":res, "nproc":uw.mpi.size, "vrms":vrms, "prms":prms }
    with open(filename+".json", 'w') as fp:
        json.dump(module_timing_data, fp,sort_keys=True, indent=4)

uw.timing.print_table(group_by="routine", output_file=filename+".txt", display_fraction=0.99)
