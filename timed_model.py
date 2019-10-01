import os
os.environ["UW_ENABLE_TIMING"] = "1"
import underworld as uw
from underworld import function as fn
import underworld.visualisation as vis
import numpy as np
import time

order = int(os.getenv("UW_ORDER","2"))
res   = int(os.getenv("UW_RESOLUTION",16))
dim   = int(os.getenv("UW_DIM",3))

itol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-10))
otol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-10))

soln_name = str(os.getenv("UW_MODEL","SolDB3d"))

do_IO  = bool(int(os.getenv("UW_ENABLE_IO","0")))

jobid = str(os.getenv("UW_JOB_ID","00000"))
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
time_launch_srun   = float(os.getenv("TIME_LAUNCH_SRUN"  ,time_post_import))/1000.
time_launch_python = float(os.getenv("TIME_LAUNCH_PYTHON",time_post_import))/1000.

uw.timing.start()

other_timing = {}
other_timing["Python_Import_Time"] = time_post_import - time_launch_python
other_timing["Container_Launch_Time"] = time_launch_python - time_launch_srun

def normalise_press(press):
    intSwarm = uw.swarm.GaussIntegrationSwarm(mesh,3)  # use 3 point gauss swarms for efficiency
    av_press = uw.utils.Integral( press, mesh, integrationSwarm=intSwarm, integrationType=None).evaluate()[0]
    
    return press - av_press

def rms_error(numeric, analytic, mesh):
    '''
    Calculates the rms error.
    
    Returns
    -------
    abs, abs_scaled: float
        The absolute and scaled absolute errors.
    '''

    delta     = analytic - numeric
    delta_dot = fn.math.dot(delta,delta)

    analytic_dot = fn.math.dot(analytic,analytic)
    
    # l2 norms
    intSwarm = uw.swarm.GaussIntegrationSwarm(mesh,3)  # use 3 point gauss swarms for efficiency
    rms_err_abs = np.sqrt(uw.utils.Integral(    delta_dot, mesh, integrationSwarm=intSwarm, integrationType=None ).evaluate()[0])
#     rms_sol_ana = np.sqrt(uw.utils.Integral( analytic_dot, mesh, integrationSwarm=intSwarm, integrationType=None ).evaluate()[0])
#     rms_err_sca = rms_err_abs / rms_sol_ana
        
    return rms_err_abs, 0. #rms_err_sca


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
solver.set_inner_rtol(itol)
solver.set_outer_rtol(otol)
if soln.nonlinear==True:
    stokes.fn_viscosity = 1.
    solver.solve()
    stokes.fn_viscosity = visc

# functions for calculating RMS velocity
vdotv = fn.math.dot(velocityField,velocityField)
v2sum_integral  = uw.utils.Integral( mesh=mesh, fn=vdotv )
volume_integral = uw.utils.Integral( mesh=mesh, fn=1. )

pdotp = fn.math.dot(pressureField,pressureField)
pressure_int = uw.utils.Integral( mesh=mesh, fn=pdotp )

# Get instantaneous Stokes solution
solver.solve()
# Calculate the RMS velocity.
vrms = math.sqrt( v2sum_integral.evaluate()[0] )
prms = math.sqrt( pressure_int.evaluate()[0] )

stats=solver.get_stats()
solver.print_stats()

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

# Create a dummy temperature field.
advdiff = uw.systems.AdvectionDiffusion(velocityField=velocityField, phiField=temperatureField, phiDotField=temperatureFieldDeriv, 
                                        fn_diffusivity=1.,conditions=advdiffBc, allow_non_q1=True)
advdiff2 = uw.systems.AdvectionDiffusion(velocityField=velocityField, phiField=temperatureField, 
                                        fn_diffusivity=1.,conditions=advdiffBc, allow_non_q1=True, method="SLCN")

# Create a swarm.
swarm = uw.swarm.Swarm( mesh=mesh, particleEscape=True)
# Create a layout object, populate the swarm with particles.
materialIndex = swarm.add_variable('int',1)
swarmLayout = uw.swarm.layouts.PerCellSpaceFillerLayout( swarm=swarm, particlesPerCell=40 )
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
dt = advector.get_max_dt()
advector.integrate(dt)
dt = advdiff.get_max_dt()
advdiff.integrate(dt)
dt = advdiff2.get_max_dt()
advdiff2.integrate(dt)


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

uw.timing.stop()
module_timing_data_orig = uw.timing.get_data(group_by="line_routine")

# write out data
filename = "{}_Res_{}_Nproc_{}_JobID_{}".format(os.getenv("UW_JOB_NAME","Job"),res,uw.mpi.size,jobid)
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

uw.timing.print_table(group_by="line_routine", output_file=filename+".txt", display_fraction=0.99)

if picklename != "None":
    velocity_key = "Velocity"
    pressure_key = "Pressure"
    if uw.mpi.rank==0:
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

    # grab copy of analytic solutions onto mesh variables. 
    # this is to avoid excessive calls into analytic solutions,
    # some of which (solH for example) are prohibitively 
    # expensive to calculate 
    #velocityFieldA = mesh.add_variable(mesh.dim)
    #pressureFieldA = mesh.subMesh.add_variable(1)
    #velocityFieldA.data[:] = soln.fn_velocity.evaluate(mesh)
    #pressureFieldA.data[:] = soln.fn_pressure.evaluate(mesh.subMesh)
    pressn = normalise_press(pressureField)
    pressa = normalise_press(soln.fn_pressure)
    #pressa = normalise_press(pressureFieldA)

    err_vel[res] = rms_error( velocityField, soln.fn_velocity, mesh )
    #err_vel[res] = rms_error( velocityField, velocityFieldA, mesh )
    err_pre[res] = rms_error( pressn,        pressa,         mesh )


    # record full state back to pickled dict
    if uw.mpi.rank==0 :
        soln_results[(soln_name,order,velocity_key)] = err_vel
        soln_results[(soln_name,order,pressure_key)] = err_pre
        with open(picklename,'wb') as f:
            import pickle
            f.write(pickle.dumps(soln_results))


