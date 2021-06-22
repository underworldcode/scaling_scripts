import os
os.environ["UW_TIMING_ENABLE"] = "1"
import underworld3 as uw
uw.tools.parse_cmd_line_options()
import numpy as np
import time
import math
from petsc4py import PETSc
from underworld3.systems import Stokes
from mpi4py import MPI

order = int(os.getenv("UW_ORDER","1"))
res   = int(os.getenv("UW_RESOLUTION",8))
dim   = int(os.getenv("UW_DIM",3))

otol  = float(os.getenv("UW_SOL_TOLERANCE",1.e-6))
max_its  = int(os.getenv("UW_MAX_ITS",-1))

jobid = str(os.getenv("PBS_JOBID",os.getenv("SLURM_JOB_ID","0000000")))

picklename = str(os.getenv("PICKLENAME","None"))

# # Find all available solutions. 
# # Use ordered dict to preserve alphabetical ordering
# import collections
# solns_avail = collections.OrderedDict()
# for _soln in dir(fn.analytic):
#     if _soln[0] == "_": continue  # if private member, ignore
#     # get soln class
#     soln = getattr(fn.analytic,_soln)
#     # check if actually soln
#     if issubclass(soln, fn.analytic._SolBase):
#         solns_avail[_soln] = soln
# soln = solns_avail[soln_name]()


time_post_import   = time.time()
time_launch_mpi    = float(os.getenv("TIME_LAUNCH_MPI"   ,time_post_import))/1000.
time_launch_python = float(os.getenv("TIME_LAUNCH_PYTHON",time_post_import))/1000.

uw.timing.start()

other_timing = {}
other_timing["Python_Import_Time"] = time_post_import - time_launch_python
other_timing[   "MPI_Launch_Time"] = time_launch_python - time_launch_mpi


options = PETSc.Options()

if max_its < 0:
    options["ksp_rtol"] =  otol
else:
    options["ksp_rtol"] =  1e-99
    options["ksp_max_it"] = max_its

boxLength      = 0.9142
ppcell         = 1
amplitude  = 0.02
offset     = 0.2
viscosityRatio = 1.0

# options["help"] = None
# options["pc_type"]  = "svd"
# options["ksp_atol"] =  1.0e-6
options["ksp_monitor"] = None
# options["snes_type"]  = "fas"
options["snes_converged_reason"] = None
options["snes_monitor_short"] = None
# options["snes_view"]=None
# options["snes_test_jacobian"] = None
# options["snes_rtol"] = 1.0e-2  # set this low to force single SNES it. 
options["snes_max_it"] = 1
# options["pc_type"] = "fieldsplit"
# options["pc_fieldsplit_type"] = "schur"
# options["pc_fieldsplit_schur_factorization_type"] = "full"
# # options["fieldsplit_pressure_ksp_rtol"] = 1e-6
# options["fieldsplit_velocity_pc_type"] = "lu"
# options["fieldsplit_pressure_pc_type"] = "jacobi" 
# options["fieldsplit_velocity_ksp_type"] = "gmres"
# sys = PETSc.Sys()
# sys.pushErrorHandler("traceback")

mesh   = uw.mesh.Mesh(elementRes=(res,)*dim,minCoords=(0.,)*dim,maxCoords=(1.,)*dim)
stokes = Stokes(mesh, u_degree=order )
swarm  = uw.swarm.Swarm(mesh)
# Add variable for material
matSwarmVar = swarm.add_variable(name="matSwarmVar",  num_components=1, dtype=PETSc.IntType)
# Note that `ppcell` specifies particles per cell per dim.
swarm.populate(ppcell=ppcell)

import numpy as np
with swarm.access():
    vel_on_particles = uw.function.evaluate(stokes.u.fn,swarm.particle_coordinates.data[0:3])
np.random.seed(0)
with swarm.access(swarm.particle_coordinates):
    factor = 0.5*boxLength/res/ppcell
    swarm.particle_coordinates.data[:] += factor*np.random.rand(*swarm.particle_coordinates.data.shape)

# define these for convenience. 
denseIndex = 0
lightIndex = 1

# material perturbation from van Keken et al. 1997
wavelength = 2.0*boxLength
k = 2. * np.pi / wavelength

# init material variable
with swarm.access(matSwarmVar):
    perturbation = offset + amplitude*np.cos( k*swarm.particle_coordinates.data[:,0] )
    matSwarmVar.data[:,0] = np.where( perturbation>swarm.particle_coordinates.data[:,1], lightIndex, denseIndex )

from sympy import Piecewise, ceiling, Abs

density = Piecewise( ( 0., Abs(matSwarmVar.fn - lightIndex)<0.5 ),
                        ( 1., Abs(matSwarmVar.fn - denseIndex)<0.5 ),
                        ( 0.,                                True ) )

stokes.bodyforce = -density*mesh.N.j

stokes.viscosity = Piecewise( ( viscosityRatio, Abs(matSwarmVar.fn - lightIndex)<0.5 ),
                                (             1., Abs(matSwarmVar.fn - denseIndex)<0.5 ),
                                (             1.,                                True ) )

# note with petsc we always need to provide a vector of correct cardinality. 
bnds = mesh.boundary
stokes.add_dirichlet_bc( (0.,0.), [bnds.TOP,  bnds.BOTTOM], (0,1) )  # top/bottom: function, boundaries, components 
stokes.add_dirichlet_bc( (0.,0.), [bnds.LEFT, bnds.RIGHT ], 0  )  # left/right: function, boundaries, components

volume_int = uw.maths.Integral( mesh, 1. )
volume = volume_int.evaluate()
v_dot_v_int = uw.maths.Integral(mesh, stokes.u.fn.dot(stokes.u.fn))

import math


# Solve time
stokes.solve()
vrms = math.sqrt(v_dot_v_int.evaluate()/volume)

dt = stokes.dt()
with swarm.access():
    vel_on_particles = uw.function.evaluate(stokes.u.fn,swarm.particle_coordinates.data)

with swarm.access(swarm.particle_coordinates):
    swarm.particle_coordinates.data[:]+=dt*vel_on_particles

if MPI.COMM_WORLD.rank==0: print(f"VRMS = {vrms}")


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
    module_timing_data["Other_data"]   = { "res":res, "nproc":uw.mpi.size, "vrms":vrms,  }
    with open(filename+".json", 'w') as fp:
        json.dump(module_timing_data, fp,sort_keys=True, indent=4)

uw.timing.print_table(group_by="routine", output_file=filename+".txt", display_fraction=0.99)
