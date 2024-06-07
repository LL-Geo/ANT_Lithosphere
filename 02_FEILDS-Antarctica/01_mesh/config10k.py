#
# Configuration file for gravity inversion for use by planeGravInv.py 

# File names
# mesh file name.  This needs to be in msh or fly format
mesh_name = "tmp.msh"

# data file name:
# See readme.md for more details.
data_name = "./Model/FAA_10k.csv"

# data file name for density
dens_name = "./Model/Crust_density.csv"

# data file name for topograpy surface
# format: x,y,h, top crustal density, surface density
surface_name = "./Model/Surface_10k.csv"

# data file name for moho surface
# format: x,y,h, bottom crustal density, top mantel density
moho_name = "./Model/Moho_10k.csv"


lab_name = "./Model/LAB.csv"

mantle_name= "./Model/Mantle_Density_PUM.csv"


# Geometrical set-up as not pulled from input file
Mesh_DataHeightAboveGround = 10000.
Mesh_DataMeshSizeVertical = 10000.
Mesh_CoreThickness= 400000.
Mesh_AirLayerThickness=2000.
Mesh_PaddingX, Mesh_PaddingY, Mesh_PaddingZ, Mesh_PaddingAir=100000., 100000., 100000., 100000.
Mesh_MeshSizeAirFactor=1
Mesh_MeshSizeCoreFactor=1
Mesh_MeshSizePaddingFactor=5
#
#
#
#
Latitude0 = -90.
Longitude0 = 0.
TrueScaleLatitude = 71.
Hemisphere = -1

subsurfaceTags = ["Crust", "CrustPadding", "Mantel", "MantelPadding"]
DataAreaTag = "DataArea"
fixBase = True

# Inversion constants:
#
# scale between misfit and regularization
mu       = 1.e-14
#
# used to scale computed density.  kg/m^3    
rho_0    = 1.
#
# IPCG tolerance *|r| <= atol+rtol*|r0|*  (energy norm)
# absolute tolerance for IPCG interations
atol     = 0.        
#
# relative tolerance for IPCG iterations   
rtol     = 1.e-2
#
# tolerance for solving PDEs 
# make sure this is not more than the square of rtol     
pdetol   = 1.e-10    
#
# maximum number of IPCG iterations
iter_max = 500 
#
# data scale.  Program assumes m/s^2, 
# converts micrometres/s^2 to m/s^2
data_scale = 1.e-5
#
#
density_scale =1.e3
#
rho_surf_scale=1.e3

Option1= True

output_name="Model_20k_mu_{0:1.3e}_".format(mu)

# Level for the verbosity of the output, "low", "medium" or "high".
# low: 
#   screen outputs:
#      data range, 
#      summaries of gravity data and final gravity
#      initial, final and difference misfits
#   file output:
#      silo of final solution
# medium: low outputs + 
#   screen outputs:
#      residual norm from the IPCG iterations
# high: medium outputs + 
#   screen outputs:
#      misfit and smoothing value at each iteration step
#   file outputs:
#      csv files for misfit and smoothing at each IPCG iteration
#      silos at misfit values of 0.05, 0.01, 0.008 and 0.005. (Initial misfit is 0.5.)

#VerboseLevel = "low"
#VerboseLevel = "medium"
VerboseLevel = "high"
