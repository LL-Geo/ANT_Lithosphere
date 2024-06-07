#!/usr/bin/python3
__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse

from esys.escript import *
from esys.finley import ReadMesh
from esys.finley import ReadGmsh
import esys.escript.unitsSI as U
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import saveVTK

#from scipy.io import netcdf_file
from esys.escript.pdetools import MaskFromBoundaryTag
#from converters import writeNetCDF, grepValuesByMask
from esys.escript.pdetools import Locator
#import matplotlib.pyplot as plt
import argparse
from tools import *
from tools2 import *

################################################################
#              CONFIG FILE 
################################################################
parser = argparse.ArgumentParser(description='Gravity forward. ENU coordinate system.', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
args = parser.parse_args()
config = importlib.import_module(args.config)
print("Configuration "+args.config+".py imported.")

################################################################
#              MESH
################################################################
domain=ReadGmsh(config.mesh_name, 3, optimize=True)
print("mesh read from ",config.mesh_name)

################################################################
#              model basics from config file
################################################################
earthModel= WGS84()
Latitude0  = config.Latitude0
Longitude0 = config.Longitude0
TrueScaleLatitude = config.TrueScaleLatitude
Hemisphere = config.Hemisphere
data_scale = config.data_scale

################################################################
#              TRANSFORMATIONS
################################################################
#
# transforming between ECEF (geocentric) and geographical
#                     (X,Y,Z)     <==>       (lon,lat, h)
transform1 = ECEFTrafo(model=earthModel) 

# transforming between geocentric and ENU (and geographical)
#                       (X, Y, Z) <==>  (U, V, W)
transform2 = TopCentricTransformation(Longitude0, Latitude0, trafo = transform1)

# stereographic projection (for gravity data)
#
projection = StereographicProjection(trueScaleLatitude=TrueScaleLatitude, hemisphere=Hemisphere, model=earthModel)

################################################################
#              COORDINATES FOR MESH CENTROIDS
################################################################
#
#  Mesh         ENU    (UVW)
#  geocentric   ECEF   (XYZ)
#  geographical        (LOLAH)
#  stereographic       (x,y)          
#
# get the coordinate of the centroid of each element in ENU (UVW)
# mesh coordinate system
UVW = ReducedFunction(domain).getX()

#   ENU   =>   ECEF
# (UVW) => (XYZ)
# 
XYZ = transform2.getENU2ECEF(UVW[0],UVW[1],UVW[2])

#  ECEF   => geographical 
#  (XYZ)  => (LOLAH)
LOLAH = transform1.getPosition(XYZ[0],XYZ[1],XYZ[2])

#  geographical to stereographic
#  (LOLAH) => (x,y)
x,y = projection.toXY(LOLAH[0],LOLAH[1])
centrePts = Data(0,(3,),ReducedFunction(domain))
centrePts[0] = x 
centrePts[1] = y
centrePts[2] = LOLAH[2]

################################################################
#              Reference Density
################################################################
# crust density
theCrust = DensityModel(config.dens_name)
print("the crust density is read from file ",theCrust.FileName)
# topography
theSurface=Surface(config.surface_name, isCellCentered=True)
print("the surface is read from file ",theSurface.FileName)
# moho
theMoho=Moho(config.moho_name, isCellCentered=True)
print("the moho read from file ",theMoho.FileName)
# mantle
theMantle=Mantle(config.mantle_name)
print("the Mantle read from file ",theMantle.FileName)
#LAB
theLAB=LAB(config.lab_name, isCellCentered=True)
print(len(theLAB), " LAB data read from ", theLAB.FileName)
print("moho data grid = ", theLAB.getGridDimensions(), theLAB.getRange())
#====
# make surface density 
mask=Scalar(0,FunctionOnBoundary(domain))
mask.setTaggedValue("EarthSurface", 1)
#print(mask)
#====
UVW_face=FunctionOnBoundary(domain).getX()
XYZ_face=transform2.getENU2ECEF(UVW_face[0],UVW_face[1],UVW_face[2])
LOLAH_face = transform1.getPosition(XYZ_face[0],XYZ_face[1],XYZ_face[2])
xy_face= projection.toXY(LOLAH_face[0],LOLAH_face[1])
#====
rho_surf=theSurface.interpolate2DCellsToEscript(xy_face, target=Scalar(0., FunctionOnBoundary(domain)), mask=mask, data='Density_Surface')
rho_surf = rho_surf*config.rho_surf_scale
print("surface density  =",rho_surf)
del UVW_face, XYZ_face, LOLAH_face, xy_face



################################################################
#              DOWN IN MESH COORDINATE SYSTEM
################################################################
# defines point on bottom surface directly below each element centroid.
# should be no division by zero because no element centroid will be on the surface.
hmin=inf(LOLAH[2])
XYZ0 = transform1.getXYZ(LOLAH[0], LOLAH[1], Scalar(hmin, ReducedFunction(domain)))
UVW0 = transform2.getECEF2ENU(XYZ0[0], XYZ0[1], XYZ0[2])

gdir = Vector(0., ReducedFunction(domain))
gdir[0] = (UVW[0]-UVW0[0])
gdir[1] = (UVW[1]-UVW0[1])
gdir[2] = (UVW[2]-UVW0[2])
lengthgdir = length(gdir)
gdir[0] = safeDiv(gdir[0],lengthgdir) 
gdir[1] = safeDiv(gdir[1],lengthgdir) 
gdir[2] = safeDiv(gdir[2],lengthgdir)

################################################################
#              DOMAIN REGIONS 
################################################################
# ground
ground_e =Scalar(0,ReducedFunction(domain))
for dirt in config.subsurfaceTags:
    ground_e.setTaggedValue(dirt, 1.)
ground_e.expand()

# ground
Inv_e2 =Scalar(0,ReducedFunction(domain))
for dirt in config.inversionTags:
    Inv_e2.setTaggedValue(dirt, 1.)
Inv_e2.expand()

# data area
dataArea = Scalar(0,ReducedFunction(domain))        
dataArea.setTaggedValue(config.DataAreaTag, 1.)
dataArea.expand() 


################################################################
#              GRAVITY DATA 
################################################################
# read in data for interpolation  
# dataArea element centres
gravData=GravityData(config.data_name, isCellCentered=True)
print(len(gravData), " gravity data read from ", gravData.FileName)
print("gravity data grid = ", gravData.getGridDimensions(), gravData.getRange())
g_data = Scalar(0., ReducedFunction(domain))
gravData.interpolateGravityToEscript(centrePts,target= g_data, mask = dataArea)
g_data = g_data * config.data_scale
print("gravity data range ", inf(g_data), " to ", sup(g_data))
#saveVTK("datastuff", g_data=g_data,x=x,y=y)


################################################################
################################################################
#              INVERSION 
################################################################
################################################################


model = curvedGravityModel(domain, fixBase=True)


g_ref_dataArea=load("g_ref_dataArea.nc",domain)

Inv_e=load("Inv_e.nc",domain)
rho_0=load("rho_0.nc",domain)





if config.Option1:
    G_DATA = g_data - g_ref_dataArea
    G_DATA_mean = integrate(G_DATA)/integrate(dataArea) 

    G_DATA = dataArea*(G_DATA - G_DATA_mean)
    print("new data gravity range", inf(G_DATA), sup(G_DATA))
    m0 = Scalar(0., ContinuousFunction(domain))

    grav = ACGravity(domain, dataArea = dataArea, g_dir=gdir, g_data = G_DATA, rho_0 = Inv_e,
                 ground=Inv_e, mu=config.mu, m0=m0, fixBase=config.fixBase, atol=config.atol, rtol=config.rtol, 
                 iter_max=config.iter_max, pde_tol=config.pdetol, 
                 output_name=config.output_name+'1', verboseLevel=config.VerboseLevel)

    m, rho_anomaly = grav.solve()
    model.setDensity(rho_anomaly)
    model.setSurfDensity(0.0)
    potential_anomaly = model.getGravityPotential()
    g_anom = model.getDownGravity(gdir)
    g_anom_dataArea = model.getDownGravityDataArea(dataArea, gdir)
    rho_final=rho_0+rho_anomaly
    rho_anomaly.dump("rho_anomaly.nc")

    rho_final.dump("rho_final.nc")
    saveDataCSV("Density_a.csv", mask=Inv_e2, X=x,Y=y,Z=LOLAH[2],
           d_anmo=rho_anomaly)
    saveDataCSV("Density_0.csv", mask=Inv_e2, X=x,Y=y,Z=LOLAH[2],
           d_all = rho_0)
    saveDataCSV("Density_all.csv", mask=Inv_e2, X=x,Y=y,Z=LOLAH[2],
           d_all = rho_final)

    saveVTK("Density", X=x,Y=y,Z=LOLAH[2],
           d_anmo=rho_anomaly,d_all=rho_final)
    #rho_anomaly.dump(args.config+"rho_ano.nc")
    model.setDensity(rho_final)
    model.setSurfDensity(rho_surf)

    potential_final = model.getGravityPotential()
    g_final = model.getDownGravity(gdir)
    g_final_dataArea = model.getDownGravityDataArea(dataArea, gdir)
    saveDataCSV("gravity_5k.csv", mask=dataArea, X=x,Y=y,Z=LOLAH[2],
           g_data=g_data, g_data_inv = G_DATA, g_comp_anomaly =g_anom,
           g_final = g_final)

    print(" --------------------------------------------------------------------")
    print(" --------------------------------------------------------------------")
