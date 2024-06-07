#!/usr/bin/python3
"""
This creates 

"""
from esys.escript import *
import argparse
import importlib, sys, os, json
#from fingal import *
import subprocess
from esys.finley import ReadGmsh
from tools import GravityData, Surface, Moho, LAB, WGS84, StereographicProjection, TopCentricTransformation, ECEFTrafo
sys.path.append(os.getcwd())
from esys.weipa import saveSilo
from esys.weipa import saveVTK


sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='Creates a mesh fly file. the gmsh mesh generator is used.', epilog="fingal by l.gross@uq.edu.auversion 21/12/2020")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--plot', '-p',  dest='plot', action='store_true', default=False, help="plot input data (in stage 1 only)")
parser.add_argument('--silo', '-s',  dest='silo', default=None, help="name of silo file to visualize mesh and subdomain tags (in stage 3 only!)")
parser.add_argument('--usePlane', '-u',  dest='usePlane', action='store_true', default=False, help="no stereographic project is applied. (in stage 2 only!)")
parser.add_argument('--stage', '-S',  dest='stage', type=int, default=None, help="stage to perform =1-create theSurface mesh, =2 - add theSurface to 2D mesh, =3 convert to fly; if None all steps are performed including mesh generation.")
parser.add_argument('--tmpfile', '-t',  dest='tmpfile', type=str, default="tmp", help="prefix for temporary files.")

args = parser.parse_args()

config = importlib.import_module(args.config)
print("Configuration "+args.config+".py imported.")

#========================================================

GEOFN=args.tmpfile
MESHFN=args.tmpfile
FLYFN=""

earthModel=WGS84()



DataHeightAboveGround = config.Mesh_DataHeightAboveGround
DataMeshSizeVertical = config.Mesh_DataMeshSizeVertical
CoreThickness= config.Mesh_CoreThickness
AirLayerThickness=config.Mesh_AirLayerThickness
PaddingX, PaddingY, PaddingZ, PaddingAir=config.Mesh_PaddingX, config.Mesh_PaddingY, config.Mesh_PaddingZ, config.Mesh_PaddingAir
MeshSizeAirFactor=config.Mesh_MeshSizeAirFactor
MeshSizeCoreFactor=config.Mesh_MeshSizeCoreFactor
MeshSizePaddingFactor=config.Mesh_MeshSizePaddingFactor
#config = importlib.import_module(args.config)
#print("configuration "+args.config+" imported.")



GEO2DTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "AboveGroundTemplate2D.geo")
GEO3DTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "AboveGroundTemplate3D.geo")
GEOFN1=GEOFN+"_2D.geo"
MSHN1=MESHFN+"_2D.msh"
MSHN2=MESHFN+"_2D_topo.msh"
GEOFN3=GEOFN+"_3D.geo"
MSHN3=MESHFN+".msh"


# read data files:
gravData=GravityData(config.data_name, isCellCentered=True)
print(len(gravData), " gravity data read from ", gravData.FileName)
print("gravity data grid = ", gravData.getGridDimensions(), gravData.getRange())

theSurface=Surface(config.surface_name, isCellCentered=True)
print(len(theSurface), " Surface data read from ", theSurface.FileName)
print("surfaces data grid = ", theSurface.getGridDimensions(), theSurface.getRange())

theMoho=Moho(config.moho_name, isCellCentered=True)
print(len(theMoho), " Moho data read from ", theMoho.FileName)
print("moho data grid = ", theMoho.getGridDimensions(), theMoho.getRange())

theLAB=LAB(config.lab_name, isCellCentered=True)
print(len(theLAB), " LAB data read from ", theLAB.FileName)
print("moho data grid = ", theLAB.getGridDimensions(), theLAB.getRange())

#assert gravData.rangeX[0] >= theSurface.rangeX[0]
#assert gravData.rangeX[1] <= theSurface.rangeX[1]
#assert gravData.rangeY[0] >= theSurface.rangeY[0]
#assert gravData.rangeY[1] >= theSurface.rangeY[1]

# transformations:
projection=StereographicProjection(config.Latitude0, model=earthModel)
cartTrafo=TopCentricTransformation(config.Longitude0, config.Latitude0, trafo=ECEFTrafo(model=earthModel))



SurfaceLevel=theSurface.getTopographyData().mean()
Topomin, Topomax=theSurface.getTopographyData().min(), theSurface.getTopographyData().max()
MohoLevel=theMoho.getMohoData().mean()
LABLevel=theLAB.getLABData().mean()

DataLevel=gravData.getHeightData().mean()

print("SurfaceLevel [km] = ", SurfaceLevel/1000.)
print("Surface  range [km]=", Topomin/1000, Topomax/1000)
print("MohoLevel  [km] = ", MohoLevel/1000.)
print("LABLevel  [km] = ", LABLevel/1000.)
print("DataLevel  [km] = ", DataLevel/1000.)

DataSpacingX, DataSpacingY = gravData.getGridSpacing()
Xmin, Xmax=gravData.getRange()[0]
Ymin, Ymax=gravData.getRange()[1]
DataThickness=min(DataSpacingX, DataSpacingY)
DataHeightAboveGround=DataLevel-DataThickness/2 #???
DataMeshSizeVertical=DataThickness

print("Data X range = ", Xmin, Xmax,' spacing = ', DataSpacingX)
print("Data Y range = ", Ymin, Ymax,' spacing = ', DataSpacingY)

print("DataThickness = ", DataThickness)
print("DataHeightAboveGround = ", DataHeightAboveGround)
print("Thickness crust = ", SurfaceLevel-MohoLevel)

projection=StereographicProjection(trueScaleLatitude=config.TrueScaleLatitude, hemisphere=config.Hemisphere, model=earthModel)
cartTrafo=TopCentricTransformation(config.Longitude0, config.Latitude0, trafo=ECEFTrafo(model=earthModel))

if DataHeightAboveGround < DataMeshSizeVertical:
    CorrectionDataHeight= DataMeshSizeVertical-DataHeightAboveGround
    DataHeightAboveGround+=CorrectionDataHeight  #This should be the cell center ?
    print(">>>>> CorrectionDataHeight is  ", CorrectionDataHeight)    
    print(">>>>> DataHeightAboveGround is increased to ", DataHeightAboveGround)
else:
    CorrectionDataHeight=0.

AirLayerThickness+=CorrectionDataHeight   

if AirLayerThickness < DataHeightAboveGround *3 :
    AirLayerThickness=DataHeightAboveGround *3
    print(">>>>> AirLayerThickness is increased to ", AirLayerThickness)
    
if SurfaceLevel-MohoLevel < DataMeshSizeVertical :
    CorrectionMohoLevel=SurfaceLevel-2*DataMeshSizeVertical-MohoLevel
    print(">>>>> MohoLevel is ", MohoLevel)
    MohoLevel+=CorrectionMohoLevel
    print(">>>>> CorrectionMohoLevel is ", CorrectionMohoLevel)
    print(">>>>> MohoLevel is reduced to ", MohoLevel)
else:
    CorrectionMohoLevel=0

if MohoLevel-LABLevel < DataMeshSizeVertical :
    CorrectionLABLevel=MohoLevel-2*DataMeshSizeVertical-LABLevel
    print(">>>>> LABLevel is ", LABLevel)
    LABLevel+=CorrectionLABLevel
    print(">>>>> CorrectionLABLevel is ", CorrectionLABLevel)
    print(">>>>> LABLevel is reduced to ", LABLevel)
else:
    CorrectionLABLevel=0

    
if MohoLevel+CoreThickness < 2*DataMeshSizeVertical:
    CoreThickness=3*DataMeshSizeVertical-MohoLevel
    print(">>>>> CoreThickness is increased to ", CoreThickness)    
    

mapping= { 
"SurfaceLevel" : SurfaceLevel,
"MohoLevel": MohoLevel,
"LABLevel": LABLevel,
"DataThickness" : DataThickness,
"DataRefX" : Xmin, 
"DataRefY" : Ymin,
"DataHeightAboveGround" : DataHeightAboveGround,
"DataSpacingX" : DataSpacingX,
"DataSpacingY" : DataSpacingY,
"DataNumX" : gravData.getGridDimensions()[1],
"DataNumY" : gravData.getGridDimensions()[0],
"DataMeshSizeVertical" : DataMeshSizeVertical,
"CoreThickness" : CoreThickness,
"AirLayerThickness" : AirLayerThickness,
"PaddingX" : PaddingX,
"PaddingY" : PaddingY,
"PaddingZ" : PaddingZ,
"PaddingAir" : PaddingAir,
"MeshSizeAirFactor" : MeshSizeAirFactor,
"MeshSizeCoreFactor" : MeshSizeCoreFactor,
"MeshSizePaddingFactor" : MeshSizePaddingFactor
}

#this generates the 2D mesh:
if args.stage ==1 or args.stage == None:

    # lets plot the data:
    if args.plot:
        import matplotlib.pyplot as plt
        plt.tricontourf(theSurface.X.reshape((-1,)), theSurface.Y.reshape((-1,)), theSurface.getTopographyData().reshape((-1,)), levels=14, cmap="RdBu_r")
        plt.title("surface topography")
        plt.colorbar()
        plt.savefig("surface.png")
        plt.clf()
        plt.tricontourf(theMoho.X.reshape((-1,)), theMoho.Y.reshape((-1,)), theMoho.getMohoData().reshape((-1,)), levels=14, cmap="RdBu_r")
        plt.title("Moho")
        plt.colorbar()
        plt.savefig("moho.png")
        plt.clf()
        plt.tricontourf(gravData.X.reshape((-1,)), gravData.Y.reshape((-1,)), gravData.getGravityData().reshape((-1,)), levels=14, cmap="RdBu_r")
        plt.title("Gravity")
        plt.colorbar()
        plt.savefig("grav.png")
        print("data have been plotted.")
    
    print("Setting: ")
    print(json.dumps(mapping, sort_keys=True, indent=4)[2:-2])
    text=open(GEO2DTEMPLATE,'r').read().format(**mapping)
    open(GEOFN1, "w").write(text)
    print("GMSH geofile has been written to ",GEOFN1)
    if args.stage ==1 :
        print("to generate 2D mesh run:")
        print("     gmsh -2 -format msh2 -o %s %s"%(MSHN1, GEOFN1))
    else:
        rp=subprocess.run(["gmsh", "-2", "-format", "msh2", "-algo", "auto", "-o", MSHN1, GEOFN1])
        print(rp.returncode)
        rp.check_returncode()
        print(">> 2D GMSH mesh file %s generated."%MSHN1)

#now the mesh is deformed according to theSurface infos:
if args.stage == 2 or args.stage == None:
    
    projection=StereographicProjection(trueScaleLatitude=config.TrueScaleLatitude, hemisphere=config.Hemisphere, model=earthModel)
    cartTrafo=TopCentricTransformation(config.Longitude0, config.Latitude0, trafo=ECEFTrafo(model=earthModel))

    Bottom=-(CoreThickness+PaddingZ)
    Top=AirLayerThickness+PaddingAir
    Depth=-(CoreThickness+PaddingZ)
    print("read 2D mesh from ", MSHN1)
    # apply topography:
    fin=open(MSHN1, 'r')
    fout=open(MSHN2,'w')
    line=fin.readline()
    while line:
        if line.startswith("$NOD") or line.startswith("$Nodes"):
                process_nodes=True
                numNodes=int(fin.readline())
                if line.startswith("$NOD"):
                    fout.write("$NOD\n")
                else:
                    fout.write("$Nodes\n")
                fout.write("%d\n"%numNodes)
                print(numNodes," nodes found.")
                cc=0
                while cc < numNodes:
                    line=fin.readline()
                    i,x,y,z=line.split(" ")
                    i, x,y,z =int(i), float(x), float(y), float(z) 
                    lab=theLAB.getLAB(x,y)+CorrectionLABLevel
                    moho=theMoho.getMoho(x,y)+CorrectionMohoLevel
                    topo=theSurface.getTopography(x,y)
                    height=gravData.getHeight(x,y)-DataLevel+DataHeightAboveGround
                    if z < -CoreThickness:
                        # nothing to do!
                        h=z
                        k=1
                    elif z < LABLevel: # between CoreThickness and Moholevel -> fit between unchanged Depth and (new) Moho 
                        h=(lab+CoreThickness)/(LABLevel+CoreThickness)*(z+CoreThickness)-CoreThickness
                        k=1
                    elif z < MohoLevel: # between CoreThickness and Moholevel -> fit between unchanged Depth and (new) Moho 
                        h=(moho-lab)/(MohoLevel-LABLevel)*(z-LABLevel)+lab
                        k=2
                    elif z < SurfaceLevel: # between Moholevel and SurfaceLevel -> fit between theMoho and (new) Topography 
                        h=(topo-moho)/(SurfaceLevel-MohoLevel)*(z-MohoLevel)+moho
                        k=3
                    elif z < DataHeightAboveGround: # between SurfaceLevel and DataHeightAboveGround -> fit between topo and (new) height 
                        h=(height-topo)/(DataHeightAboveGround-SurfaceLevel)*(z-SurfaceLevel)+topo
                        k=4
                    elif z < DataHeightAboveGround+DataThickness: # between DataHeightAboveGround and DataHeightAboveGround+DataThickness -> fit between height  and height + DataThickness
                        h=z-DataHeightAboveGround+height
                        k=5
                    elif z < AirLayerThickness: # between DataHeightAboveGround+DataThickness and AirLayerThickness -> fit between height + DataThickness and AirLayerThickness
                        h=(AirLayerThickness-height - DataThickness)/(AirLayerThickness-DataHeightAboveGround-DataThickness)*(z-DataHeightAboveGround-DataThickness)+height+DataThickness
                        k=6
                    else:
                        h=z
                        k=7
                    if args.usePlane:
                        X,Y,Z=x,y,h
                    else:
                        p=projection.toPosition(x,y)
                        X,Y,Z =cartTrafo.getPosition2ENU(p[0], p[1],h) 
                    fout.write("%d %e %e %e\n"%(i, X,Y,Z))
                    cc+=1
        else:
                fout.write(line)
        line=fin.readline()

    fin.close()
    fout.close()
    print("updated 2D mesh written to ",MSHN2)
    text=open(GEO3DTEMPLATE,'r').read().format(MESHFILE=MSHN2)
    open(GEOFN3, "w").write(text)
    print("3D gmsg geo file written to ",GEOFN3)
    if args.stage ==2 :
        print("to generate 3D mesh run:")
        print("     gmsh -3 -format msh2 -algo frontal -o %s %s"%(MSHN3, GEOFN3))
    else:
        # now we are ready to generate the 3D mesh:  
        rp=subprocess.run(["gmsh", "-3",  "-algo", "frontal", "-o", MSHN3, GEOFN3])

        rp.check_returncode()
        print(">> GMSH mesh file %s generated."%MSHN3)

    print("*** STAGE 2 COMPLETED ***")

    print("*** STAGE 3 COMPLETED ***")