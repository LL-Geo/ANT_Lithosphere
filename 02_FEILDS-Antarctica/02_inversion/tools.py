import numpy as np
from esys.escript import *
from scipy.interpolate import RegularGridInterpolator
# update 24/04/2024 density interpolation bug due to z is negtive 
def atan2fix(arg0, arg1): # y,x
   """
   Returns inverse tangent of argument ``arg0`` over ``arg1`` 
   """
   z1=whereZero(arg1, rtol=EPSILON)
   #m2=whereNegative(arg1*arg0)
   #print("atan2fix =",m)
   #s=(atan(arg0/(arg1+z1)) + (wherePositive(arg0)-whereNegative(arg0))*whereNegative(arg1)*np.pi )*(1-z1)+(wherePositive(arg0)-whereNegative(arg0))*whereNegative(arg1)*np.pi/2*z1
   s=(atan(arg0/(arg1+z1)) + sign(arg0)*whereNegative(arg1)*np.pi )*(1-z1)+sign(arg0)*whereNegative(arg1)*np.pi/2*z1
   return s

class WGS84(object):
    EquatorialRadius =6378137.0
    f=1./298.257223563
    def __init__(self):
        self.PolarRadius=self.EquatorialRadius *(1-self.f)
        self.e2=(self.EquatorialRadius **2-self.PolarRadius**2)/self.EquatorialRadius **2
        self.e=sqrt(self.e2)
        self.ee2=(self.EquatorialRadius **2-self.PolarRadius**2)/self.PolarRadius**2 
        
class ECEFTrafo(object):

    def __init__(self, model=WGS84()):
        self.model=model
        
    def getXYZ(self, lon, lat,  h):
        p=np.deg2rad(1)*lat
        l=np.deg2rad(1)*lon
        v=self.model.EquatorialRadius /sqrt(1-self.model.e2*sin(p)**2)
        X=(v+h)*cos(p)*cos(l)
        Y=(v+h)*cos(p)*sin(l)
        Z=(v*(1-self.model.e2)+h)*sin(p)
        return X,Y,Z

    def getPosition(self, X, Y, Z):        
        P=sqrt(X**2+Y**2)
        t=atan(Z*self.model.EquatorialRadius /(P*self.model.PolarRadius))
        p=atan((Z+self.model.ee2*self.model.PolarRadius*sin(t)**3)/(P-self.model.e2*self.model.EquatorialRadius*cos(t)**3))
        v=self.model.EquatorialRadius /sqrt(1-self.model.e2*sin(p)**2)
        l=atan2fix(Y,X)#l=atan(Y/X)
        h=P/cos(p)-v
        return np.rad2deg(1)*l, np.rad2deg(1)*p, h
    

    
class TopCentricTransformation(object):
    def __init__(self, Longitude0=0, Latitude0=0, trafo=ECEFTrafo()):
        self.ECEF=trafo
        self.Latitude0=Latitude0
        self.Longitude0=Longitude0
        self.origin=self.ECEF.getXYZ(Longitude0, Latitude0,0)
        self.Latitude0=np.deg2rad(Latitude0)
        self.Longitude0=np.deg2rad(Longitude0)
        self.matrix= np.array([ [-sin(self.Longitude0)               ,  cos(self.Longitude0)               , 0.            ],
                              [-sin(self.Latitude0)*cos(self.Longitude0), -sin(self.Latitude0)*sin(self.Longitude0), cos(self.Latitude0)],
                              [ cos(self.Latitude0)*cos(self.Longitude0), cos(self.Latitude0)*sin(self.Longitude0),  sin(self.Latitude0)] ] )
        self.matrixinvers=np.linalg.inv(self.matrix)

    def getECEF2ENU(self, X, Y, Z):
        """
        mapping Global Cartesian to local East/North/Up
        """       

        U=self.matrix[0,0]*(X-self.origin[0])+self.matrix[0,1]*(Y-self.origin[1])+self.matrix[0,2]*(Z-self.origin[2])
        V=self.matrix[1,0]*(X-self.origin[0])+self.matrix[1,1]*(Y-self.origin[1])+self.matrix[1,2]*(Z-self.origin[2])
        W=self.matrix[2,0]*(X-self.origin[0])+self.matrix[2,1]*(Y-self.origin[1])+self.matrix[2,2]*(Z-self.origin[2])
        
        return U,V,W

    def getENU2ECEF(self, U, V, W):
        """
        mapping local East/North/Up to Global Cartesian
        """
        X=self.matrix[0,0]*U+self.matrix[0,1]*V+self.matrix[0,2]*W+self.origin[0]
        Y=self.matrix[1,0]*U+self.matrix[1,1]*V+self.matrix[1,2]*W+self.origin[1]
        Z=self.matrix[2,0]*U+self.matrix[2,1]*V+self.matrix[2,2]*W+self.origin[2]
        
        return X, Y, Z
    
        
    def getPosition2ENU(self, lon, lat,  h):
        """
        mapping global position  and height  to local E/N/U
        """     
        return self.getECEF2ENU(*self.ECEF.getXYZ(lon, lat, h))

    def getENU2Position(self, U, V, W):
        """
        mapping local E/N/U to global position
        """    
        return self.ECEF.getPosition( self.getENU2ECEF(U, V, W) )
 

class StereographicProjection(object):
    """
    modified from https://github.com/nsidc/polar_stereo/blob/main/source/polar_convert.py
    """
    def __init__(self, trueScaleLatitude=70, hemisphere=1, model=WGS84(), range360=False):
        """
         self.trueScaleLatitude (float): true-scale latitude in degrees
         hemisphere (1 or -1): 1 for Northern hemisphere, -1 for Southern
        """
        self.model=model
        self.trueScaleLatitude=abs(trueScaleLatitude)
        if hemisphere > 0:
            self.hemisphere=1
        else:
            self.hemisphere=-1
        self.range360=range360
        
    def toPosition(self, x, y):
        """Convert from Polar Stereographic (x, y) coordinates to
        geodetic longitude and latitude.
        Args:
            x (float): X coordinate(s) in m
            y (float): Y coordinate(s) in m
        Returns: [longitude, latitude] in degs.
        """
        e=self.model.e
        re=self.model.EquatorialRadius
    
        e2 = e * e
        slat = np.deg2rad(self.trueScaleLatitude)
        rho = sqrt(x ** 2 + y ** 2)

        if abs(self.trueScaleLatitude - 90.) < 1e-5:
            t = rho * sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e)) / (2 * re)
        else:
            cm = cos(slat) / sqrt(1 - e2 * (sin(slat) ** 2))
            t = np.tan((np.pi / 4) - (slat / 2)) / ((1 - e * sin(slat)) / (1 + e * sin(slat))) ** (e / 2)
            t = rho * t / (re * cm)
        chi = (np.pi / 2) - 2 * atan(t)
        lat = chi + ((e2 / 2) + (5 * e2 ** 2 / 24) + (e2 ** 3 / 12)) * sin(2 * chi) + \
                    ((7 * e2 ** 2 / 48) + (29 * e2 ** 3 / 240)) * sin(4 * chi) + \
                    (7 * e2 ** 3 / 120) * sin(6 * chi)

        lat = np.rad2deg(1)* self.hemisphere * lat 
        lon = atan2fix(self.hemisphere * x, -self.hemisphere * y)

        lon = np.rad2deg(1) * self.hemisphere * lon

        #lon = lon + np.less(lon, 0) * 360
        if self.range360:
            lon = lon + whereNegative(lon) * 360
        return lon, lat  


    def toXY(self, longitude, latitude):
        """
        Convert from geodetic longitude and latitude to Polar Stereographic (X, Y) coordinates
        Args: 
            longitude (float): longitude or longitude array in degrees
            latitude (float): latitude or latitude array in degrees (positive)
        Returns:
                x, y in m.
        """

        e=self.model.e
        re=self.model.EquatorialRadius
    
        lat = abs(latitude) *  np.deg2rad(1)
        lon = longitude * np.deg2rad(1)
        slat = np.deg2rad(self.trueScaleLatitude)

        e2 = e * e

        # Snyder (1987) p. 161 Eqn 15-9
        t = tan(np.pi / 4 - lat / 2) / ((1 - e * sin(lat)) / (1 + e * sin(lat))) ** (e / 2)

        if abs(90 - self.trueScaleLatitude) < 1e-5:
            # Snyder (1987) p. 161 Eqn 21-33
            rho = 2 * re * t / sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e))
        else:
            # Snyder (1987) p. 161 Eqn 21-34
            tc = np.tan(np.pi / 4 - slat / 2) / ((1 - e * sin(slat)) / (1 + e * sin(slat))) ** (e / 2)
            mc = cos(slat) / sqrt(1 - e2 * (sin(slat) ** 2))
            rho = re * mc * t / tc

        x = rho * self.hemisphere * sin(self.hemisphere * lon)
        y = -rho * self.hemisphere * cos(self.hemisphere * lon)
        return x, y


class RasterData(object):
    GRIDTOL=1e-6
    def __init__(self, FN, X,Y, isCellCentered=False, **data):
        self.FileName=FN

        # find grid:
        
        nx=np.argmax( abs(X[1:]-X[0])>self.GRIDTOL*abs(X).max())+1
        print(nx)
        if nx == 1:
            nx=np.argmax( abs(Y[1:]-Y[0])>self.GRIDTOL*abs(X).max())+1          
            #raise ValueError("Grid identification along X failed. Code is missing here!")
            self.X=X.reshape((-1, nx)).T
            self.Y=Y.reshape((-1, nx)).T
            self.data={ d : data[d].reshape((-1, nx)).T for d in data }
        else:
            self.X=X.reshape((-1, nx))
            self.Y=Y.reshape((-1, nx))
            self.data={ d : data[d].reshape((-1, nx)) for d in data }
        self.interpolator={ d : None for d in data }
        self.dimensions=self.X.shape
        self.isCellCentered=isCellCentered
        if self.isCellCentered:
           self.spacingX=(self.X.max()-self.X.min())/(self.dimensions[0]-1)
           self.spacingY=(self.Y.max()-self.Y.min())/(self.dimensions[1]-1)
           self.rangeX=(self.X.min()-self.spacingX/2, self.X.max()+self.spacingX/2)
           self.rangeY=(self.Y.min()-self.spacingY/2, self.Y.max()+self.spacingY/2)
        else:
           self.spacingX=(self.rangeX[1]-self.rangeX[0])/(self.dimensions[0]-1)
           self.spacingY=(self.rangeY[1]-self.rangeY[0])/(self.dimensions[1]-1)
           self.rangeX=(self.X.min(), self.X.max()+self.spacingX/2)
           self.rangeY=(self.Y.min(), self.Y.max()+self.spacingY/2)

        
    def __len__(self):
        return self.X.size
    def getGridDimensions(self):
        return self.dimensions
    def getRange(self):
        return [self.rangeX, self.rangeY]
    def getGridSpacing(self):
        return self.spacingX, self.spacingY
    
    def getData(self, x,y, dataset):
        if self.interpolator[dataset] == None:
            self.interpolator[dataset] = RegularGridInterpolator((self.X[:,0], self.Y[0,:]), self.data[dataset],  method='linear', bounds_error=False, fill_value=np.mean(self.data[dataset]))

        return self.interpolator[dataset]((x,y))

    def interpolate2DCellsToEscript(self, xy, target, mask, data):
        """
        this function interpolates 
        """
        data=self.data[data]
        X=xy[0]
        Y=xy[1]
        assert self.isCellCentered
        assert X.getFunctionSpace() == target.getFunctionSpace()
        assert Y.getFunctionSpace() == target.getFunctionSpace()
        assert mask.getFunctionSpace() == target.getFunctionSpace()
        
        dx=self.spacingX
        dy=self.spacingY
        nx=self.dimensions[0]
        ny=self.dimensions[1]
        xmin=self.rangeX[0]
        ymin=self.rangeY[0]
        
        for p in range( target.getNumberOfDataPoints() ):

            if mask.getTupleForDataPoint(p)[0] > 0:
                x,y=X.getTupleForDataPoint(p)[0] ,Y.getTupleForDataPoint(p)[0] 
                ix=min(nx-1,max(0,int((x-xmin)/dx)))
                iy=min(ny-1,max(0,int((y-ymin)/dy)))
                #print(x,y,ix,iy)
                target.setValueOfDataPoint(p, data[ix, iy])
                #print(x,y,ix,iy,data[iy, ix])
        return target        
    
class GravityData(RasterData):
    def __init__(self, FN,  isCellCentered=True, delimiter=',', skiprows=1):
        Long, Lat, X,Y, Grav, Height = np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        super().__init__(FN, X,Y, isCellCentered=isCellCentered, Gravity=Grav, Height=Height)
    
    def getGravityData(self):
        return self.data['Gravity']
    def getHeightData(self):
        return self.data['Height']
    def getHeight(self, x, y):
        return self.getData(x,y, 'Height')
    def interpolateGravityToEscript(self, xy, target, mask):
        return self.interpolate2DCellsToEscript(xy, target, mask, 'Gravity')
    
class Surface(RasterData):
    def __init__(self, FN,  isCellCentered=True, delimiter=',', skiprows=1):
        X,Y, Topography, Roh_c, M_surf  = np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        super().__init__(FN, X,Y,  isCellCentered=isCellCentered, Topography=Topography, Density_Crust=Roh_c, Density_Surface=M_surf)
    def getTopographyData(self):
        return self.data['Topography']
    def getTopography(self, x,y):
        return self.getData(x,y, 'Topography')
    
    def getDensityCrustData(self):
        return self.data['Density_Crust']
    def getDensityCrust(self, x, y):
        return self.getData(x,y,'Density_Crust')

    def getDensitySurfaceData(self):
        return self.data['Density_Surface']
    def geDensitySurface(self, x, y):
        return self.getData(x,y,'Density_Surface')
    
class Moho(RasterData):
    def __init__(self, FN,  isCellCentered=True, delimiter=',', skiprows=1):
        X,Y, Moho, Roh_c, Roh_m  = np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        super().__init__(FN, X,Y, isCellCentered=isCellCentered, Moho=Moho, Density_Crust=Roh_c, Density_Mantel=Roh_m)

    def getMohoData(self):
        return self.data['Moho']
    def getMoho(self, x, y):
        return self.getData(x,y,'Moho')
    
    def getDensityCrustData(self):
        return self.data['Density_Crust']
    def getDensityCrust(self, x, y):
        return self.getData(x,y,'Density_Crust')    

    def getDensityMantelData(self):
        return self.data['Density_Mantel']
    def getDensityMantel(self, x, y):
        return self.getData(x,y,'Density_Mantel') 

class LAB(RasterData):
    def __init__(self, FN,  isCellCentered=True, delimiter=',', skiprows=1):
        X,Y, LAB= np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        super().__init__(FN, X,Y, isCellCentered=isCellCentered, LAB=LAB)

    def getLABData(self):
        return self.data['LAB']
    def getLAB(self, x, y):
        return self.getData(x,y,'LAB')

        
class Mantle(RasterData):
    GRIDTOL=1e-6
    def __init__(self, FN, delimiter=',', skiprows=1):
        tab=np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        self.FileName=FN
        X=tab[0]
        Y=tab[1]
        
        nx=np.argmax( abs(X[1:]-X[0])>self.GRIDTOL*abs(X).max())+1
        print(nx)
        if nx == 1:
            #raise ValueError("Grid identification along X failed. Code is missing here!")
            tab = tab[:, tab[0].argsort(kind='mergesort')] # First sort doesn't need to be stable.
            X=tab[0]
            Y=tab[1]
            nx=np.argmax( abs(X[1:]-X[0])>self.GRIDTOL*abs(X).max())+1
            print(nx)
            
            self.X=X.reshape((-1, nx)).T
            self.Y=Y.reshape((-1, nx)).T
        else:
           self.X=X.reshape((-1, nx))
           self.Y=Y.reshape((-1, nx))
        self.dimensions=self.X.shape
        self.h=tab[2:].reshape((-1,)+self.dimensions)[::2]
        self.data=tab[2:].reshape((-1,)+self.dimensions)[1::2]/1000
        self.rangeX=(self.X.min(), self.X.max())
        self.rangeY=(self.Y.min(), self.Y.max())
        self.spacingX=(self.rangeX[1]-self.rangeX[0])/(self.dimensions[0]-1)
        self.spacingY=(self.rangeY[1]-self.rangeY[0])/(self.dimensions[1]-1)    
        
        self.lenh=np.argmax(np.isnan(self.h), axis=0)
        self.lenh=np.where(self.lenh==0, self.h.shape[0], self.lenh )

    def getDensity(self, x, y):
        X=x;
        Y=y;
        #assert X.getFunctionSpace() == target.getFunctionSpace()
        #assert Y.getFunctionSpace() == target.getFunctionSpace()
        #assert mask.getFunctionSpace() == target.getFunctionSpace()
        
        dx=self.spacingX
        dy=self.spacingY
        nx=self.dimensions[0]
        ny=self.dimensions[1]
        xmin=self.rangeX[0]
        ymin=self.rangeY[0]
        
        ix=min(nx-1,max(0,int((x-xmin)/dx)))
        iy=min(ny-1,max(0,int((y-ymin)/dy)))
        
        ds=self.data[:self.lenh[ix,iy],ix,iy][::-1]
        return ds
        
    def getDepth(self, x, y):
        X=x;
        Y=y;
        #assert X.getFunctionSpace() == target.getFunctionSpace()
        #assert Y.getFunctionSpace() == target.getFunctionSpace()
        #assert mask.getFunctionSpace() == target.getFunctionSpace()
        
        dx=self.spacingX
        dy=self.spacingY
        nx=self.dimensions[0]
        ny=self.dimensions[1]
        xmin=self.rangeX[0]
        ymin=self.rangeY[0]
        
        ix=min(nx-1,max(0,int((x-xmin)/dx)))
        iy=min(ny-1,max(0,int((y-ymin)/dy)))
        
        hs=self.h[:self.lenh[ix,iy],ix,iy][::-1]
        return hs
    
class DensityModel(object):
    GRIDTOL=1e-6
    def __init__(self, FN, delimiter=',', skiprows=1):
        tab=np.loadtxt(FN, delimiter=delimiter, skiprows=skiprows, unpack=True)
        self.FileName=FN
        X=tab[0]
        Y=tab[1]
        
        nx=np.argmax( abs(X[1:]-X[0])>self.GRIDTOL*abs(X).max())+1
        if nx == 1:
            #raise ValueError("Grid identification along X failed. Code is missing here!")
            self.X=X.reshape((-1, nx)).T
            self.Y=Y.reshape((-1, nx)).T
        else:
           self.X=X.reshape((-1, nx))
           self.Y=Y.reshape((-1, nx))
        self.dimensions=self.X.shape
        self.h=tab[2:].reshape((-1,)+self.dimensions)[::2]
        self.data=tab[2:].reshape((-1,)+self.dimensions)[1::2]
        self.rangeX=(self.X.min(), self.X.max())
        self.rangeY=(self.Y.min(), self.Y.max())
        self.spacingX=(self.rangeX[1]-self.rangeX[0])/(self.dimensions[0]-1)
        self.spacingY=(self.rangeY[1]-self.rangeY[0])/(self.dimensions[1]-1)    
        
        self.lenh=np.argmax(np.isnan(self.h), axis=0)
        self.lenh=np.where(self.lenh==0, self.h.shape[0], self.lenh )
    
    def interpolateToEscript(self, xyh, target, topography, moho, mantle):
        X=xyh[0]
        Y=xyh[1]
        H=xyh[2]

        assert X.getFunctionSpace() == target.getFunctionSpace()
        assert Y.getFunctionSpace() == target.getFunctionSpace()
        #assert mask.getFunctionSpace() == target.getFunctionSpace()
        
        dx=self.spacingX
        dy=self.spacingY
        nx=self.dimensions[0]
        ny=self.dimensions[1]
        xmin=self.rangeX[0]
        ymin=self.rangeY[0]
        

         
        for p in range( target.getNumberOfDataPoints() ):
            x,y, h=X.getTupleForDataPoint(p)[0] ,Y.getTupleForDataPoint(p)[0], H.getTupleForDataPoint(p)[0]
            if h > topography.getTopography(x,y):
                dd=0.
            elif h < moho.getMoho(x,y):
                #dd=moho.getDensityMantel(x,y)
                dd=np.interp(h, mantle.getDepth(x,y), mantle.getDensity(x,y), left=3.3, right=moho.getDensityMantel(x,y), period=None)

            else:
                # CC and SSB Density
                rho_c_top=topography.getDensityCrust(x,y)
                # CC Density
                rho_c_bottom=moho.getDensityCrust(x,y)
                ix=min(nx-1,max(0,int((x-xmin)/dx)))
                iy=min(ny-1,max(0,int((y-ymin)/dy)))
                hs=self.h[:self.lenh[ix,iy],ix,iy][::-1]
                ds=self.data[:self.lenh[ix,iy],ix,iy][::-1]
                #dd=np.interp(h, self.h[:self.lenh[iy,ix],iy,ix], self.data[:self.lenh[iy,ix],iy,ix], right=0., left=self.data[self.lenh[iy,ix]-1,iy,ix], period=None)
                dd=np.interp(h, hs, ds, left=rho_c_bottom, right=rho_c_top, period=None)

                #print(h, dd, moho.getMoho(x,y), hs, ds, rho_c_top, rho_c_bottom)
                #print(dd, h, self.lenh[iy,ix], hs, ds)
            target.setValueOfDataPoint(p, dd)
        return target   


    def maskToEscript(self, xyh, target, lab):
        X=xyh[0]
        Y=xyh[1]
        H=xyh[2]

        assert X.getFunctionSpace() == target.getFunctionSpace()
        assert Y.getFunctionSpace() == target.getFunctionSpace()
        #assert mask.getFunctionSpace() == target.getFunctionSpace()
        
        dx=self.spacingX
        dy=self.spacingY
        nx=self.dimensions[0]
        ny=self.dimensions[1]
        xmin=self.rangeX[0]
        ymin=self.rangeY[0]
        

         
        for p in range( target.getNumberOfDataPoints() ):
            x,y, h=X.getTupleForDataPoint(p)[0] ,Y.getTupleForDataPoint(p)[0], H.getTupleForDataPoint(p)[0]
            if h < lab.getLAB(x,y):
                dd=0.
                target.setValueOfDataPoint(p, dd)
        return target   
