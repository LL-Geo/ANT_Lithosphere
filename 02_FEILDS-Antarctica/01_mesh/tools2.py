#!/usr/bin/python3
__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
sys.path.insert(0, os.getcwd())

from esys.escript import *
from esys.finley import ReadMesh
import esys.escript.unitsSI as U
import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import saveSilo
from esys.weipa import saveVTK
from esys.escript.pdetools import MaskFromBoundaryTag
from esys.escript.pdetools import Locator
from tools import *

class curvedGravityModel(object):
    """
    This class is a simple wrapper for a 3D gravity PDE model
    with earth curvature.  
    It solves PDE
        - div (grad u) = -4 pi G rho
    where 
       G  is the gravitational constant and 
       rho is density
       u  is computed anomaly potential
 
    Possible boundary conditions included in the class are Dirichlet on 
       - top (default) or  
       - top and base.
       
    It has a function to set ground property 
       - setDensity
    get ground property
       - getDensity
    and functions to output solutions
       - getgravityPotential : u
       - getDownGravity : grad(u) . gdir
       - getGravityVector : grad (u) .
    """
    
    def __init__(self, domain, fixBase = False):
        """
        Initialise the class with domain and boundary conditions.
        Setup PDE and density.
        : param domain: the domain 
        : type domain: `Domain`
        : param fixBase: if true the gravitational potential at the bottom is set to zero.
        : type fixBase: `bool`
        : param fixVert: if true the magnetic field on all vertical sudes is set to zero.
        : type fixBase: `bool`
        : if fixBase is True then gravity field is set to zero at base and top surfaces.
        """
        self.domain = domain
        self.fixBase=fixBase
        assert self.domain.getDim() == 3
        self.pde=self.__createPDE()
        self.setDensity()
        

    def __createPDE(self):
        """
        Create the PDE and set boundary conditions.
        """
        pde = LinearSinglePDE(self.domain, isComplex=False)
        domdim = self.domain.getDim()       
        optionsG = pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        pde.setSymmetryOn()
        pde.setValue(A = kronecker(domdim))
        if self.fixBase:
            bounds=MaskFromBoundaryTag(self.domain, "SurfaceTop","SurfaceBottom")
        else:
            bounds=MaskFromBoundaryTag(self.domain, "SurfaceTop")
        pde.setValue(q=bounds)
        if hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
        return pde
             
    def setDensity(self, rho=0):
        """
        set density
        : param rho: density
        : type rho: `Data` or `float` 
        """
        self.rho = rho
        self.reset = True

    def getDensity(self):
        """
        returns density
        : returns: rho
        """
        return self.rho

    def setSurfDensity(self, rho_surf=0):
        """
        set density
        : param rho: density
        : type rho: `Data` or `float` 
        """
        self.rho_surf = rho_surf
        self.reset = True

    def getSurfDensity(self):
        """
        returns density
        : returns: rho
        """
        return self.rho_surf
        
    def getGravityPotential(self):
        """
        get the potential of the density anomaly
        """
        if self.reset:
            self.pde.setValue(y = -4.0*np.pi*U.Gravitational_Constant*self.rho_surf) 
            self.pde.setValue(Y = -4.0*np.pi*U.Gravitational_Constant*self.rho)
        return self.pde.getSolution()

    def getDownGravity(self, down=-kronecker(3)[2]):
        """
        get Bouger gravity in down direction 
        """
        return inner(self.getGravityVector(), down)

    def getDownGravityDataArea(self, dataArea, down=-kronecker(3)[2]):
        """
        get Bouger gravity in down direction in dataArea
        """
        return dataArea * self.getDownGravity(down)
 
    def getGravityVector(self):
        """
        get the Bouger gravity vector
        """
        return grad(self.getGravityPotential(), ReducedFunction(self.pde.getDomain()))

################################################################
################################################################

class ACGravity(object):
    def __init__(self, domain, dataArea, g_dir, g_data, rho_0, ground, mu, m0, fixBase,
         atol=1.0, rtol=1.0, iter_max=100, pde_tol=1e-8, output_name='solutions', verboseLevel="low"):

        self.domain = domain   
     
        self.g_dir = g_dir                   # gravity direction vector everywhere
        self.dataArea = dataArea             # 1 in data area
        self.w = self.dataArea*self.g_dir    # gravity direction vector only in data area 
        self.g_data =g_data                  # data for inversion

        self.rho_0 = rho_0                   # reference density for inversion     
        self.beta = 4.0*np.pi*U.Gravitational_Constant*ground

        self.m0   = m0       
        self.atol = atol
        self.rtol = rtol
        self.iter_max = iter_max
        self.pdetol = pde_tol
        self.verboseLevel = verboseLevel
        self.output_name = output_name

        # set up gravity pde
        self.pdeu = self.setupPDE()
        if fixBase:
            bounds=MaskFromBoundaryTag(domain, "SurfaceTop","SurfaceBottom")
        else:
            bounds=MaskFromBoundaryTag(domain, "SurfaceTop")
        self.pdeu.setValue(q=bounds)

        self.gg   = integrate(self.g_data*self.g_data)                                 
        self.mu   = mu*self.gg

    def setupPDE(self):
        pde=LinearSinglePDE(self.domain, isComplex=False)
        pde.setValue(A=kronecker(3))  
        pde.setSymmetryOn()
        options=pde.getSolverOptions()
        options.setTolerance(self.pdetol)
        options.setSolverMethod(SolverOptions.PCG)
        options.setPackage(SolverOptions.TRILINOS)
        options.setPreconditioner(SolverOptions.AMG)
        options.setTrilinosParameter("max levels", 10)  
        options.setTrilinosParameter("problem: symmetric", True)
        options.setTrilinosParameter("reuse: type", "full")
        return pde

    def RHS(self): 
        self.pdeu.setValue(X = self.g_data*self.g_dir, Y = Scalar(0., Solution(self.domain)))
        u = self.pdeu.getSolution()
        return ArithmeticTuple (-(self.beta/self.mu)*u, Vector(0., Solution(self.domain)))         

    def Aprod(self,m):
        self.pdeu.setValue(Y = -self.beta*m, X = Vector(0., Solution(self.domain)))   
        u = self.pdeu.getSolution()                            
        self.pdeu.setValue(Y = Scalar(0., Solution(self.domain)), X = self.w*inner(grad(u),self.w))   
        u = self.pdeu.getSolution()
        return ArithmeticTuple(-(self.beta/self.mu)*u, grad(m))     
        
    def Msolve(self,r):
        self.pdeu.setValue(X=r[1], Y=r[0])
        return self.pdeu.getSolution()       

    def bilinearform(self, m, r):
        return integrate(inner(grad(m), r[1])+ m*r[0])

    def myPCG(self, x,r,itermax,rtol):
       mfs = []
       smooths = []
       piter=0
       rhat=self.Msolve(r)
       d = rhat
       rhat_dot_r = self.bilinearform(rhat, r)
       if rhat_dot_r<0: print("negative norm.")
       norm_r0=np.sqrt(rhat_dot_r)
       atol2=self.rtol*norm_r0
       if atol2<=0:
          print("Non-positive tolarance.")
       print(("PCG: initial residual norm = %e (absolute tolerance = %e)"%(norm_r0, atol2)))
       mfold=0.5*integrate(self.g_data**2)/self.gg
       smooths.append(0.)
       mfs.append(mfold)
       donep05=False
       donep01=False
       donep008=False
       donep005=False
       while not np.sqrt(rhat_dot_r) <= atol2:
           piter+=1
           if piter  >= self.iter_max: 
               print("maximum number of %s steps reached."%iter_max)
               break 
           q=self.Aprod(d)
           alpha = rhat_dot_r / self.bilinearform(d, q)
           x += alpha * d
           r += q * (-alpha)      
           rhat = self.Msolve(r)
           rhat_dot_r_new = self.bilinearform(rhat, r)
           beta = rhat_dot_r_new / rhat_dot_r
           rhat += beta * d
           d = rhat

           rhat_dot_r = rhat_dot_r_new
           if rhat_dot_r<0: print("negative norm.")
           self.pdeu.setValue(Y = -self.beta*x, X = Vector(0., Solution(self.domain)))
           u = self.pdeu.getSolution()
           g_vector = grad(u)
           g_comp = inner(g_vector, self.g_dir)       
           g_comp_dataArea = self.dataArea*g_comp
           mf=0.5*integrate((g_comp_dataArea - self.g_data)**2)/self.gg
           if mf<0.05 and not donep05:
               print("0.05",piter)
               saveVTK(self.output_name+"p05", gravity=g_comp, m=x*self.rho_0)
               donep05=True
           if mf<0.01 and not donep01:
               print("0.01",piter)
               saveVTK(self.output_name+"p01", gravity=g_comp, m=x*self.rho_0)
               donep01=True
           if mf<0.008 and not donep008:
               print("0.008",piter)
               saveVTK(self.output_name+"p008", gravity=g_comp, m=x*self.rho_0)
               donep008=True
           if mf<0.005 and not donep005:
               print("0.005",piter)
               saveVTK(self.output_name+"p005", gravity=g_comp, m=x*self.rho_0)
               donep005=True
           mfs.append(mf)
           smooth=integrate(inner(grad(x),grad(x)))
           smooths.append(smooth)
           print('mf ', mf,' smooth ',smooth)
           print(("PCG: iteration step %s: residual norm = %e"%(piter, np.sqrt(rhat_dot_r))))
       print(("PCG: tolerance reached after %s steps."%piter))
       smooths=np.array(smooths)
       mfs=np.array(mfs)
       np.savetxt(self.output_name+"smooths.csv",smooths, delimiter=",")   #
       np.savetxt(self.output_name+"mfs.csv",mfs,delimiter=",")
       return x,

    def solve(self):
        r = self.RHS()
        if self.verboseLevel=="low":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=False)
        elif self.verboseLevel=="medium":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=True)
        elif self.verboseLevel == "high":
            m = self.myPCG(self.m0, r,self.iter_max,self.rtol)

        rho_final = self.rho_0 * m[0]
        print("density range",inf(rho_final),sup(rho_final))
         
        return m[0], rho_final











