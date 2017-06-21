# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 11:03:48 2017

@author: Karin
"""

#-----------------------------IMPORTERING--------------------------------------
import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu
import cfext as tg
import visvis as vv

from PyQt5 import QtGui, QtCore
from PyQt5.QtGui import QPixmap



#-----------------------------INPUTDATA----------------------------------------
class InputData(object):
    """Klass för att definiera modellens indata."""
    
    def __init__(self):
        """Värden för beräkning om inget annat anges"""
    
        # --- Inputdata tab 0       
        self.b = 2  
        self.B = 4*self.b
        self.q = 30e3
        self.elSize = 0.5
        self.elSizenr = 5

        self.perm = True
        
        self.hn = 1                                       #Antal lager
        self.hvec = 15*np.ones([self.hn,1])          #Vektor med lagertjocklekar                    
        self.kmat = 1e-9*np.ones([self.hn,2])
        self.Mvec = 900e3*np.ones([self.hn,1])

        self.rhow = 10e3
        
        
        # --- Time approximation tab 1
        self.T = 2          #Max time  
        self.dt = 0.025       # Time step in approx.
        self.Tnr = 8       #Number of screenshots
        self.Tv = 0         #Number to plot
        
        # --- 1D plot
        self.Xnr = 10       #Nr of stops slider x    
        self.Xv = 0         #Nr to plot
        self.Ynr = 10       #Nr of stops slider y
        self.Yv = 0          #Nr to plot

        self.warnings = {}
        self.comboindex = 0
        

        
    def geometry(self):
        """Skapa en geometri-instans baserad för definierade parametrar"""

        # --- Geometri-instans för geometribeskrivning
        g = cfg.Geometry()
        hn = self.hn        
        hvec = self.hvec
        b = self.b
        B = self.B
        
        # --- Punkter

        g.point([b, 0])
        g.point([B, 0])
        g.point([0, 0])
        
        h = 0
        for i in range(hn):
            h=h - float(hvec[i])
            g.point([B, h])
            g.point([0, h])

        # --- Linjer
        g.spline([1, 0], marker = 2) #markyta
        g.spline([0, 2], marker = 3) #under last
       
        a = 0
        b = 500
        c = 600
        d = 700
        for i in range (hn):
            g.spline([a + 2, a + 4], marker = b) #CL 
            g.spline([a + 4, a + 3], marker = c) #underkant  
            g.spline([a + 3, a + 1], marker = d)
            a = a + 2
            b = b + 1
            c = c + 1
            d = d + 1 
          
        # --- Yta
        g.surface([0, 1, 2, 3, 4], marker = 100)
        
        svec=np.array([3, 5, 6, 7])
        a = 101
        for i in range(hn - 1):
            g.surface(svec, marker = a)
            svec = svec + 3
            a = a + 1
            
        # --- Returnera skapade geometrin
        return g



#-------------------------OUTPUTDATA-------------------------------------------        
class OutputData(object):
    """Klass för att lagra resultaten från beräkningen."""
    
    def __init__(self):
        """ Definiera tomma variabler innan de tilldelats värden."""
        
        self.geometry = None
        
        #Mesh
        self.coords = None  
        self.edof = None
        self.bdofs = None
        self.dofs = None    
        self.elType = None        
        self.dofsPerNode = None   
        self.elementmarkes = None
        self.ex = None
        self.ey = None



        #FEM
        self.K = None
        self.C = None
        self.a0 = None
        self.bc = None
        self.bcVal = None
        
        
        #Store values to plot when slider changed
        self.Tvec = None
        self.avec = None
        self.xudict = None        
        self.xlvec = None
        self.yudict = None
        self.ylvec = None

#-------------------------SOLVER----------------------------------------------- 
class Solver(object):
    """Klass för att sköta beräkningarna."""
    
    def __init__(self, inputData, outputData):
        """ Hämta inputdata och outputdata."""
        
        self.inputData = inputData
        self.outputData = outputData

    def executeMesh(self):
        """ Utför FE-beräkningen."""
        # --- Överför modellvariabler till lokala referenser
        geometry = self.inputData.geometry()
        elSize = self.inputData.elSize
                     
        
        # --- Skapa mesh
        dofsPerNode = 1
        elType = 2      # <-- Trenodselement 
        meshGen = cfm.GmshMeshGenerator(geometry)
        meshGen.elSizeFactor = elSize
        meshGen.elType = elType
        meshGen.dofsPerNode = dofsPerNode
        coords, edof, dofs, bdofs, elementmarkers = meshGen.create()
        
        # --- Bestäm ex, ey               
        ex, ey = cfc.coordxtr(edof, coords, dofs)        
        
        
        #Returnera framtagen data
        self.outputData.coords = coords  
        self.outputData.edof = edof
        self.outputData.dofs = dofs
        self.outputData.bdofs = bdofs    
        self.outputData.elType = elType
        self.outputData.dofsPerNode = dofsPerNode
        self.outputData.elementmarkers = elementmarkers
        self.outputData.ex = ex
        self.outputData.ey = ey
        


        
    def executeFem(self):
        
        q = self.inputData.q
        B = self.inputData.B
        b = self.inputData.b
        hvec = self.inputData.hvec
        Mvec = self.inputData.Mvec
        kmat = self.inputData.kmat
        rhow = self.inputData.rhow
        hn = self.inputData.hn
        
        coords = self.outputData.coords
        edof = self.outputData.edof
        bdofs = self.outputData.bdofs
        elementmarkers = self.outputData.elementmarkers 
        ex = self.outputData.ex
        ey = self.outputData.ey
        
        ndof=edof.max()
        
        # --- Definera randvillkor
        bc = np.array([],int)
        bcVal = np.array([],int)
        
        bc,bcVal = cfu.applybc(bdofs, bc, bcVal, 3, 0) 
        bc,bcVal = cfu.applybc(bdofs, bc, bcVal, 2, 0)
        
    
        
        if self.inputData.perm == True:
            bc,bcVal = cfu.applybc(bdofs, bc, bcVal, (599+hn), 0)    
        
        # --- Definiera materialmatris       
        ep=[(1)]
        D=np.zeros([2,2]) 
        
        # --- Skapa tom K-matris och f-vektor
        
        K = np.zeros([ndof,ndof])   
        C = np.zeros([ndof,ndof])


        
        for eltopo, elx, ely, elMarker in zip(edof, ex, ey, elementmarkers):
            Dn = elMarker - 100
            D[0,0] = kmat[Dn,0]*60*60*24*365.25
            D[1,1] = kmat[Dn,1]*60*60*24*365.25

            M = Mvec[Dn]
            
            Ke = cfc.flw2te(elx, ely, ep, D)
            Ce = tg.flwtec(elx,ely,rhow,M)
            cfc.assem(eltopo, K, Ke)
            cfc.assem(eltopo,C,Ce)
    
        a0 = tg.strdist(coords, ndof, q, B, b, hvec)
        a0 = np.asmatrix(a0)
        
        self.outputData.K = K
        self.outputData.C = C
        self.outputData.a0 = a0
        self.outputData.bc = bc
        self.outputData.bcVal = bcVal
    
        

        
    def executeTime(self):

        dt = self.inputData.dt
        T = self.inputData.T
        Tnr = self.inputData.Tnr  
        

        
        K = self.outputData.K
        C = self.outputData.C
        a0 = self.outputData.a0
        bc = self.outputData.bc
        bcVal = self.outputData.bcVal 
        

       #Beräkna avec =vektor med samtliga noders u-värde vid tiderna i Tvec

        avec, Tvec = tg.step(K, C, a0, bc, bcVal, dt, T, Tnr)
        
               
        self.outputData.avec = avec
        self.outputData.Tvec = Tvec 
        
    def executeExtr(self):
        B = self.inputData.B
        hvec = self.inputData.hvec
        Xnr = self.inputData.Xnr  
        Ynr = self.inputData.Ynr
        


        
        ex = self.outputData.ex
        ey = self.outputData.ey
        edof = self.outputData.edof
        avec = self.outputData.avec
        ndof = edof.max()       
        xlvec = np.linspace(0, B, Xnr)
        ylvec = np.linspace(0, np.sum(hvec), Ynr)
        surfx = np.array([B])
        surfy = np.cumsum(hvec)
        
        
        xudict = {}
        dictnr = 0
        for x in xlvec:
            yvec, uvec = tg.xlinextr(x, surfx, ex, ey, ndof, edof, avec)
            idx = np.argsort(yvec)
            yvec = yvec[idx]
            uvec = uvec[idx] 
            xudict[dictnr] = yvec, uvec
            dictnr = dictnr + 1
            
        yudict = {}
        dictnr = 0
        for y in ylvec:     
            xvec, uvec = tg.ylinextr(y, surfy, ex, ey, ndof, edof, avec)
            idx = np.argsort(xvec)
            xvec = xvec[idx]
            uvec = uvec[idx] 
            yudict[dictnr] = xvec, uvec
            dictnr = dictnr + 1         
            

        self.outputData.xudict = xudict
        self.outputData.xlvec = xlvec
        self.outputData.yudict = yudict
        self.outputData.ylvec = ylvec
        
        

    
    
#-------------------------VISUALISATION-------------------------------------------      
class Visualisation(object):
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData

    def show(self):

        coords = self.outputData.coords
        edof = self.outputData.edof
        dofsPerNode = self.outputData.dofsPerNode
        elType = self.outputData.elType
        avec = self.outputData.avec   
        plotnr = self.inputData.Tv
        
        
        #Plot excess pore water pressure 2D
        cfv.figure(self.figPress.nr) 
        cfv.clf()
        cfv.drawNodalValues(avec[:,plotnr], coords, edof, dofsPerNode, elType, clim=(0,self.inputData.q), doDrawMesh = False, title = "Excess porewater pressure")



        xudict = self.outputData.xudict
        yudict = self.outputData.yudict  
        comboindex = self.inputData.comboindex
  
        PP = np.shape(avec)[1]
        
        
        #If vertical 1 D plot
        if comboindex == 0:
            Xv = self.inputData.Xv
            
            coordvec = xudict[Xv][0]
            uvec = xudict[Xv][1]
            
            cfv.figure(self.figuvec.nr) 
            cfv.clf()
            for i in range(PP):
                vv.plot(uvec[:, i], coordvec, ls = ":", lc = "c") 
            vv.plot(uvec[:,plotnr], coordvec, ls = "-", lc = "b")
            vv.title("u_e - Excess porewater pressure 1D")
            
            
            
            cfv.figure(self.figsigvec.nr) 
            cfv.clf()
            for i in range(PP):
                vv.plot(uvec[:,0]-uvec[:, i], coordvec, ls = ":", lc = "c") 
            vv.plot(uvec[:,0]-uvec[:,plotnr], coordvec, ls = "-", lc = "b")
            vv.title("\Delta\sigma' - Effective stress 1D")
            
            
        else:
            Yv = self.inputData.Yv  
            coordvec = yudict[Yv][0]
            uvec = yudict[Yv][1]
            
            cfv.figure(self.figuvec.nr) 
            cfv.clf()
            for i in range(PP):
               vv.plot(coordvec, uvec[:, i], ls=":", lc = "c") 
            vv.plot(coordvec, uvec[:, plotnr], ls="-", lc = "b")
            vv.title("u_e - Excess porewater pressure 1D")
            
            cfv.figure(self.figsigvec.nr) 
            cfv.clf()
            for i in range(PP):
                vv.plot(coordvec, uvec[:, 0]-uvec[:, i], ls = ":", lc = "c") 
            vv.plot(coordvec, uvec[:, 0]-uvec[:, plotnr], ls = "-", lc = "b")
            vv.title("\Delta\sigma' - Effective stress 1D")
   


    def wait(self):
        """Denna metod ser till att fönstren hålls uppdaterade och kommer att returnera
        När sista fönstret stängs"""

        cfv.showAndWait()
        



      