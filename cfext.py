# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 11:39:05 2017

@author: Karin
"""


import numpy as np
import calfem.core as cfc
        
"""Ny funktion för att beräkna spänningsfördelning under en utbredd linjelast"""
def strdist(coords, ndof, q, B, b, hvec):
    """
    Function to compute stress distrubution in given points (x,y) under a strip load.
    
    Parameters:
    
        coords =       nodal coordinates (x,y)
        ndof =          number of d.o.f.
        q =              magnitude load [N/m2]
        B =              length of the soil to study, measured from center line [m]
        b =              width of the strip load
        hvec =          vector containing the tickness of the different soil layers
        Returns:
    
        sigmaz              Vectrical stress component in each node, dim(sigmaz)= ndof x 1
        
    """

    sigmaz = np.zeros([ndof,1])
    for i in range(ndof):
        x = coords[i, 0]
        y = coords[i, 1]
        A1 = x + b
        B1 = x - b
        C1 = - y

        if C1 != 0:
            beta = np.arctan(B1/C1)
            alpha = np.arctan(A1/C1)-beta
            sigmaz[i] = q/np.pi*(alpha+np.sin(alpha)*np.cos(alpha+2*beta))
        else:
            if x < b:
                beta = 2*np.pi
                alpha = 0
                sigmaz[i]= q
            else:
                beta = 2*np.pi
                alpha = 0
                sigmaz[i]= 0
    return sigmaz

"""Ny funktion för att bestämma kapacitetsmatrisen C"""
def flwtec(ex, ey, rho, M):
    
    """
    Compute the capacity matrix C for a 2D-element with 3 nodes.
    
    Parameters:
    
        rho =           Desity of pore fluid     
        M =             Compression modulus
        ex =            Vector x-coordinates
        ey =            Vector y-coordinates
        
        Returns:
    
        Ce              Capacity matrix, dim(Ce)= 3 x 3
        
    """

    exm = np.asmatrix(ex)
    eym = np.asmatrix(ey)
    ci = np.hstack((np.ones([3,1]), np.transpose(exm), np.transpose(eym)))
    A = 0.5*np.linalg.det(ci)
    ce = np.matrix([[1/6, 1/12, 1/12],
                   [1/12, 1/6, 1/12],
                   [1/12, 1/12, 1/6]])
    Ce=rho/float(M)*A*ce
    return Ce
    
    
"""Ny funktion för att göra en tidsberoende analys i python, jmf med step i CALFEM för matlab """    

def step(K, C, a, bc, bcval, dt, T, Tnr):
    """
    Compute the capacity matrix C for a 2D-element with 3 nodes.
    
    Parameters:
    
        K =             Stiffness matrix   
        C =             Capacity matrix
        a =             Stress distrubution at time t=0
        bc =            Nodes with prescribed values
        bcval =        Prescribed nodal values
        dt =            Time increment size
        Tvec =        Vector containing times for which the stress distrubution  should be displayed
        
        Returns:
    
        avec              Excess pore water pressure at given times (Tvec) in all nodes, dim(avec)= ndof x Tvec
        
    """

    n = int(T/dt)
    ndt = np.floor(T/dt/Tnr)
    Tvec = np.arange(ndt, n, ndt)
    Tvec = np.append(Tvec, n)

    theta = 1
    avec = a
    
    Kt = K
    Ct = C

    

    for i in range (n):

        Kt = Ct + dt*Kt
        f=(Ct+(1-theta)*Kt)*a
        a, r = cfc.solveq(Kt, f, bc)
        if (i+1) in Tvec:
            avec=np.hstack((avec, a))
    return avec, Tvec



def xlinextr(xl, surfx, ex,ey, ndof, edof,avec):
    
    
    """
    xl - x-koordinat där ett snitt ska göras
    ex - x-koordinater för element [nel,3]
    ey - ey-koordinater för element [nel,3]
    edof - [ndof, 3]
    avec - [ndof, tstep]
    FI
    """    

    xls = xl
    
    if xls == 0:
        xls = 0.01
    elif xls in surfx:
        xls = xls - 0.01

    #Antal element att söka över
    nel = int(np.shape(edof)[0])
    elnr = []
    
    #Plocka ut element som skär x-linje
    for i in range(nel):
        if np.min(ex[i,:]) < xls  < np.max(ex[i,:]):
            #Spara elementnummer
            elnr = np.append(elnr, i)
    
    #Antal element som skär linjen samt antal tidssteg i avec
    nelc = int(np.shape(elnr)[0])
    nstep = int(np.shape(avec)[1])
    
 
    exdof = np.zeros([nelc,3])
    eydof = np.zeros([nelc,3])
    eu1dof = np.zeros([nelc,nstep])
    eu2dof = np.zeros([nelc,nstep])
    eu3dof = np.zeros([nelc,nstep])
    
    a = 0
    for i in elnr:
        #Spara x och y koordinater för elementet
         exdof[a,:] = ex[i,:]
         eydof[a,:] = ey[i,:]

         #Plocka ut frihetsgrader och u-väden för dessa vid samtliga tidpunkter
         dofs = edof[i,:]-1

         eu1dof[a,:] = avec[dofs[0],:]
         eu2dof[a,:] = avec[dofs[1],:]
         eu3dof[a,:] = avec[dofs[2],:]
         a = a +1

    
    uvec = np.zeros([1,nstep])
    yvec = []
    for elex, eley, eleu1, eleu2, eleu3 in zip(exdof, eydof, eu1dof, eu2dof, eu3dof):
        eu = np.vstack((eleu1,eleu2,eleu3))

        exv = []
        eyv = []
        euv = np.zeros([1,nstep])

        
        exh = []
        eyh = []
        euh = np.zeros([1,nstep])

    
        for i in range(3):
            
            #Dela upp i två vektorer beroende på om koordinaterna är till vänster eller höger
            if elex[i] < xls:
                exv = np.append(exv, elex[i])
                eyv = np.append(eyv, eley[i])
                euv = np.vstack((euv, eu[i,:]))
    
            else:
                exh = np.append(exh, elex[i])
                eyh = np.append(eyh, eley[i])
                euh = np.vstack((euh, eu[i,:]))
    
            
        euv = euv[1::,:]
        euh = euh[1::,:]

            
        #Om en koordinat till vänster
        if np.shape(exv)[0] == 1:

            # xv-koord konstanta för de två interpoleringarna
            xv = exv[0]
            yv = eyv[0]
            uv = euv[0,:]
            
            
            #interpolera två gånger, varierande xh- värden (två punker)
            for i in range(2):
                xh = exh[i]
                yh = eyh[i]
                uh = euh[i,:]
                
                y = yv + (yh - yv) * (xl - xv) / (xh - xv)
                u = uv + (uh - uv) * (xl - xv) / (xh - xv)
                #spara värden
                yvec = np.append(yvec, y)
                uvec = np.vstack((uvec, u))
        #Om två koordinater till vänster
        elif np.shape(exv)[0] == 2:
            # xh-koord konstanta för de två interpoleringarna
            xh = exh[0]
            yh = eyh[0]
            uh = euh[0,:]

            
            #interpolera två gånger, varierande xv- värden (två punker)
            for i in range(2):
                xv = exv[i]
                yv = eyv[i]
                uv = euv[i,:]

                y = yv + (yh - yv) * (xl - xv) / (xh - xv)
                u = uv + (uh - uv) * (xl - xv) / (xh - xv)
                
                #Spara värden
                yvec = np.append(yvec, y)
                uvec = np.vstack((uvec, u))

    uvec = uvec[1::,:]
    return yvec, uvec
        




def ylinextr(yl, surfy, ex, ey, ndof, edof,avec):
    """
    yl - y-koordinat för slider
    ex - x-koordinater för element [nel,3]
    ey - y-koordinater för element [nel,3]
    edof - [ndof, 3]
    avec - [ndof, tstep]
    FORTSÄTT
    
    """
    
    yl = - yl #Bytt tecken vid byte x ->y
    yls = yl

    if yls ==0:
        yls = - 0.001
        
    #Antal element att söka över
    elif -yls in surfy: 
        yls = yls + 0.001
    
    nel = int(np.shape(edof)[0])
    elnr = []
    
    #Plocka ut element som skär x-linje
    for i in range(nel):
        if np.min(ey[i,:]) < yls  < np.max(ey[i,:]):
            #Spara elementnummer
            elnr = np.append(elnr, i)
    
    #Antal element som skär linjen samt antal tidssteg i avec
    nelc = int(np.shape(elnr)[0])
    nstep = int(np.shape(avec)[1])
    
 
    exdof = np.zeros([nelc,3])
    eydof = np.zeros([nelc,3])
    eu1dof = np.zeros([nelc,nstep])
    eu2dof = np.zeros([nelc,nstep])
    eu3dof = np.zeros([nelc,nstep])
    
    a = 0
    for i in elnr:
        #Spara x och y koordinater för elementet
         exdof[a,:] = ex[i,:]
         eydof[a,:] = ey[i,:]

         #Plocka ut frihetsgrader och u-väden för dessa vid samtliga tidpunkter
         dofs = edof[i,:]-1

         eu1dof[a,:] = avec[dofs[0],:]
         eu2dof[a,:] = avec[dofs[1],:]
         eu3dof[a,:] = avec[dofs[2],:]
         a = a +1

    
    uvec = np.zeros([1,nstep])
    xvec = []
    for elex, eley, eleu1, eleu2, eleu3 in zip(exdof, eydof, eu1dof, eu2dof, eu3dof):
        eu = np.vstack((eleu1,eleu2,eleu3))

        exv = []
        eyv = []
        euv = np.zeros([1,nstep])

        
        exh = []
        eyh = []
        euh = np.zeros([1,nstep])

    
        for i in range(3):
            
            #Dela upp i två vektorer beroende på om koordinaterna är till vänster eller höger
            if eley[i] > yls:  #Bytt olikhet vid x ->y
                exv = np.append(exv, elex[i])
                eyv = np.append(eyv, eley[i])
                euv = np.vstack((euv, eu[i,:]))
    
            else:
                exh = np.append(exh, elex[i])
                eyh = np.append(eyh, eley[i])
                euh = np.vstack((euh, eu[i,:]))
    
            
        euv = euv[1::,:]
        euh = euh[1::,:]

            
        #Om en koordinat till vänster
        if np.shape(exv)[0] == 1:

            # xv-koord konstanta för de två interpoleringarna
            xv = exv[0]
            yv = eyv[0]
            uv = euv[0,:]
            
            
            #interpolera två gånger, varierande xh- värden (två punker)
            for i in range(2):
                xh = exh[i]
                yh = eyh[i]
                uh = euh[i,:]
                
                x = xv + (xh - xv) * (yl - yv) / (yh - yv)
                u = uv + (uh - uv) * (yl - yv) / (yh - yv)
                #spara värden
                xvec = np.append(xvec, x)
                uvec = np.vstack((uvec, u))
        #Om två koordinater till vänster
        elif np.shape(exv)[0] == 2:
            # xh-koord konstanta för de två interpoleringarna
            xh = exh[0]
            yh = eyh[0]
            uh = euh[0,:]

            
            #interpolera två gånger, varierande xv- värden (två punker)
            for i in range(2):
                xv = exv[i]
                yv = eyv[i]
                uv = euv[i,:]

                x = xv + (xh - xv) * (yl - yv) / (yh - yv)
                u = uv + (uh - uv) * (yl - yv) / (yh - yv)
                #spara värden
                xvec = np.append(xvec, x)
                uvec = np.vstack((uvec, u))

    uvec = uvec[1::,:]
    return xvec, uvec
    
    
    
    