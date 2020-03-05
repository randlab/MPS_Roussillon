#!/usr/bin/python3
#-*- coding: utf-8 -*-

#2019 
#Valentin Dall'alba

###################
###################

import numpy as np
import pickle
import pandas as pd

from geone import img
import geone.imgplot as imgplt
import geone.customcolors as ccol
import geone.deesseinterface as dsi



def create3DGrid(topLayer, bottomLayer, pathGSLIB):
    '''
    Function that creates a 3D grid based on top and bottom img values.
    Inputs : 
    -------------
    topLayer : top topography of the grid, Img object.
    bottomLayer : bottom topography of the grid, Img object.
    pathGSLIB : path where to store the gslib/pickle grid file.
    
    Outputs : 
    -------------
    One Img object with two variables (the grid and the transformed grid).
    One Img object with one variable (the shift value for each cell).
    '''
    
    nx, ny     = topLayer.nx, topLayer.ny
    ox, oy     = topLayer.ox, topLayer.oy
    sx, sy, sz = topLayer.sx, topLayer.sy, topLayer.sz
    
    minDepth = np.nanmin(bottomLayer.val)
    maxDepth = np.nanmax(topLayer.val)
    nz       = np.ceil((maxDepth-minDepth)/sz)
    
    grid3D = img.Img(nx=int(nx), ny=int(ny), nz=int(nz),
                  sx=sx, sy=sx, sz=sz,
                  ox=ox, oy=oy, oz=minDepth,val='nan',
                  nv=2,varname=['Pliocene','Transformation'])
    
    grid2D_transfoInfo = img.Img(nx=int(nx), ny=int(ny), nz=1,
                             sx=sx, sy=sy, sz=sz,
                             ox=ox, oy=oy, oz=0,val='nan',
                             nv=1,varname='nb_layer_transfo')
    oz = grid3D.oz
    
    for i in range(ny):
        for j in range(nx):
        #We only take cell where both top and bottom layers are informed.
        #We also want that top>bottom.
        
            if str(bottomLayer.val[0,0,i,j])!='nan' and topLayer.val[0,0,i,j]!='nan' and topLayer.val[0,0,i,j]-bottomLayer.val[0,0,i,j]>=0: 
                #We calculate the index of the bottom and top layer of the 3D grid,
                mur = bottomLayer.val[0,0,i,j]
                #murLayer = int((nz+mur)/sz)
                #murLayer = int(mur-(oz*sz))
                
                toit = topLayer.val[0,0,i,j]
                #toitLayer = int((nz+toit)/sz)
                #toitLayer = int(toit-(oz*sz))
                
                thickness = toit-mur
                thicknessLayer = int(np.ceil(thickness/sz)) #nb of layer thick
                murLayer =int(np.trunc((mur-oz)/sz)) #bottom layer 
                
                #We assigne the one value of the cells composing the 3D grid
                grid3D.val[0,murLayer:(murLayer+thicknessLayer+1),i,j] = 1
                
                #We assigne the transform value to the 3D grid
                grid3D.val[1,0:int(thicknessLayer+1),i,j] = 1
                
                #The murLayer info corresponds to the number of layer of which the cell has been transpose
                grid2D_transfoInfo.val[0,0,i,j] = murLayer
                
    
    with open(pathGSLIB+'grid3D.pickle','bw') as file:
        pickle.dump(grid3D, file, pickle.HIGHEST_PROTOCOL)
        
    with open(pathGSLIB+'grid3D_info.pickle','bw') as file:
        pickle.dump(grid2D_transfoInfo, file, pickle.HIGHEST_PROTOCOL)
    
    return grid3D, grid2D_transfoInfo