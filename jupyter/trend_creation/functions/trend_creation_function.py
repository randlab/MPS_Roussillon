#!/usr/bin/python3
#-*- coding: utf-8 -*-


########
#The following functions are used to create simply diffusive map, where boundary condition can be specify per zones,
#The functions have to be used in their presented order.
#Valentin Dall'alba 2020
########

import flopy
from flopy.utils.triangle import Triangle as Triangle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import geopandas as gp
import shutil
from scipy.interpolate import griddata
from geone import img
import geone.imgplot as imgplt
import geone.customcolors as ccol
import geone.deesseinterface as dsi
import matplotlib.ticker as plticker


########
#1
########
def create_grid(path_shp_file, tolerance=1, hull=False, plot=True):
    '''
    Function that create the grid from a shp file.
    The tolerance value allow to simplify the grid geometry.
    The hull parameter decide if we wand the geometry to be create by a convex hull polygone.
    The plot parameter controls if the created grid is plotted.
    
    Inputs : 
    -----------
    path_shp_file : path to the shape file.
    tolerance : >=1 to simplify the geometry
    hull : convex hull polygone
    plot : return a plot of the create grid (boolean).
    
    Outputs:
    -----------
    layer : geopanda layer geometry
    '''
    
    
    layer_shp = gp.read_file(path_shp_file)
    
    if hull is False:
        layer = layer_shp.geometry.simplify(tolerance=tolerance)
    else:
        layer = layer_shp.geometry.convex_hull
        
    if plot==True:
        layer.plot()
    
    return layer


########
#2
########
def create_mesh(layer, max_area=10000000, max_angle=30, saveFig=False):
    '''
    Function that create a triangular mesh based on a layer input.
    The created mesh is then plot.
    The max_area and max_angle parameter control the shape of the triangular shape.
    The saveFig parameter if True save the create mesh to pdf.
    The path to the triangular.exe must be change inside the function regarding the location.
    
    Inputs : 
    -----------
    layer : geopanda layer geometry, created with the create_grid function
    max_area : maximum area of a cell
    max_angle : maximum angle of a cell
    saveFig : store the to pdf the created mesh (boolean).
    
    Outputs : 
    -----------
    layer_pts : coordinates of the points geometry.
    mesh : triangular mesh    
    '''
    
    #Extract the points that compose the geometry
    layer_pts = []
    for x,y in zip(layer.geometry[0].boundary.xy[0],layer.geometry[0].boundary.xy[1]):
        layer_pts.append((float(x),float(y)))
    
    #Create a directory to store the mesh/mf6 files
    path_ws = './mesh'
    if os.path.exists(path_ws):
        shutil.rmtree(path_ws)
    os.makedirs(path_ws)
    
    #Build the mesh
    path_tri = './linux_bin/triangle'
    mesh     = Triangle(maximum_area=max_area, angle=max_angle,
                   model_ws=path_ws, 
                   exe_name=path_tri)
    mesh.add_polygon(layer_pts)
    mesh.build()
    
    #Plot the mesh
    fig = plt.figure(figsize=(5,5))
    ax  = plt.subplot(1, 1, 1, aspect='equal')
    pc  = mesh.plot()
    
    start, end = ax.get_xlim()
    locx = plticker.MultipleLocator(base=(end-start)/15) # this locator puts ticks at regular intervals
    start, end = ax.get_ylim()
    locy = plticker.MultipleLocator(base=(end-start)/20) # this locator puts ticks at regular intervals

    ax.xaxis.set_major_locator(locx)
    ax.yaxis.set_major_locator(locy)
    
    plt.xticks(rotation=90)
    ax.grid(c='lightblue',which='minor', alpha=0.8)
    ax.grid(c='tan',which='major')
    plt.tight_layout()
    
    if saveFig==True:
        fig.savefig('mesh.pdf') 
        
    return layer_pts, mesh


########
#3
########
def define_cst_heads(layer_pts, mesh, zones_list=False):
    '''
    Create the constant head list base on the zones coordinnates passed as input.
    
    Inputs : 
    -----------
    layer_pts : coordinates of the points geometry, created with the function create_mesh.
    mesh : triangular mesh, created with the function create_mesh.
    zones_list : coordinnate of the zone, can be defined manually by the user if False.
    
    Outputs:
    -----------
    chdlist : constant head list for the edge cells.
    
    '''
    
    #Get the edge cells Id
    edgenodes = []
    for iedge in range(len(layer_pts)):
        nodes = mesh.get_edge_cells(iedge)
        for n in nodes:
            if n not in edgenodes:
                edgenodes.append(n)
    
    cells       = mesh.get_cell2d()
    cst_heads   = []            
    
    #Read the coordinnate of the zones
    if zones_list is not False:
        nb_cst_head = len(zones_list)
        
        for i in range(len(zones_list)):
            cst_head = []
            for zone in zones_list[i]:
                xMin, xMax = zone[0], zone[1]
                yMin, yMax = zone[2], zone[3]
            
                for node in edgenodes:
                    if cells[node][1]>=xMin and cells[node][1]<xMax and cells[node][2]>=yMin and cells[node][2]<yMax:
                        cst_head.append(node)
                        
            cst_heads.append(cst_head)
            
    #Manually define the zones
    else: 
        #Define the cells cst head value
        nb_cst_head = int(input('How many cst head groups ? :'))

        for i in range(nb_cst_head):
            print('Define zone for cst head group {} :'.format(i))
            zones = int(input('How many zones ? :'))
            cst_head = []
        
            for zone in range(zones):
                print('Please define the coordinnate of the zone :')
                xMin, xMax = float(input('xMin :')), float(input('xMax :'))
                yMin, yMax = float(input('yMin :')), float(input('yMax :'))
            
                for node in edgenodes:
                    if cells[node][1]>=xMin and cells[node][1]<xMax and cells[node][2]>=yMin and cells[node][2]<yMax:
                        cst_head.append(node)
                    
            cst_heads.append(cst_head)
            
    #Print how many edge cells have been attribuate a value       
    for i in range(len(cst_heads)):
        print(len(cst_heads[i]))

    #We create the cst head list
    chdlist = []
    for i in range(nb_cst_head):
        cst_hd_value = float(input('What is the value of the cst head group {} ?'.format(i)))
        for icpl in edgenodes:
            if icpl in cst_heads[i]:
                chdlist.append([(0, icpl), cst_hd_value])
    
    return chdlist


########
#4
########
def run_simulation(mesh, chdlist):
    '''
    Run the simulation.
    The path to the mf6 exe file must be change inside the function.
    
    Inputs : 
    -----------
    mesh : triangular mash, created with the function create_mesh.
    chdlist : constant head cell list, created with the function define_cst_heads.
    
    Outputs :
    -----------
    The simulation are automatically stored in a simulation_mf6 folder.
    '''
    
    #Create the mesh element
    cell2d   = mesh.get_cell2d()
    vertices = mesh.get_vertices()
    xcyc     = mesh.get_xcyc()
    ncpl     = mesh.ncpl
    nvert    = mesh.nvert
    nlay,hk  = 1, 0.8
    top,botm = 1, [0.]
    name     = 'mf'
    path     = './mesh/'
    exe_path = './linux_bin/mf6'
    
    #Create the input of the Modflow model
    sim  = flopy.mf6.MFSimulation(sim_name=name, version='mf6',
                             exe_name=exe_path,
                             sim_ws='simulation_mf6')
    
    tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS',
                             perioddata=[[1.0, 1, 1.]])
    
    gwf  = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    
    ims  = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='complex', 
                           outer_hclose=1.e-5, inner_hclose=1.e-4)
    
    dis  = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, nvert=nvert,
                               top=top, botm=botm,
                               vertices=vertices, cell2d=cell2d)
    
    npf  = flopy.mf6.ModflowGwfnpf(gwf, k=hk,xt3doptions=[True], 
                              save_specific_discharge=True)
    
    ic   = flopy.mf6.ModflowGwfic(gwf)
    
    chd  = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdlist)
    
    rch  = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0)               
    
    oc   = flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='{}.cbc'.format(name),
                            head_filerecord='{}.hds'.format(name),
                            saverecord=[('HEAD', 'LAST'),
                                        ('BUDGET', 'LAST')],
                            printrecord=[('HEAD', 'LAST'),
                                         ('BUDGET', 'LAST')])
    
    #Run the simulation
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(report=True, silent=True)
    print('Simulation is a success ? : {}'.format(success))
    
    return


########
#5
########
def get_head(path='./simulation_mf6/',plotFig =False, saveFig=False):
    '''
    Get the head output of the simulation.
    The head output can be plot or save as pdf.
    
    Inputs:
    ---------
    path : path ot the simulation folder.
    plotFig : (boolean) 
    saveFig : (boolean)
    
    Outputs:
    ---------
    head : simulated head valus.
    '''
      
    fname = os.path.join(path, 'mf' + '.hds')                #Create the path to get the head file
    hdobj = flopy.utils.HeadFile(fname, precision='double')  #Read the head file output
    head  = hdobj.get_data()                                 #Extract an array of the head outputs
    
    if plotFig == True:
        #Plot the head
        fig = plt.figure(figsize=(5,5))
        ax  = plt.subplot(1, 1, 1, aspect='equal')
        h   = grid.plot(ax=ax, a=head[0, 0, :], cmap='winter', alpha=.9)
        
    if saveFig ==True:
        fig.savefig('head_output.pdf')
    
    return head


########
#6
########
def mf_to_geone(mesh, head, mask=False, xMin=664328.1865, xMax=664737.1865,
                                        yMin=6153000.2413, yMax=6153512.2413,
                                                           nx=409, ny=512, sx=100, sy=100):
    '''
    Interpolate the head output to a cartesian grid.
    The interpolate grid is then transformed to an Img object.
    
    Inputs:
    ---------
    mesh : triangular mesh, created with the create_mesh function.
    head : head value, created with the get_head function.
    
    Outputs:
    ---------
    trend : Geone Img of the head values interpolate on the simulation grid.
    '''
    if mask != False:
        nx, ny     = mask.nx, mask.ny
        sx, sy     = mask.sx, mask.sy 
        xMin, yMin = mask.ox, mask.oy
        xMax, yMax = xMin+(sx*nx), yMin+(sy*ny)
        
    #Get the coordinate of the grid
    x, y = mesh.get_xcyc()[:,0],mesh.get_xcyc()[:,1]
    head_val = head[0,0,:]
    
    #Create the new grid 
    xG   = np.linspace(xMin ,xMax, nx)
    yG   = np.linspace(yMin ,yMax, ny)
    X, Y = np.meshgrid(xG, yG)

    #Interpolate the values
    interpo = griddata((x,y),head_val,(X,Y),method='cubic')

    #Mask values
    trend = interpo[np.newaxis,np.newaxis,:]
    if mask!= False:
        trend[mask.val!=1]=np.nan
    
    #Create Img Geone
    trend = img.Img(nx=nx,ny=ny,nz=1, 
                sx=sx,sy=sy,sz=1, 
                ox=xMin,oy=yMin,oz=0,
                nv=1,val=trend)
    
    return trend