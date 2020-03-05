#!/usr/bin/python3
#-*- coding: utf-8 -*-

#2019 
#Valentin Dall'alba
from geone import img
import geone.imgplot as imgplt
import geone.customcolors as ccol


#################
#################

def txtToGslib_GIS(pathTXT, pathGSLIB, nz=1, sx=100, sy=100, sz=2, nanV=-999):
    '''
    Function that convert a ASCII raster file from QGIS to Gslib file.
    Inputs : 
    -------------
    pathTXT : path to the ascii file to transform.
    pathGSLIB : path where to store the convert file.
    Dimensions : sx, sy ,sz dimension of the Img cells size to create.
    nanV : no value of the ascii file.
    
    Outputs :
    -------------
    Img object.
    '''
    
    #read the ascii
    with open(pathTXT,'r')as textASCII:
        lines = textASCII.readlines()
        oy = (lines[1].split()[1])
        ox = (lines[3].split()[1])
        nx = int(lines[5].split()[1])
        ny = int(lines[4].split()[1])
        dataASCII = lines[6:]
        data = []
        dataASCII.reverse()
        for line in dataASCII:
            listData = line.split()
            for elt in listData:
                data.append(elt)
                
    #write the gslib
    with open(pathGSLIB,'w') as textGSLIB:
        textGSLIB.write('{} {} {} {} {} {} {} {} 0\n'.format(nx,ny,nz,sx,sy,sz,ox,oy))
        textGSLIB.write('1\n')
        textGSLIB.write('altitude\n')
    
        for value in data:
            if value==str(nanV):
                textGSLIB.write('nan')
            else:
                textGSLIB.write(value)
            textGSLIB.write('\n')    
            
    #read the gslib to Img        
    imageRead = img.readImageGslib(pathGSLIB)
    
    print('*** txt to Gslib Done***')
    
    return imageRead