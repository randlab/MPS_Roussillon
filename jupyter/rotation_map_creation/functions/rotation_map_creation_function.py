#!/usr/bin/python3
#-*- coding: utf-8 -*-

#07-2019 
#Valentin Dall'alba

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.spatial import distance
from scipy.spatial.distance import squareform, pdist
from scipy.linalg import lu_factor, lu_solve, solve


############################################
############################################

def cloud(x, y, v):
    """
    Computes cloud variogram.
    Arguments:
      x, y : vectors containing the acoordinates of the points
      v : vector containing the values of the measurements at those locations
    Results:
      hc : vector of contains the lag values (distances)
      gc : vector of squared differences of measurements
    """
    X = np.hstack((x[:, np.newaxis], y[:, np.newaxis])
                  )  # alternative: x.reshape((-1,1))  #transforme une liste en colonne, hstack crée une matrice 2*len(x)
    hc = distance.pdist(X) #calcule la distance euclédienne par pair de point
    gc = 0.5 * distance.pdist(v[:, np.newaxis]) ** 2 #return 0.5*diff²
    return (hc, gc)


############################################
############################################

def experimental(hc, gc, lag, nlag):
    """
    Computes experimental variogram from cloud variogram data.
    Arguments:
      hc, gc : the points of the variogram cloud
      lag,nlag : the lag distance and number of lags (number of classes)
    Results:
      he : vector of distances
      ge : vector of values of the experimental variogram
    """

    variogram_range = nlag * lag #must represent your entire domain
    he = np.linspace(lag / 2, variogram_range - lag / 2, nlag) #x axis creation

    # sum of gamma for a given interval
    gamma_cum = np.zeros(nlag)
    
    # number of points in a given interval
    num_points = np.zeros(nlag)
    N = hc.shape[0] #return le nombre de point dans le vecteur distance
    
    # assign values to bins
    for i in np.arange(0, N): #np.arange environ = à range()
        if (hc[i] < variogram_range): #normalement tjs le cas
            if hc[i] < 0:             #normalement tjs le cas
                raise ValueError('Input must be a positive array')
                
            class_index = int((hc[i] / variogram_range) * nlag) #defini dans qu'elle intervalle on place la valeur 
            num_points[class_index] += 1 #on incrémente l'intervalle
            gamma_cum[class_index] += gc[i] #la somme des différence 
    
    #the value of the experimental variogramme
    ge = np.zeros(nlag)
    
    # compute mean, avoid division by 0
    for j in np.arange(0, nlag):
        if num_points[j] != 0:
            ge[j] = gamma_cum[j] / num_points[j] #somme des diff / par le nombre de point dans l'intervalle
        else:
            ge[j] = 0.0
            
    return he, ge


############################################
############################################

def gaussian(h, sill, range):
    """
    Gaussian variogram model
    """
    return sill * (1 - np.exp(-3 * (h / range)**2))


############################################
############################################

def exponential(h, sill, range):
    """
    Exponential variogram model
    """
    return sill * (1 - np.exp(-3 * np.abs(h) / range))


############################################
############################################

def sinus_cardinal(h, sill, range):
    """
    Sinus cardinal variogram model
    """
    return sill * (1 - (range / h) * np.sin(h / range))


############################################
############################################

def hyperbolic(h, sill, range):
    """
    Hyperbolic variogram model
    """
    return sill / (1 + h / range)


############################################
############################################

def nugget(h, nugget):
    """
    Pure nugget variogram model
    """
    return (h == 0) * 0 + (h > 0) * nugget


############################################
############################################

def spherical(h, sill, range):
    """
    Spherical variogram model
    """
    return sill * (3 * h / (2 * range) - 0.5 * (h / range) **
                   3) * (h < range) + sill * (h >= range)


############################################
############################################

def linear(h, sill, range):
    """
    Linear variogram model
    """
    return sill * h / range


############################################
############################################

def stable(h, sill, range):
    """
    Stable variogram model
    """
    return sill * (1 - np.exp(-3 * np.sqrt(h / range)))


############################################
############################################

def ordinary(x, y, v, xi, yi, model_function):
    """
    Ordinary kriging implementation.
    Arguments:
      x,y,v : the data points
      xi,yi : the point where a kriging interpolation is requested
      model_function: variogram model function
    Results:
      v_est : the estimated value at location (xi,yi)
      v_var : kriging variance at location (xi,yi)
    """
    #G*lambda=g
    G = _G_matrix(x, y, model_function)
    g = _g_vector(x, y, xi, yi, model_function)
    
    #on estime les weight du krigeage
    lambda_vec = np.linalg.solve(G, g)
    
    #le krigeage est une somme de points pondérés
    v_est = np.sum(lambda_vec[0:-1] * v)
    v_var = np.sum(-lambda_vec * g)
    
    return v_est, v_var


############################################
############################################

def _G_matrix(x, y, model_function):
    """
    Builds G kriging matrix
    """
    n = x.shape[0]
    X = np.hstack((x[:, np.newaxis], y[:, np.newaxis]))
    G = np.ones((n + 1, n + 1))
    G[0:-1, 0:-1] = -model_function(squareform(pdist(X)))
    G[-1, -1] = 0
    return G


############################################
############################################

def _g_vector(x, y, xi, yi, model_function):
    """
    Builds g kriging vector
    """
    n = x.shape[0]
    g = np.ones(n + 1)
    g[0:-1] = -model_function(np.sqrt((xi - x)**2 + (yi - y)**2))
    return g


############################################
############################################

def ordinary_mesh(x, y, v, xi, yi, model_function):
    """
    Ordinary kriging implementation returning kriging values on a mesh
    Arguments:
      x,y,v : the data points
      xi,yi : the point where a kriging interpolation is requested
      model_function: variogram model function
    Results:
      v_est : array of estimated values at locations (xi,yi)
      v_var : array of kriging variances at locations (xi,yi)
    """
    G = _G_matrix(x, y, model_function)
    lu, piv = lu_factor(G)
    nb_points = xi.shape[0]
    v_est = np.zeros(nb_points)
    v_var = np.zeros(nb_points)
    for i in np.arange(nb_points):
        g = _g_vector(x, y, xi[i], yi[i], model_function)
        lambda_vec = lu_solve((lu, piv), g)
        v_est[i] = np.sum(lambda_vec[0:-1] * v)
        v_var[i] = np.sum(-lambda_vec * g)
    return v_est, v_var


############################################
############################################

def simple(x, y, v, xi, yi, covmodel, mu):
    """
    Simple kriging implementation
    Arguments:
        x, y, v: data points
        xi, yi: point where kriging interpolation is requested
        covmodel: covariance model function
        mu: mean
    Results:
        v_est : array of estimated values at locations (xi,yi)
        v_var : array of kriging variances at locations (xi,yi)
    """
    n = x.shape[0]
    if n==0:
        d = 0;
    else:
        X = np.hstack((x[:, np.newaxis], y[:, np.newaxis]))
        d = squareform( pdist(X))
    C = covmodel(d);
    c = covmodel( np.sqrt( (xi-x)**2 + (yi-y)**2 ) )
    l = solve(C,c)
    
    v_est = np.sum(l*(v-mu)) + mu
    v_var = covmodel(0)-np.sum(l*c)
    return v_est, v_var

