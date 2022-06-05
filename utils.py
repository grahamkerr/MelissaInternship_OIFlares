import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import copy 
from palettable.colorbrewer.sequential import Blues_9, PuBu_9
from colour import Color
from skimage.exposure import equalize_adapthist, equalize_hist

#############################################################################
#############################################################################
#Graham Kerr
#NASA/GSFC
#Nov 4th 2019
#
#NAME:            mgaussian_model_6term.py
#
#PURPOSE:         Defines an n-component Gaussian 
#
#INPUTS:          p -- the parameters defining the Gaussian(s)
#                      p[0:2] = Amplitude; mean; sigma for 1st comp
#                      p[3:5] = Same for second comp etc., (if present)
#                      p[-1] = background level 
#                      p[-2] = linear comp of background (back = p[-2]*X + p[-1])
#                 x -- the independent variable
#             ncomp -- the number of Gaussian components
#          
#OUTPUTS:        res -- the gaussian model
#
#NOTES:
#
#
def mgaussian_model_6term(p, x, ncomp):
   #-----------------------------------------------------------------------
   # This describes the model and its parameters for which we want to find
   # the best fit. 'p' is a sequence of parameters (array/list/tuple).
   #-----------------------------------------------------------------------

   y = 0.0
   zerolev = p[-1]   # Last element
   for i in range(ncomp):
      A, mu, sigma = p[i*3:(i+1)*3]
      lincomp = p[-2]*x
      y += A * np.exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))
   return y + zerolev + lincomp
          
#############################################################################
#############################################################################
# Graham Kerr
# NASA/GSFC
# 1st July 2019
# 
# NAME:               closest_ind.py
#
# PURPOSE:            To find the index in an array X such that X(ind) is
#                     closest to y
#
# INPUTS:             x -- an array 
#                     y -- the value in x to find
#
# OPTIONAL
# INPUTS:             
#                    
#
# OUTPUTS:            ind -- the index such that X(ind) is closest to y
#                     diff -- the difference between X(ind) and y
#
#
# NOTES:              
#
def closest_ind(x, y):

    ind = (np.abs(x - y)).argmin()

    diff = x[ind] - y

    return ind, diff

#############################################################################
#############################################################################
# Graham Kerr
# NASA/GSFC
# 2nd July 2019
# 
# NAME:               _1gaussian.py
#
# PURPOSE:            A Gaussian fn defined
#
# INPUTS:             
#
# OPTIONAL
# INPUTS:             
#                    
#
# NOTES:              
#
def _1gaussian(x, amp1,cen1,sigma1):

    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2)))

#############################################################################
############################################################################# 

def planck_fn(wvls = [], tg = [], *args):
    
    '''
    Calculates the Planck fn in units of 
    [ergs /s / cm^2 / sr / ang] given wavelength
    in angstrom and temperature in Kelvin. 

    e.g. to produce the Planck Fn at 5800K 
    from 2000 to 10000 A every 1 angstrom:
        wvls = np.arange(2000,10001,1,dtype=float)
        bb = planck_fn(wvls, tg=5800.00)

    Parameters
    __________

    wvls : float
           the wavelength(s)
    tg : float
         gas temperature 
   
    Graham Kerr, Feb 18th 2020

    '''


    # Convert to np array in case it is in input 
    # as a regular list
    wvls = np.array(wvls)
    tg = np.array(tg)

    # Convert to wavelength in cm
    w = wvls / 1.0e8

    # Constants appropriate to cgs units.
    c1 =  3.7417749e-05     # =2*!DPI*h*c*c       
    c2 =  1.4387687e0       # =h*c/k

    bbflux = np.zeros([len(wvls),len(tg)],dtype=float)
    for i in range(len(tg)):
        bbflux[:,i] =  ([c1 / ( x**5 * ( np.exp( c2/tg[i]/x )-1.e0 ) ) for x in w])

    bbint = bbflux*1.0e-8/np.pi

    return bbint 

#############################################################################
############################################################################# 

def nearest(items, pivot):
    '''
    Calculates the nearest value in an array

    Parameters
    __________

    items : the list of items to search over
    pivot : the value to search items for
   
    Graham Kerr, July 2020

    '''


    return min(items, key=lambda x: abs(x - pivot))  

def plotsetup(font_size = 22):

  '''
  Graham Kerr
  NASA/GSFC & CUA
  Feb 2022

  NAME: plotsetup

  PURPOSE: Creates a dictionary with my preferred properties
           for MPL plots. 

  INPUTS: font_size -- flt
                       default = 22

  OUTPUTS: figprops -- a dictionary that is used with rcParams.update()

  NOTES:  
    
    
    '''

  font = {'family': 'Avenir LT Std',
        'color':  'black',
        'weight': 'medium',
        'size': font_size,
        }

  plot_params = {'ytick.direction': 'in', 
               'xtick.direction': 'in', 
               'xtick.minor.visible': True,
               'ytick.minor.visible': True,
               'xtick.major.size': 10, 'xtick.minor.size': 5,
               'ytick.major.size': 10, 'ytick.minor.size': 5,
               'ytick.right': False,
               'xtick.top': False,
               'ytick.left':True,
               'xtick.bottom':False,
               'ytick.major.width': 1.5,
               'xtick.major.width': 1.5,
               'ytick.minor.width': 1.5,
               'xtick.minor.width': 1.5,
               'axes.linewidth': 1.5,
               'axes.spines.top': False,
               'axes.spines.bottom': True,
               'axes.spines.left': True,
               'axes.spines.right': False,
               'axes.titlepad' : 18 }

  plot_lg_params = {'legend.frameon': False}

  plot_dict = {'font.size':font['size'], 
                 'font.family':font['family'], 
                 'font.weight':font['weight'],
                 'ytick.direction': plot_params['ytick.direction'],
                 'xtick.direction': plot_params['xtick.direction'],
                 'xtick.minor.visible': plot_params['xtick.minor.visible'],
                 'ytick.minor.visible': plot_params['ytick.minor.visible'],
                 'ytick.major.size':  plot_params['ytick.major.size'], 
                 'ytick.minor.size':  plot_params['ytick.minor.size'],
                 'xtick.major.size':  plot_params['xtick.major.size'],                                
                 'xtick.minor.size':  plot_params['xtick.minor.size'],
                 'ytick.right': plot_params['ytick.right'],
                 'xtick.top': plot_params['xtick.top'],
                 'ytick.major.width': plot_params['ytick.major.width'],
                 'xtick.major.width': plot_params['xtick.major.width'],
                 'ytick.minor.width': plot_params['ytick.minor.width'],
                 'xtick.minor.width': plot_params['xtick.minor.width'],                    
                 'axes.linewidth': plot_params['axes.linewidth'],
                 'axes.spines.top' : plot_params['axes.spines.top'],
                 'axes.spines.bottom' : plot_params['axes.spines.bottom'],
                 'axes.spines.left' : plot_params['axes.spines.left'],
                 'axes.spines.right' : plot_params['axes.spines.right'],
                 'axes.titlepad' : plot_params['axes.titlepad'],
                 'legend.frameon': plot_lg_params['legend.frameon']}

  return plot_dict

####################################################################
####################################################################

def plotsetup_image(font_size = 22):

  '''
  Graham Kerr
  NASA/GSFC & CUA
  Feb 2022

  NAME: plotsetup_image

  PURPOSE: Creates a dictionary with my preferred properties
           for MPL plots, this version for images (so all spines are included)

  INPUTS: font_size -- flt
                       default = 22

  OUTPUTS: plot_dict -- a dictionary that is used with rcParams.update()

  NOTES:  
    
    
    '''

  font = {'family': 'Avenir LT Std',
        'color':  'black',
        'weight': 'medium',
        'size': font_size,
        }

  plot_params = {'ytick.direction': 'inout', 
               'xtick.direction': 'inout', 
               'xtick.minor.visible': True,
               'ytick.minor.visible': True,
               'xtick.major.size': 10, 'xtick.minor.size': 5,
               'ytick.major.size': 10, 'ytick.minor.size': 5,
               'ytick.right': True,
               'xtick.top': True,
               'ytick.left':True,
               'xtick.bottom':True,
               'ytick.major.width': 1.5,
               'xtick.major.width': 1.5,
               'ytick.minor.width': 1.5,
               'xtick.minor.width': 1.5,
               'axes.linewidth': 1.5,
               'axes.spines.top': True,
               'axes.spines.bottom': True,
               'axes.spines.left': True,
               'axes.spines.right': True,
               'axes.titlepad' : 18 }

  plot_lg_params = {'legend.frameon': False}

  plot_dict = {'font.size':font['size'], 
                 'font.family':font['family'], 
                 'font.weight':font['weight'],
                 'ytick.direction': plot_params['ytick.direction'],
                 'xtick.direction': plot_params['xtick.direction'],
                 'xtick.minor.visible': plot_params['xtick.minor.visible'],
                 'ytick.minor.visible': plot_params['ytick.minor.visible'],
                 'ytick.major.size':  plot_params['ytick.major.size'], 
                 'ytick.minor.size':  plot_params['ytick.minor.size'],
                 'xtick.major.size':  plot_params['xtick.major.size'],                                
                 'xtick.minor.size':  plot_params['xtick.minor.size'],
                 'ytick.right': plot_params['ytick.right'],
                 'xtick.top': plot_params['xtick.top'],
                 'ytick.major.width': plot_params['ytick.major.width'],
                 'xtick.major.width': plot_params['xtick.major.width'],
                 'ytick.minor.width': plot_params['ytick.minor.width'],
                 'xtick.minor.width': plot_params['xtick.minor.width'],                    
                 'axes.linewidth': plot_params['axes.linewidth'],
                 'axes.spines.top' : plot_params['axes.spines.top'],
                 'axes.spines.bottom' : plot_params['axes.spines.bottom'],
                 'axes.spines.left' : plot_params['axes.spines.left'],
                 'axes.spines.right' : plot_params['axes.spines.right'],
                 'axes.titlepad' : plot_params['axes.titlepad'],
                 'legend.frameon': plot_lg_params['legend.frameon']}

  return plot_dict
