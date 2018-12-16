# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:21:03 2018

Plots guassian function used in simulations.

@author: fionnlagh
"""
import numpy as np
import matplotlib.pyplot as plt
import math


def erf(z):
    t = 1.0 / (1.0 + 0.5 * abs(z))
    # use Horner's method
    ans = 1-t*math.exp(-z*z-1.26551223 +
                       t*(1.00002368 +
                          t*(0.37409196 +
                             t*(0.09678418 +
                                t*(-0.18628806 +
                                   t*(0.27886807 +
                                      t*(-1.13520398 +
                                         t*(1.48851587 +
                                            t*(-0.82215223 +
                                               t*(0.17087277))))))))))
    if z >= 0.0:
        return ans
    else:
        return -ans


def theta(x):
    theta = 0.5*(1.0+erf(x))
    return theta

T = True
F = False

skewed = F
norm = F
tanh = T

fontsize = 18
linewidth = 5
markersize = 200

cm_2_km = 1.0
A = 9
Unit_l = 1
jet_w = 0.05/cm_2_km/Unit_l  # no units
jet_plot = [-jet_w, jet_w]

# x information
domain_nx1 = 40  # number of grid pts
AMR_lvl = 4
domain_nx1_max = domain_nx1*(2**AMR_lvl)  # max res
xprobmin1 = -1.5
xprobmax1 = 1.5
x0 = 0.0
# @ Contians jet in guass dist as 3 sigma guass dist is zero.
deltax = jet_w*2

# y information
domain_nx2 = 20
xprobmin2 = 0.0
xprobmax2 = 1.5
y0 = 0.0  # (xprobmax2-xprobmin2)/2
deltay = (xprobmax2-xprobmin2)/domain_nx2
# create grid points
y = np.linspace(xprobmin2, xprobmax2, domain_nx2)
x = np.linspace(xprobmin1, xprobmax1, domain_nx1_max)

rhobc = 10

if norm:
    dist_test_x = A*np.exp(-((x-x0)/deltax)**2)
    dist_test_y = A*np.exp(-((y-y0)/deltay)**2)
# endif

if skewed:
    alpha = 80.0
    dist_test_x = np.exp(-((x-x0)/deltax)**2)
    theta_dist = np.zeros(len(x))
    for i in range(len(x)):
        theta_dist[i] = theta(x[i]*alpha)
    thi = (1/np.sqrt(2*np.pi))*dist_test_x
    skewed_dis = 2*thi*theta_dist
    dist_test_x = skewed_dis
# endif

if tanh:
    trwidth = jet_w  # transition width
    PoI = 0.5*(-jet_w+(-jet_w+trwidth))  # point of inflection
    smoothness = trwidth*0.2
    dist_test_x = []
    dist_x_pts = []
    for i in range(len(x)):
        if x[i] > -jet_w and x[i] < 0:
            dist_test_x.append(A*0.5*(1+np.tanh((x[i]-PoI)/smoothness)))
            dist_x_pts.append(x[i])
        if x[i] < jet_w and x[i] > 0:
            dist_test_x.append(A*0.5*(1-np.tanh((x[i]+PoI)/smoothness)))
            dist_x_pts.append(x[i])
        #endif
    #endfor
#endif


if skewed or norm or tanh:
#    plt.xlim(-0.2, 0.2)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.xlabel('x [Mm]', fontsize=fontsize)
    plt.ylabel('Velocity [km s-1]', fontsize=fontsize)
    plt.scatter(jet_plot, np.zeros(len(jet_plot)), alpha=0.9, s=markersize,
                zorder=1, color='r')
    if tanh:
        plt.plot(dist_x_pts, dist_test_x, linewidth=linewidth, zorder=-1)
    else:
        plt.plot(x, dist_test_x, linewidth=linewidth, zorder=-1)
    plt.grid()
    plt.show()
