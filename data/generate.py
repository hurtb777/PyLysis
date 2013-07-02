# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 12:32:08 2013

@author: BHurt
"""
import numpy as np
from numpy.random import randn
from special import discrete_laguerre


def SISO_with_Laguerres(num_points=1024, memory=40, num_laguerres=3, laguerre_coef=(-0.5, 1., -1.5), alpha=0.5):
    """
    Function used to generate SISO time series data. Default parameters are based off 'test_LET1' in the Lysis7.2 MATLAB
    package.

    :param num_points: Number of data points to generate. len(x) & len(y)=num_points
    :type num_points: Integer, default _1024_
    :param memory: Number of lags in series
    :type memory: Integer, default _40_
    :param num_laguerres: Number of generalized Laguerre functions to use in generating output.
    :type laguerre_coef: Coefficients/weights in front each laguerre to formulate output
    :param alpha: Scaling factor of the generalized Laguerre function
    :type alpha: [-1, 1], default _0.5_

    Returns
    :param x: input data
    :type x: Numpy array
    :param y: output data
    :type y: Numpy array

    """
    #assert isinstance(num_points, (int, float))
    #assert isinstance(memory, (int, float))
    #assert isinstance(num_laguerres, (int, float))

    N = int(num_points)
    M = int(memory)
    L = int(num_laguerres)

    laguerre_coef = laguerre_coef[:num_laguerres] if len(laguerre_coef) >= num_laguerres else np.ones(num_laguerres)
    #lags = range(M)  # generate index memory

    # psuedo gaussian white noise
    x = randn(N)
    dlf = discrete_laguerre(alpha, L, M)
    h = (np.array(laguerre_coef) * dlf).sum(axis=1)  #
    v = np.convolve(x, h)
    y = v[:N] + np.power(v[:N], 2)

    return x,y