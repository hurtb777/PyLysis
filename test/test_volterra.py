# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:28:16 2013

@author: BHurt
"""

import sys
import numpy as np
sys.path.append('C:\Users\BHurt\Documents\math')
sys.path.append('C:\Users\BHurt\Documents\math\Non_Linear_Dynamics')
from Non_Linear_Dynamics.models.non_parametric.volterra_series  \
    import VolterraLaguerre
    
from Non_Linear_Dynamics.data.generate import SISO_with_Laguerres

X, Y = SISO_with_Laguerres(num_points=1024, num_laguerres=3, memory=40)
v1 = VolterraLaguerre(X, Y, num_laguerres=3)
v1.fit()
#y_est = v1.estimate_y()
v1.plot_kernels()
v1.plot_output()

