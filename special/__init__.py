__author__ = 'BHurt'

"""
========================================
Special functions (:mod:`scipy.special`)
========================================

.. module:: scipy.special

Nearly all of the functions below are universal functions and follow
broadcasting and automatic array-looping rules. Exceptions are noted.

Error handling
==============

Errors are handled by returning nans, or other appropriate values.
Some of the special function routines will emit warnings when an error
occurs.  By default this is disabled.  To enable such messages use
``errprint(1)``, and to disable such messages use ``errprint(0)``.

Example:

    >>> print scipy.special.bdtr(-1,10,0.3)
    >>> scipy.special.errprint(1)
    >>> print scipy.special.bdtr(-1,10,0.3)

.. autosummary::
   :toctree: generated/

   errprint

Available functions
===================

Airy functions
--------------

.. autosummary::
   :toctree: generated/

   airy     -- Airy functions and their derivatives.
   airye    -- Exponentially scaled Airy functions
   ai_zeros -- [+]Zeros of Airy functions Ai(x) and Ai'(x)
   bi_zeros -- [+]Zeros of Airy functions Bi(x) and Bi'(x)
"""

import numpy as np
from specfun import *

__all__ = ['discrete_laguerre']