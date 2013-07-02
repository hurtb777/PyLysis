.. _outline:

**************
About pyLysis
**************

PyLysis is Python module for non-parametric non-linear modeling using Volterra-Laguerre Series and Principal Dynamic Modes (PDM). 

Introduction
=============

A Volterra series, similar to a Taylor series expansion in structure, has the added benefit of capturing 'memory' effects. However direct estimation of the higher order kernels is challenging for a pure Volterra series because kernel coefficients generally covary with one another so traditional regression techniques do poorly in parameter estimations. One way to simplify kernel estimation is through expanding them on a orthogonal Laguerre basis. 

Once the kernels are estimated, PDM analysis can be used to characterize and explain the dynamics of the underlying system. PDM analysis has been used in neuroscience and electrical engineering because of the timeseries-like nature of the measurements in these fields. 

.. seealso::

   Marmarelis, V. Z.
   Nonlinear Dynamic Modeling of Physiological Systems
   Wiley-IEEE Press, 2004.

Model Structure
================

In discrete time, the general input :math:`x(n)` - output :math:`y(n)` relation of a stable (finite-memory) nonlinear time-invariant dynamic system is given by the discrete-time Volterra series

.. math::

    y(n) = k^{(0)} + \sum_{m_1=0}^{M} k^{(1)}_{m_1}x(n-m_1) + \sum_{m_1=0}^{M_1}\sum_{m_2=0}^{M_2} k^{(2)}_{m_1, m_2}x(n-m_1)x(n-m_2)+\cdots
..

where :math:`k^{(n)}` are known as the Volterra kernels which describe the dynamics of the system at each order of nonlinearity, and constitute a complete and canonical representation of the nonlinear system. Estimating these kernels in practice is challenging because of the exponential scaling of parameters as the order of the model grows (for a model of order :math:`O` and finite memory :math:`M` the number of parameters to estimate are :math:`\sum_{i=0}^{O} M^i` (:math:`M\sim 50`). Furthermore the kernel coefficients can be estimated in practice by using linear regression, however in general kernel coefficients co-vary so further mathematical treatment is necessary.


A useful practice is to transform the model into an orthogonal projection (which assures zero-covariance between parameters). Wiener suggested expanding the kernels by using :math:`L` orthogonal using Laguerre Basis functions :math:`\left\{ b_{j}(m) \right\}` which also have a built-in exponential term making them suitable for physical systems. An added benefit from a statistical perspective is that this transformation reduces the number of free paramters significantly( :math:`\sum_{i=0}^{O} L^i`, where :math:`L<10` ). For this module a Volterra-Laguerre series is used which is given by

.. math::

    y(n) = c^{(0)} + \sum_{j_1=0}^{L} c^{(1)}_{j_1}v_{j_1}(n) + \sum_{j_1=0}^{L_1}\sum_{j_2=0}^{L_2} c^{(2)}_{j_1,j_2}v_{j_1}(n)v_{j_2}(n)+\cdots
..

where

.. math::

    v_{j}(n)=\sum_m b_j(m)x(n-m)
..

and :math:`c^{(1)}_j, c^{(2)}_{j_1,j_2}, \cdots` represent the expansion coefficients of the respective kernels. Estimating the coefficients :math:`c_i` can be done using linear regression or other statistical estimation methods. Recapturing the kernels from these coefficients is then straightforward.

Principal Dynamic Modes
========================

Principal Dynamic Mode (PDM) analysis is a non-parametric framework to estimate a system's nonlinear dynamics. Built on top of a Volterra series expansion, PDM analysis is a way to segregate the linear from the nonlinear dynamics of a system. Performing such an analysis opens a new window to view the structure of the system - something traditional kernel & non-parametric approaches typically lack.

Essentially a signal processing technique, PDM analysis uses a Volterra-Laguerre expansion framework to model the system's output.

*******************
Installing PyLysis
*******************

TODO

Dependences
===========

PyLysis depends on

* Python 2.7
* Numpy 1.6.2
* Scipy 0.11
* matplotlib 1.1.0

Download PyLysis
==================

TODO

Licence
=========

PyLysis is licensed under the terms of the GNU General Public License as
published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Click `here <#>`_ to download PyLysis.

PyLysis would like to be a collaborative project. Any help is welcome, from
a bug notification to a more demanding collaboration. Please refer to my
`webpage <#>`_ for contacts.

