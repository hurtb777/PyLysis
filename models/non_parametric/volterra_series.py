__author__ = 'BHurt'

import numpy as np
import numpy.ma as ma
from numpy.linalg import pinv
from scipy.signal import convolve
import scipy as sp
from scipy.special import genlaguerre
import settings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from special import discrete_laguerre


class Volterra(object):
    """ Volterra Model

        :param x: input data
        :type x: array-like


    """

    def __init__(self, x, y=None, kernel=None, order=2, memory=None, **kwargs):
        self.x = x  # input data
        self.y = y  # actual output values

        if not isinstance(self.x, (list, tuple, np.ndarray, np.matrix)):
            raise ValueError('Input data x must be defined.')
        self.N = len(self.x)  # length of input data (also maximum memory)
        if isinstance(self.x, (list, tuple, np.ndarray, np.matrix)):  # actual output values
            assert self.N == len(self.y)
        self._initialize_kernel(kernel, order, memory)

    def _initialize_kernel(self, kernel=None, order=2, memory=None):
        if kernel:  # check if kernal input is consistent with data
            assert len(kernel) in [1, 2]  # Limit Kernal order to 1st or 2nd order
            self.kernel = kernel
            self.Q = len(self.kernel)  # model order
            self.M = len(self.kernel[1])  # finite memory lag

        else:
            self.Q = order
            self.set_memory(memory)
            self.kernel = []
            self.kernel.append(np.zeros(1))
            for q in range(1, self.Q + 1):
                self.kernel.append(np.zeros(q * [self.M]))
            self.kernel = np.array(self.kernel)

    def estimate_y(self, x=None, y=None):
        """
        Method estimates output values

        """
        if x is not None:
            self.x = x
            self.N = len(self.x)
        if y is not None:
            self.y = y

        dc_offset = self.kernel[0]
        fir_order_conv = convolve(self.kernel[1], self.x)[0:self.N]  # 1-D array of size N
        sec_order_conv = self._calc_second_order_convolution() if self.Q == 2 else np.zeros(self.N)  # 1-D array
        y_est = dc_offset + fir_order_conv + sec_order_conv  # vectorized computation of y_est
        self.y_est = np.array(y_est.T)  # clean up; turn into numpy array consistent with self.x
        self.errors = (self.y.T - self.y_est.T).T
        self.invalid = np.less(range(self.N), self.N * [self.M])
        self.errors = ma.masked_array(self.errors, mask=self.invalid)
        self.rms = np.power(self.errors, 2).mean() / np.mean(self.y)
        return self.y_est, self.errors, self.rms

    def train(self, y=None, order=2, method='LS'):
        """
        Train method estimates the kernals of the system

        :param y: training output. If `None` it looks if self.y is array-like.
        :type y: array-like
        :param order: training output. If `None` it looks if self.y is array-like.
        :type order: array-like
        :param method: training output. If `None` it looks if self.y is array-like.
        :type method: _LS_

        """

        raise NotImplementedError()  # method meant to be overloaded

    def set_memory(self, memory):

        if isinstance(memory, (int, float)):
            self.M = int(memory) if int(memory)<=self.N else self.N
        else:
            self.M = self.N

    def _calc_second_order_convolution(self):
        sec = np.zeros(self.N, dtype=float)
        for n in range(self.M, self.N + 1):
            x_rev = self.x[n - self.M:n][::-1]
            sec[n - 1] = np.dot(np.dot(self.kernel[2], x_rev), x_rev)
        return sec

    def plot(self, kind='x-data', ax=None, **kwargs):
        """
        Basic plotting method for the model capable of plotting the input/output data, output estimations, and kernels.

        :param kind: Data to plot.
        :type kind: string. 'x-data', 'y-data', 'y-est',...
        :param ax: axes to plot data on. If omitted a new axes will be created in a new figure.
        :type ax: _None_ or valid axes instance.

        TODO
        """
        pass

    def plot_data(self):
        """
        TODO
        """
        fig = plt.figure()

        ax1 = plt.subplot(2, 1, 1)
        ax1.plot(self.x, 'b', label='X_1')
        ax1.set_title('Input Data')

        ax2 = plt.subplot(2, 1, 2)
        ax2.plot(self.y, 'r', label='Y_1')
        ax2.set_title('Output Data')
        ax2.set_xlabel('Tick #')

    def plot_kernels(self):
        """
        Plots the 1st and 2nd order kernals

        TODO
        """
        fig = plt.figure()

        # Plot k1
        ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=1)
        ax1.plot(self.kernel[1], '.-b', linewidth=2, label='k_1')
        ax1.set_title('First Order Kernel, $k_1$\nMemory: %d Lags' % self.M)
        ax1.set_xlabel('$m_1$')
        ax1.set_ylabel('$k_1(m_1)$')
        ax1.grid(True)

        # Plot k2
        # TODO: Adjust to an appropriate view of 2nd kernel.
        proj_alpha = 0.5
        subplotspec = GridSpec(3, 1).new_subplotspec((1, 0), rowspan=2)  # creates the larger, lower axis
        ax2 = fig.add_subplot(subplotspec, projection='3d')
        X, Y = np.meshgrid(range(self.M), range(self.M))
        Z = self.kernel[2]
        ax2.plot_surface(X, Y, Z, rstride=2, cstride=2, alpha =0.7)
        ax2.set_zlim(-1.01, 1.01)
        ax2.set_title('Second Order Kernel, $k_2$')
        ax2.set_xlabel('$m_1$')
        ax2.set_ylabel('$m_2$')
        ax2.set_zlabel('$k_2(m_1, m_2)$')
        ax2.set_zlim(-10, 100)
        cset = ax2.contourf(X, Y, Z, zdir='z', offset=2 * np.min(Z), cmap=cm.coolwarm, alpha=proj_alpha)
        #cset = ax2.contourf(X, Y, Z, zdir='x', offset=0, cmap=cm.coolwarm, alpha=proj_alpha)
        #cset = ax2.contourf(X, Y, Z, zdir='y', offset=self.M - 1, cmap=cm.coolwarm, alpha=proj_alpha)

    def plot_output(self, **kwargs):

        # get estimation
        y_est, errors, rms = self.estimate_y()

        # Plot output
        fig = plt.figure()
        ax1 = plt.subplot2grid((2, 1), (0, 0), rowspan=1)
        ax1.plot(self.y, '.-b', linewidth=1, label='$y$', alpha=0.5)
        ax1.plot(self.y_est, '.-g', linewidth=1, label='$y_{est}, m_1=%d$' % self.M, alpha=0.5)
        ax1.set_xlim([0, self.N])
        ax1.set_title('First Order Kernel, $k_1$')
        ax1.set_xlabel('$t$')
        ax1.set_ylabel('$y_1(m_1)$')
        ax1.grid(True)
        ax1.legend()

        # Plot errors
        subplotspec = GridSpec(2, 1).new_subplotspec((1, 0), rowspan=1)  # creates the larger, lower axis
        ax2 = fig.add_subplot(subplotspec)
        ax2.plot(self.errors, 'r',
                 label=r"$\epsilon , RMS = $\displaystyle\frac{\overline{(y-y_{est})^2}}{\overline{y}} = %0.2e$" % rms)
        ax2.set_xlim([0, self.N])
        ax2.set_title('Absolute Estimation Errors')
        ax2.set_xlabel('$t$')
        ax2.set_ylabel('$\epsilon (t)$')
        ax2.grid(True)
        ax2.legend()


class VolterraLaguerre(Volterra):
    """
    VolterraLaguerre is an extension of the Volterra Model. Instead of direct estimation of the kernals, an expansion
     is performed about each kernal using Laguerre orthogonal basis functions.

    :param x: training input.
    :type x: array-like
    :param kernal: training output. If `None` it looks if self.y is array-like.
    :type kernal: array-like

    """

    def __init__(self, x, y=None, kernel=None, order=2, alpha=0.5, memory='from_alpha', num_laguerres=5, **kwargs):
        Volterra.__init__(self, x, y, kernel, order, memory, **kwargs)
        self.set_memory_from_alpha(alpha) if memory == 'from_alpha' else None
        self.alpha = alpha
        self.num_laguerres = num_laguerres

    def set_memory_from_alpha(self, alpha):
        self.M = int(np.ceil((-30 - np.log(1 - alpha)) / np.log(alpha)))
        self.alpha = alpha

    def train(self, y=None, order=2, method='LS'):
        """
        Train method estimates the kernals of the system using orthogonal laguerre basis functions. Training is
         carried out with the following:

        * L basis functions are generated and evaluated - yields B matrix
        * Convolve each ith Laguerre function, b_i, with input data x(t).  - yields v_i(t) and V matrix
        * Solve Y = C*V, where C is the coefficient matrix of the Laguerre functions. Solve using `method`
        * Calculate kernels:
            * k_0 = c_0
            * k_1(m) = c_1* b(m)
            * k_2(m_1, m_2) = b(m)^T * c_2 * b(m)

        :param y: training output. If `None` it looks if self.y is array-like.
        :type y: array-like
        :param order: training output. If `None` it looks if self.y is array-like.
        :type order: array-like
        :param method: Method used to find the Laguerre coefficients, C
        :type method: _LS_, _(more to be added)_
        """

        try:
            self.y = y if y is not None else self.y
            assert len(self.y) == len(self.x)
        except TypeError():
            print "output y must be defined and of length len(x)"

        # setup volterra matrix: B
        self._eval_laguerre_fxns()
        # setup convolved voltera matrix with data: conv(B,x)V
        self._eval_laguerre_conv()
        # perform regression to estimate the coefficients C for each laguerre
        self._est_laguerre_coef()
        # calc kernels
        self._calculate_kernels()

    def _eval_laguerre_fxns(self):
        """
        Evaluate matrix B
         Rows: Memory size (self.M)
         Cols: Number of Laguerre Fxns (self.num_laguerres)
        """
        self.B = discrete_laguerre(self.alpha, self.num_laguerres, self.M)

    def _eval_laguerre_conv(self):
        """
        Evaluate matrix V
        Rows: Length of input data (self.N)
        Cols: Number of Laguerre Fxns (self.num_laguerres)
        """
        V = np.zeros((self.N, self.num_laguerres))

        for l in range(self.num_laguerres):
            V[:, l] = np.convolve(self.B[:, l], self.x)[0:self.N]
        self.V = V

    def _est_laguerre_coef(self):
        """
        Method used to estimate the coefficients in front of each laguerre convolution, v

        """
        num_rows = self.N
        L = self.num_laguerres
        num_cols = 1 + L + L * L  # mod cols: 1+L+L(L+1)/2
        V = np.matrix(np.zeros((num_rows, num_cols)))  # modified V matrix to solve Y=V*C
        V[:, 0] = 1
        V[:, 1:L + 1] = self.V
        col = L
        for i in range(L):
            for j in range(L):
                col += 1
                V[:, col] = np.matrix(self.V[:, i] * self.V[:, j]).T

        self.VV = V

        self.VVinv = pinv(V.T * V)
        self.Cc1 = self.VVinv * V.T
        self.Cc2 = self.Cc1 * np.matrix(self.y).T
        C_ = self.Cc2

        self.C = [C_[0]]
        self.C.append(C_[1:L + 1])
        #print L
        #print C_.shape
        #print len(C_[L + 1:])
        #for i in range(L):
        #    for j in range(L):

        self.C.append(C_[L + 1:].reshape((L , L)))

    def _calculate_kernels(self):
        """
        Kernel calculation using the Laguerre functions and its coefficients.

        Populates the kernel attribute. self.kernel.size = (kernel order,) for the k0, k1, k2, ..., kQ.
        Currently `Q=1` or `Q=2` are only supported.
        """
        self.kernel[0] = self.C[0]
        self.kernel[1] = np.squeeze(np.asarray(np.dot(self.C[1].T, self.B.T)))  # matrix multiplication
        self.kernel[2] = np.zeros((self.M, self.M))
        for i in range(self.M):
            for j in range(self.M):
                self.kernel[2][i, j] = np.matrix(self.B[i, :]) * self.C[2] * np.matrix(self.B[j, :]).T
