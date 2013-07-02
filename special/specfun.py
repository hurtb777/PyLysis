__author__ = 'BHurt'
import numpy as np


def discrete_laguerre(alpha=0.5, L=5, M=1):
    """
    Computes and evaluates the discrete laguerre basis functions.

    :param alpha: Scaling parameter of the laguerre functions
    :type alpha: number between 0 and 1 (_0.5_)
    :param L: Number of Laguerre functions to compute and evaluate
    :type L: Greater than 0 integer (_5_)
    :param M: Number of integer points along the semi-infinite line to evaluate (starting from 0)
    :type M: Greater than 1 integer (_1_)

    Returns:
        2-D array where the primary index (row) is the memory and the secondary index (column) are the laguerres

        >>> discrete_laguerre()

    """
    assert L >= 1
    assert M >= 1
    assert alpha > 0.
    assert alpha < 1.

    m = np.arange(0, M)  # sequence of numbers from 0 to M-1
    beta = 1 - alpha
    ralpha = np.sqrt(alpha)

    dlf = np.zeros((M, L))
    temp = np.zeros(M)
    L_0_N = np.sqrt(np.power(alpha, m) * beta)  # 0th order Laguerre
    dlf[:, 0] = L_0_N
    if L == 0:
        return dlf
    for l in range(1, L):
        for n in m:
            temp[n] = L_0_N[n]
            if n == 0:
                L_0_N[n] = ralpha * L_0_N[n]
            else:
                L_0_N[n] = ralpha * (L_0_N[n - 1] + temp[n]) - temp[n - 1]
        dlf[:, l] = L_0_N.T

    return dlf




