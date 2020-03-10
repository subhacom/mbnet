# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 12:37:54 2016

@author: Subhasis
"""
from __future__ import print_function
import numpy as np



def nhpp_thinning(rate_fn, tmax, delta, lbound=None):
    """Nonhomogeneous Poisson process with intensity function `rate_fn` for 
    time range (0, tmax) using the algorithm by Lewis and Shelder 1978.
    
    rate_fn: a function `f(t)` of one variable `t` that returns a finite non negative
    value for `t` in trange.

    tmax: right bound of time (the event times will be for (0, t] time interval)

    delta: interval for evaluating the rate_fn.
    
    lbound: upper bound on lambda. This is used as the rate of the HPP to be 
    thinned. If unspecified then use the maximum value of rate_fn evaluated at 
    points of rate change (0, delta, 2 delta ....)
    
    """
    trange = np.arange(0, tmax+delta, delta)
    if lbound is None:
        lbound = max(rate_fn(t) for t in trange) * 1.0
    st = [0]
    while st[-1] < trange[-1]:
        isi = np.random.exponential(1/lbound, size=len(trange))
        st_ = st[-1] + np.cumsum(isi)
        st += list(st_)
    st = np.array(st[1:])  # remove the dummy 0
    st = st[st <= tmax].copy()    
    if len(st) == 0:
        return np.empty(0)
    accept_prob = np.random.uniform(0, lbound, size=len(st))
    intensity = np.array([rate_fn(t) for t in st])
    return st[accept_prob <= intensity].copy()
        
    
def test_nhpp_thinning():    
    rate_fn = lambda t: 100 * (np.sin(t * np.pi) + 1)
    tmax = 10.0
    delta = 0.1
    plt.subplot(211)        
    spikes = []
    for i in range(10):
        st = nhpp_thinning(rate_fn, tmax, delta)
        plt.plot(st, np.ones_like(st) * (i+1), '|')
        spikes.append(st)
    plt.subplot(212)
    t = np.arange(0, tmax+0.1, 0.1)
    plt.hist(np.concatenate(spikes), bins=t)
    plt.plot(t, rate_fn(t))
    plt.show()


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    test_nhpp_thinning()
