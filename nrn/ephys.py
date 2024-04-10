# ephys.py ---
#
# Filename: ephys.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Fri May  6 14:45:36 2016 (-0400)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Apr 10 16:54:48 2024 (+0530)
#           By: Subhasis Ray
#     Update #: 930

# Commentary:
#
#
#
#

# Code:
"""Utilities for electrophysiology simulation experiments."""
from __future__ import print_function

import sys
import numpy as np
from config import Q_, ur, h
from scipy.optimize import minimize


def calc_em(elist, glist):
    """Takes a list of the reversal potential of the different ions and a
    list of the corresponding ionic conductances.
    """
    return np.dot(elist, glist) / np.sum(glist)


def balance_gbar(em, elist, glist):
    """Balance the gbar of `index`-th entry in glist such that we get em
    as steady state potential"""
    def fun(x):
        ret = np.abs(calc_em(elist, x) - em)
        return ret

    bounds = [(0, None) for ii in glist]
    res = minimize(fun, glist, method='SLSQP', bounds=bounds, tol=1e-6)
    return res


class Mechanism(object):
    """This snippet is modified from that by Andrew Davison.

    http://www.davison.webfactional.com/notes/hoc-to-python-bulbnet/

    """
    def __init__(self, name, **parameters):
        self.name = name
        self.parameters = parameters

    def insert_into(self, section):
        print('A' * 10, 'inserting', self.name, 'into', section.name())
        section.insert(self.name)
        for name, value in self.parameters.items():
            for segment in section:
                mech = getattr(segment, self.name)
                setattr(mech, name, value)


def set_mech_param(section, mech_name, param_name, param_value):
    """Set parameter value for mechanism in every segment in section"""
    for segment in section:
        mech = getattr(segment, mech_name)
        setattr(mech, param_name, param_value)


def create_mechs(param_dict):
    """Create mechanisms from parameter dictionary. It should be a dict of
    dict like: {name: {'gbar': gbar_quantity}} except for passive conductance
    `pas`, which must have {'g': gpas_quantity, 'e': epas_quantity}

    Returns: a dict of {name: ephys.Mechanism}
    """
    ret = {}
    for name, params in param_dict.items():
        if name == 'pas':
            gpas = Q_(params['g'])
            epas = Q_(params['e'])
            mech = Mechanism(name, g=gpas.to('S/cm**2').m,
                             e=epas.to('mV').m)
        else:
            mech = Mechanism(name, gbar=Q_(params['gbar']).to('S/cm**2').m)
        ret[name] = mech
    return ret


def create_cable(name, diameter, length, nseg, RA,
                 CM=Q_('1.0uF/cm**2'), mechs=[], ek=Q_('-80.0mV'),
                 ena=Q_('60.0mV'), eca=Q_('60.0mV'), verbose=False):
    """Create a cable with specified dimeter, length and number of
    segments and insert the mechanisms in mechs list.

    Parameters:

    diameter (um):  cable diameter

    length (um): cable length

    nseg: number of segments to divide the cable into

    RA (ohm-cm): specific axial resistance

    CM (uF/cm2): specific membrane capacitance.

    mechs (Mechanism): list of mechanisms to be inserted into the
    section. You may want to update the properties later using
    `set_mech_param`.

    ek (mV): K+ reversal potential

    ena (mV): Na+ reversal potential

    eca (mV): Ca2+ reversal potential.

    Returns:

    cable: a new Section with all the properties and mechanisms.

    """
    cable = h.Section(name=name)
    cable.diam = diameter.to(ur.um).m
    cable.L = length.to(ur.um).m
    cable.cm = CM.to('uF/cm**2').m
    cable.Ra = RA.to('ohm*cm').m
    cable.nseg = nseg
    print('%' * 100, 'here')
    
    for mech in mechs:
        print('#' * 10, mech)
        mech.insert_into(cable)
    if hasattr(cable, 'ek'):
        cable.ek = ek.to('mV').m
    if hasattr(cable, 'ena'):
        cable.ena = ena.to('mV').m
    if hasattr(cable, 'eca'):
        cable.eca = eca.to('mV').m
    if verbose:
        print('Created cable', name, 'with', nseg, 'segments')
        for x in cable:
            print('Segment at', x.x, 'diam=', x.diam, 'area=', x.area(),
                  'area computed=', np.pi * x.diam * cable.L / cable.nseg,
                  'ri=', x.ri(), 'computed ri=',
                  0.01 * cable.Ra * (cable.L / 2 / cable.nseg) /
                  (np.pi * (x.diam / 2) ** 2))
            for mech in mechs:
                m = getattr(x, mech.name)
                print('  Has mechanism', m.name())
    sys.stdout.flush()
    return cable


def setup_sec_rec(cable, field, t=None):
    """Create and return a list of vectors to record `field` from every segment of
    the cable"""
    vecs = []
    for seg in cable:
        vec = h.Vector()
        # # https://www.neuron.yale.edu/phpBB2/viewtopic.php?f=8&t=1442 tells how nrn counts node pos
        # pos = (0.5 + n) / cable.nseg
        # print('Measuring Vm of segment', n, 'at', pos)
        # if pos > 1.0:
        #     pos = 1.0
        if t is None:
            vec.record(getattr(seg, '_ref_{}'.format(field)))
        else:
            vec.record(getattr(seg, '_ref_{}'.format(field)), t)
    vecs.append(vec)
    return vecs


def setup_mech_rec(cable, mechname, field):
    """Record the specified field of the mechanism `mechname` from all
    segments of the cable"""
    vecs = []
    for n in range(cable.nseg):
        vec = h.Vector()
        # https://www.neuron.yale.edu/phpBB2/viewtopic.php?f=8&t=1442 tells how nrn counts node pos
        pos = (0.5 + n) / cable.nseg
        print('Measuring {}.{} of segment {} at {}'.format(mechname, field, n,
                                                           pos))
        if pos > 1.0:
            pos = 1.0
        mech = getattr(cable(pos), mechname)
        fieldref = getattr(mech, '_ref_{}'.format(field))
        vec.record(fieldref)
        vecs.append(vec)
    return vecs


def setup_current_clamp(cable, delay=100.0 * ur.ms, duration=10.0 * ur.ms,
                        amplitude=100 * ur.nA,
                        pos=0.0):
    """Insert a current clamp electrode to deliver `amplitude` nA current
    at `pos` fraction of length for `duration` ms starting at `delay`
    ms"""
    trode = h.IClamp(cable(pos))
    # avoid conflict with python keyword `del`
    trode.delay = delay.to(ur.ms).m
    trode.dur = duration.to(ur.ms).m
    trode.amp = amplitude.to(ur.nA).m
    print('Inserted current clamp electrode at', pos,
          'fraction of cable to deliver', amplitude,
          'nA current starting at', delay, 'ms for', duration, 'ms')
    return trode


def setup_voltage_clamp(cable, vhold, thold, vpre, tpre, vclamp,
                        tclamp, pos=0.0, rs=Q_('0.001Mohm')):
    """Create a voltage clamp at `pos` location on `cable`.

    vhold: holding voltage (with unit, converted to mV)
    thold: holding period (with unit, converted to ms)
    vpre: prepulse voltage (with unit, converted to mV)
    tpre: prepulse period (with unit, converted to ms)
    vclamp: clamping voltage (with unit, converted to mV)
    tclamp: clamping period (with unit, converted to ms)
    rs: series resistance with clamp circuit  (with unit, converted to Mohm)
    """
    clamp = h.SEClamp(pos, sec=cable)
    clamp.rs = rs.to(ur.Mohm).m
    clamp.amp1 = vhold.to(ur.mV).m
    clamp.dur1 = thold.to(ur.ms).m
    clamp.amp2 = vpre.to(ur.mV).m
    clamp.dur2 = tpre.to(ur.ms).m
    clamp.amp3 = vclamp.to(ur.mV).m
    clamp.dur3 = tclamp.to(ur.ms).m
    return clamp


def update_section(sec, args):
    """Update section properties.

    args: {name: value} for parameters.

    The names `length`, `dia` and `RA` are mapped to L, diam and Ra
    respectively.

    `length`: section.L
    `dia`: section.diam
    `nseg`: section.nseg

    for all mechanisms, the entries should be:

    name: {field: value, ...}

    """
    sec.L = args.pop('length', sec.L)
    sec.diam = args.pop('dia', Q_(sec.diam, 'um')).to('um').m
    sec.Ra = args.pop('RA', Q_(sec.Ra, 'ohm*cm')).to('ohm*cm').m
    sec.nseg = args.pop('nseg', sec.nseg)
    for mech, params in args.items():
        for pname, pvalue in params.items():
            if pname.startswith('g'):
                val = pvalue.to('S/cm**2').m
            elif pname.startswith('e'):
                val = pvalue.to('mV').m
            else:
                raise Exception('Do not know how to handle {} {}'.format(pname, pvalue))
            set_mech_param(sec, mech, pname, val)


#
# ephys.py ends here
