# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 10:09:45 2021

@author: maria
"""

from neuron import h
#from neuron.units import ms, mV
h.load_file('stdrun.hoc')

import matplotlib.pyplot as plt

class Cell:
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0                     # <-- NEW
        h.define_shape()
        self._rotate_z(theta)                            # <-- NEW        
        self._set_position(x, y, z)  
        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
    # everything below here is NEW
    
    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z

    # def _setup_morphology(self):
    #     self.soma = h.Section(name='soma', cell=self)
    #     self.dend = h.Section(name='dend', cell=self)
    #     self.dend.connect(self.soma)
    #     self.soma.L = self.soma.diam = 12.6157
    #     self.dend.L = 200
    #     self.dend.diam = 1

    def _rotate_z(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))

class BallAndStick(Cell):
    name = 'BallAndStick'
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('hh')                                          
        for seg in self.soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')                 
        for seg in self.dend:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -65    # Leak reversal potential mV    

my_cell = BallAndStick(0, 0, 0, 0, 0)
#%%
# my_other_cell = BallAndStick(1)

# del my_other_cell


# h.PlotShape(False).plot(plt)
# ps = h.PlotShape(True)
# ps.show(0)

stim = h.IClamp(my_cell.dend(1))
stim.delay = 5
stim.dur = 1
stim.amp = 0.1

soma_v = h.Vector().record(my_cell.soma(0.5)._ref_v)
dend_v = h.Vector().record(my_cell.dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
h.finitialize(-65) # * mV)

h.continuerun(25) # * ms)



f1 = plt.figure()
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.plot(t, soma_v, linewidth=2)
plt.show(f1)

f = plt.figure()
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
amps = [0.075 * i for i in range(1, 5)]
colors = ['green', 'blue', 'red', 'black']
for amp, color in zip(amps, colors):
    stim.amp = amp
    for my_cell.dend.nseg, width in [(1, 2), (101, 1)]:
        h.finitialize(-65)
        h.continuerun(25)
        plt.plot(t, list(soma_v),
               linewidth=width,
               label='amp=%g' % amp if my_cell.dend.nseg == 1 else None,
               color=color)
        plt.plot(t, list(dend_v), '--',
               linewidth=width,
               color=color)
plt.legend()
plt.show(f)