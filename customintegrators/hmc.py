#!/usr/bin/env python

"""
Hybrid Monte Carlo integrator implemented using CustomIntegrator

"""

from simtk import openmm, unit
from openmmtools import testsystems

timestep = 1.0 * unit.femtoseconds
temperature = 300.0 * unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature

# Integrator initialization
integrator = openmm.CustomIntegrator(timestep)
integrator.addPerDofVariable("xold", 0) # old positions
integrator.addGlobalVariable("Eold", 0) # old total energy
integrator.addPerDofVariable("x1", 0) # positions prior to constraint correction
integrator.addPerDofVariable("sigma", 0) # Gaussian standard deviation for Maxwell-Boltzmann distribution
integrator.addGlobalVariable("accept", 0) # determines whether move is accepted or not
integrator.addGlobalVariable("kT", kT.value_in_unit_system(unit.md_unit_system)) # thermal energy

# Draw a new velocity from the Maxwell-Boltzmann distribution
integrator.addComputePerDof("v", "sigma*gaussian; sigma = sqrt(kT/m)")
integrator.addConstrainVelocities() # constrain velocities (RATTLE)

integrator.addComputeSum("ke", "0.5*m*v*v") # compute kinetic energy
integrator.addComputeGlobal("Eold", "ke + energy") # compute old total energy
integrator.addComputePerDof("xold", "x") # store old positions

# Velocity Verlet step
integrator.addUpdateContextState() # allow other forces (e.g. MonteCarloBarostat) to update system
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("x", "x+dt*v") # position full kick
integrator.addComputePerDof("x1", "x") # store positions before position constraints are applied
integrator.addConstrainPositions() # apply position constraints (CCMA)
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("v", "v+(x-x1)/dt") # correct velocities for position constraints (RATTLE)
integrator.addConstrainVelocities() # constrain velocities (RATTLE)

# Accept/reject step
integrator.addComputeSum("ke", "0.5*m*v*v") # compute kinetic energy
integrator.addComputeGlobal("accept", "step(exp(-(Enew-Eold)/kT) - uniform); Enew = ke + energy")
integrator.addComputePerDof("x", "x*accept + xold*(1-accept)") # keep old or new configuration

# Run dynamics
testsystem = testsystems.AlanineDipeptideImplicit()
context = openmm.Context(testsystem.system, integrator)
context.setPositions(testsystem.positions)
print "initial potential: " + str(context.getState(getEnergy=True).getPotentialEnergy())
integrator.step(10)
print "final potential:   " + str(context.getState(getEnergy=True).getPotentialEnergy())
