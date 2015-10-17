#!/usr/bin/env python

"""
Andersen thermostat wiht per-particle collisions

WARNING: Use cution when using this integrator with systems with constraints.

"""

from simtk import openmm, unit
from openmmtools import testsystems

timestep = 1.0 * unit.femtoseconds
temperature = 300.0 * unit.kelvin
collision_rate = 5.0/unit.picosecond
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature

# Integrator initialization.
integrator = openmm.CustomIntegrator(timestep)
integrator.addGlobalVariable("kT", kT)  # thermal energy
integrator.addGlobalVariable("p_collision", timestep * collision_rate)  # per-particle collision probability per timestep
integrator.addPerDofVariable("sigma_v", 0)  # velocity distribution stddev for Maxwell-Boltzmann (computed later)
integrator.addPerDofVariable("collision", 0)  # 1 if collision has occured this timestep, 0 otherwise
integrator.addPerDofVariable("x1", 0)  # for constraints

# Andersen thermostat stochastic velocity update
integrator.addComputePerDof("sigma_v", "sqrt(kT/m)")
integrator.addComputePerDof("collision", "step(p_collision-uniform)")  # if collision has occured this timestep, 0 otherwise
integrator.addComputePerDof("v", "(1-collision)*v + collision*sigma_v*gaussian")  # randomize velocities of collided particles

# Velocity Verlet
integrator.addUpdateContextState() # allow other forces (e.g. MonteCarloBarostat) to update system
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("x", "x+dt*v") # position full kick
integrator.addComputePerDof("x1", "x") # store positions before position constraints are applied
integrator.addConstrainPositions() # apply position constraints (CCMA)
integrator.addComputePerDof("v", "v+0.5*dt*f/m") # velocity half kick
integrator.addComputePerDof("v", "v+(x-x1)/dt") # correct velocities for position constraints (RATTLE)
integrator.addConstrainVelocities() # constrain velocities (RATTLE)

# Run dynamics
testsystem = testsystems.AlanineDipeptideImplicit()
context = openmm.Context(testsystem.system, integrator)
context.setPositions(testsystem.positions)
print "initial potential: " + str(context.getState(getEnergy=True).getPotentialEnergy())
integrator.step(10)
print "final potential:   " + str(context.getState(getEnergy=True).getPotentialEnergy())
