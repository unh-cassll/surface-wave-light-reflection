# surface-wave-light-reflection
Code to simulate ocean surface wave fields and to model the light reflected from the synthetic air-sea interface.

## Overview

Much of this exists---in more complete and well-featured forms---elsewhere, but I wanted to try my own MATLAB implementation.

### How to use

Call _produce_simulated_sea_surface_modeled_reflection.m_ with wind speed and spatial/temporal parameters (grid size, spatial resolution, temporal sampling interval, number of frames desired). Then make yourself a snack and wait

### Remaining updates

This is a work in progress.

It does not account for volumetric effects in the water column.

Nor does it account for the polarization of sky-leaving radiance.

The code has not been optimized to run efficiently via GPU acceleration.

## Part 1: simulation of a surface wave field via random phase approach

Tessendorf, J. (2001)
"Simulating ocean water"
Simulating nature: realistic and interactive techniques.
SIGGRAPH, 1(2), 5.

## Part 2: modeling of polarized light reflected off that surface

Mobley, C. D. (2015)
"Polarized reflectance and transmittance properties of windblown sea surfaces"
Applied optics, 54(15), 4828-4849.