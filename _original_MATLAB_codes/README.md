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

Following the approach laid out by

Tessendorf, J. (2001)
"Simulating ocean water"
Simulating nature: realistic and interactive techniques.
SIGGRAPH, 1(2), 5.

I used the wavenumber-directional spectrum of

Elfouhaily, T., Chapron, B., Katsaros, K., & Vandemark, D. (1997)
"A unified directional spectrum for long and short wind‐driven waves"
Journal of Geophysical Research: Oceans, 102(C7), 15781-15796.

## Part 2: modeling of polarized light reflected off that surface

I followed the general framework of

Mobley, C. D. (2015)
"Polarized reflectance and transmittance properties of windblown sea surfaces"
Applied optics, 54(15), 4828-4849.

But did not flesh out the sky-leaving radiance bits.

In order to get a better handle on the relationship between the surface normal vector **n** and the vectors **h** and **s** computed through Mobley's framework, I wrote a Live Script (*surface_reflection_demonstration.mlx*) that lets you tinker with surface slope and camera incidence/azimuth angle and examine the resulting vector orientations.