# SAR Sub-Aperture Processor

This repository contains Python utilities for Synthetic Aperture Radar
(SAR) sub-aperture processing.\
It allows the generation of sub-aperture images from SLC products and
export of both complex and intensity data to GeoTIFF with down-sampling
options in both range and azimuth.\
It is designed to be part of a larger SAR processing workflow.

------------------------------------------------------------------------

## Features

-   Azimuth sub-aperture decomposition using a rolling window on a
    normalized Doppler frequency axis \[-0.5, 0.5\]\
-   Full control over window size and step increment, enabling
    configurable overlap between sub-apertures\
-   Flexible output formats: 8-bit, 16-bit or 32-bit GeoTIFF
    with optional georeferencing\
-   Optional GIF animation export for visualization\
-   Support for both IQ (In-phase / Quadrature) and Amplitude/Phase SLC
    formats

------------------------------------------------------------------------

## Use Cases

-   RFI mitigation\
    ![](https://github.com/mikhelif/sar-subaperture-processor/blob/main/examples/20260130_vv_02_01.gif)
    Loop thourgh the generated sub-aperture.
    Identify Doppler intervals contaminated by RFI.
    Generate a final sub-aperture using the maximum clean bandwidth.
    
-   Increase the number of training samples for ATR models\
-   Along-track interferometry\

------------------------------------------------------------------------

## Processing Overview

1.  Read focused SLC SAR data\
2.  FFT along azimuth transforming SLC data to the frequency domain\
3.  Apply a raised-cosine bandpass filter at different frequency
    positions\
4.  Inverse FFT back to spatial domain to obtain a sub-aperture stack\
5.  Compute intensity (optional)\
6.  Downsample and export GeoTIFFs

------------------------------------------------------------------------

## Sub-Aperture Generation Method

Sub-apertures are generated using a rolling window on a normalized
Doppler frequency axis \[-0.5, 0.5\].

It is controlled by 2 parameters:

-   `win_frac` \[0--1\] : Controls the window width, i.e., the fraction
    of the total azimuth bandwidth used per sub-aperture\
-   `step_frac` \[0--1\] : Step between consecutive sub-apertures as a
    fraction of the total azimuth bandwidth

These parameters also allow control over overlap between consecutive
sub-apertures:

overlap = max(0, win_frac - step_frac)

For example:

If you want to split the SLC into 2 sub-bands around the Doppler
centroid, both parameters must be set to 0.5\
(assuming the SLC was properly deramped and demodulated so that Doppler
centroid = 0).

This results in:

-   No overlap\
-   Two symmetric windows centered at ±0.25\


------------------------------------------------------------------------

## Resolution and Downsampling

Reducing azimuth bandwidth degrades azimuth resolution:

ρ_az,sub = ρ_az,full / win_frac

Azimuth downsampling can be applied proportionally to the
bandwidth reduction:

downsample factor ≈ 1 / win_frac

Optional range downsampling is also supported

------------------------------------------------------------------------

## Input Requirements

SLC data must be provided in GeoTIFF format.

Two-band structure:

**IQ format:** - Band 1 = In-phase\
- Band 2 = Quadrature

**Amplitude/Phase format:** - Band 1 = Amplitude\
- Band 2 = Phase

The processor converts both formats internally to complex data before
processing.

