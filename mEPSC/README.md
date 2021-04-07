# mEPSC Analysis

Collection of mEPSC analysis scripts, see below for a brief description

## Clipper.m

*Clips a portion off of the beginning of the recording and saves as .phy*

## conc_all.m

*Concatenates all recordings in a filepath into a single .phy for loading into StimFit*

## Count_sep.m

*Seperates the events by the wave number using cumsum from the Eventer output*

## Noise.m

*Scans through eventer output and detects the noise level of the rscomp recording*

## overlay_plotter.m

*plots neat looking overlay plots of ensemble averages and condition median*

## property_extractor.m

*Scans through the ml_out output and extracts the values other than synaptic, ie Input, Rm, Cap etc*

## rolling_plot.m

*Gives the illusion of a rolling live plot with selectable window*

## trainingset.m

*Allows the user to generate a training set for use in creating a machine learning model*

## Trimmer.m

*Script that merely opens up the file with ephysIO, trims if necessary, then resaves it. Definitely worth combining with Clipper, i'm just lazy ...*

# Suggested processing pathway

**Data cleaning**

1.  Raw mEPSC input as .tdms

2.  Trim/Clip with `trimmer.m` or `clipper.m`

3.  Compensate for series resistance change across the recording with `RS_Comp.m`

**Early properties (to be automated)**

1.  Generate a pearsons fit with Eventer

2.  Scan through this calculate the noise levels with `noise.m`

3.  Extract other relevent properties with `property_extractor.m`

**Model training**

1.  Calculate a training set with `trainingset.m`

2.  Train model with eventer and extract events using this model into a new output

**Process for Stimfit**

1.  Seperate the events using `count_sep.m`

2.  Concatenate ensemble averages into single array for StimFit with `conc_all.m`

**Stimfit**

1.  Perform stimfit analysis making sure to save average trace as .h5

2.  Plot neat overlay traces with `overlay_plotter.m`
