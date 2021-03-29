# SHPasteroids
Code, input and output files developed for the Senior Honours Project: Save the Earth from the Killer Asteroids!

UKST_all.py - Reads in entire UKST catalog.txt and puts out precovery plates for stated asteroids 

UKST_SC.py, POSS_II_SC.py, ESO_SC.py - Reads in SuperCOSMOS survey plates of UKST catalog.txt, POSS-II catalog.txt and ESO catalog.txt and puts out precovery plates for stated asteroids. ESO_SC.py also needs plate numbers listed in ESO survey numbers.txt

plate_download.py - reads in code output files found in output files folder to download images from the SuperCOSMOS web interface

plots.py - reads in code output files to produce histograms and scatter plots of the results

sentry.txt, ESA.txt, MPC emoid sorted.txt - lists of potentially hazardous asteroids that are read in the main codes

output files - named by plate catalogue and asteroid list, '<20' references magnitudes below 20, 'low del' contain only positional uncertainties within the plate boundaries, deltas contain the uncertainties in degrees, needed to plot the scatter plot in plots.py
