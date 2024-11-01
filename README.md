# Transient Interactions
This repository currently contains the necessary scripts for measuring the distance distribution from a given set of labelling distances and then calculating the fraction of transient interactions and polymer scaling exponent using a previuosly established adjusted polymer model. 

## polymer.py
Reference file containing functions for the adjusted polyer model. Is imported as a module in the fit_pr.py file. Users should not have to edit polymer.py, simply ensure it is located in the desired working directory

## fit_pr.py
File which takes dye.py as an argument and calculates the resulting p(r) distribution, fraction of transient interactions, and polymer scaling exponent. The final two are displayed as an output and saved to .npy files.

## dye.npy
example labelling distance file. Contains the distance between two labelling positions of a model peptide for 1/100 of every frames from coarse-grained simulation. Units of distance are in nanometers.

## get_lcp.py
File which takes any sequnece file as an argument and calculates the highest and lowest values of effective patch length (l_{cp}) as well as their position along the sequence and the size of the charge patch in residues.

## seq.dat
Example sequence file containing the alphabetic sequence for the disordered region of the protein SLCA1. Can be replaced by any target protein sequence.
