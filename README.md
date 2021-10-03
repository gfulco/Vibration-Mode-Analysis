# Vibration-Mode-Analysis
A simple matlab script that helps the user to analyze vibration mode up to 1000Hz. 

notchfilter.m -> This file contains the script for the set of notch filters used to remove the fundametal modes of vibrations found in the beam. 

movav.m -> It is a moving average filter(low pass filter) used to remove high frequencies mode that are not of interests.

main.m -> main file of the script.

peaks.m -> a simplified and custom version of the matlab function findpeaks() used for in the scope of this project.

hammer_data_x.txt -> set of sample data to run the script.

The presentation folder contains the LaTeX files for a pdf report of the script
