# EEG signal cleaning

![MATLAB](https://img.shields.io/badge/MATLAB-R2021b-blue.svg)
[![GitHub Issues](https://img.shields.io/github/issues/mateuszschab/EEG-signal-cleaning.svg)](https://github.com/mateuszschab/EEG-signal-cleaning/issues)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)

## Basic Overview

The data came from an experiment in which the user was asked to blink the indicated eye a 
a specified number of times. The user was given 3 seconds for each blink. For further analysis 
data from epoch number 5 was selected and the first two seconds were removed. The EEG signal was recorded 
was using 19 electrodes and a sampling frequency of 500 Hz.

EEG signal from 19 electrodes:
![Map1](https://github.com/mateuszschab/EEG-blink-recognition/blob/main/img_project/EEG_signal.PNG)

In the first step, the signal was processed with a **Butterworth filter**. It allows 
to remove signal components of a certain frequency. Here a 4th order filter was used, 
high-pass with a cutoff frequency of 5 Hz. Delta waves are waves with a frequency of 1 to 3 Hz, which 
are generated by the brain during sleep. In the study, users were asked to perform specific tasks 
thus, information on low-frequency components are redundant and contain noise, so they 
they were filtered out. The Butterworth filter is characterized by flat amplitude characteristics in the 
frequency response. The row of the filter, on the other hand, determines the slope of the characteristic in the stop band.
The filter was used for each channel separately.

Another method used is the calculation of the **percentile** for the entire signal. It determines the value of the 
of the signal that reaches a set percentage level relative to the whole signal. The entire signal 
was divided into one-second windows. The windows were analyzed for the presence of 
values equal to or greater than a percentile. If such a value is present in a given window, it is 
equivalent to the occurrence of a signal peak. Peaks with high probability, are an 
undesirable artifact (e.g., blinking). In detecting artifacts, only the channel was used 
Fp1, because in it blinks have the highest reflection. However, the time window was 
removed from all 19 hundred channels.

**ICA** consists of transforming the recorded set of signals to independent components. 
For the matrix of components, the power and standard deviation of these values were calculated, so that in a 
further step remove artifacts automatically (based on a certain value) instead of 
manually selecting the data.

The last method of signal cleaning is **spatial filtering** . In this method, multiple channels are taken into 
consideration, compared to previous methods, multiple channels are taken into account simultaneously. CAR filtering 
(common average reference) involves averaging samples from all channels and using this 
average as a reference for each channel. The disadvantage of this method is that it assumes an identical level of 
noise in all channels. However, each channel may have different 
noise characteristics and dominant signal amplitude values. 

**Acknowledgements**
---

+ [West Pomeranian University of Technology in Szczecin](https://www.wi.zut.edu.pl/en/) Faculty of Computer Science and Information Technology - place of data collection and processing
