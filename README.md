# Determination of Diacetyl in E-cigarette liquid

This repository contains the scripts used to generate the standard curve to determine 2,3-butanedione (diacetyl) content in vap-juice/e-cigarette liquid. This project is a part of Stanford University's CHEM 134 class. 

### Scripts
`calculator.py` : script that calculates the linear regression, limit of detection and limit of quantification for the data. It also plots the data and linear regression
`LinReg.py` : defines the custom LinReg class, which contains methods for calculation of the linear regression, limit of detection and limit of quantification.

### Data
Raw data is located in the `\raw_chromatographs` folder. Files are in .mnova format, which can be opened and modified using MestraNova. Integrated peak areas for the raw data can be found in `peak_area.pdf`.

Contributors: Nicolas Quach, Colette Brannan