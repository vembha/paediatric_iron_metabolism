# **Code repository for the publication:**
# Systems physiological modelling of paediatric iron metabolism: application to formula-based supplementation 
### URL: [bioRxiv version](https://www.biorxiv.org/content/10.1101/2023.10.13.562306v1) need to edit

#### Numerical packages and languages used: Julia v1.7.3
<br/>
<br/>

This repository contains four files: a Julia file and three CSV files. The Julia file has two functions:
1. `model_iron_homeostasis!` encodes the ODE-based model in the manuscript.
2. `Hb_estimator!` estimates the haemoglobin level under normal, non-anaemic conditions, using the levels of red blood cells.

The CSV files contain the parameter values used for infants in different age categories, named accordingly.
