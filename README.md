# FastDecayReco

## Overview

`FastDecayReco.cc` implements **fast hadron reconstruction** in the sPHENIX software framework. The module is intended for inclusive and identified hadron analyses and provides a lightweight, analysis-oriented reconstruction chain.

The code is **based on and adapted from `KshortReconstruction.cc`**, originally developed by **Anthony Frawley**, with extensions and modifications to support generic hadron reconstruction and analysis use cases.

## Features

The `FastDecayReco` module includes:

- **Track selection**
  - Quality cuts based on tracking states and kinematic variables
- **Vertexing and DCA information**
  - Access to track–vertex quantities for analysis-level selections
- **Calorimeter matching**
  - Projection of tracks to CEMC, HCALIN, and HCALOUT
  - Matching to calorimeter clusters and towers
- **dE/dx information**
  - Extraction of TPC dE/dx for particle identification studies
- **Trigger information**
  - Storage of relevant trigger bits for event selection
- **QA histograms**
  - ROOT histograms for monitoring track quality, matching performance, and kinematic distributions
- **Output ntuples**
  - Analysis-ready TNtuples containing per-track and per-event variables

## Repository Structure

The main components of this branch are organized as follows:

```text
src/
  FastDecayReco.cc
  FastDecayReco.h

condor/
  submit scripts, job wrappers, and configuration files

analysis/
  ROOT macros and analysis scripts
