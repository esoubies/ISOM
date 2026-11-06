# ISOM

This repository contains the code to reproduce the Fourier analysis and the reconstruction
algorithms presented in the paper:

<a href="https://hal.science/hal-04969247" target="_blank">Fourier Analysis of Interference Scanning Optical Probe Microscopy</a>,  <br />
IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025. <br />
E. Soubies and W. Bacsa.

## Organization
- `Model_Analysis.m`: reproduces the illustrations of the Fourier analysis.
- `Demo.m`: provides an example of use of the two proposed reconstruction approaches.
- `Exp_Resolution/`: contains the scripts to reproduce the numerical evaluation of the resolution
  limit (Figures 5 and 6 of the paper).
- `Exp_ModelMismatches/`: contains the scripts to reproduce the robustness study with respect to
  model mismatches (Figure 8 of the paper).

## Getting started 

Download or clone this repository and run the script **Install.m** which will download the required <a href="https://biomedical-imaging-group.github.io/GlobalBioIm/">GlobalBioIm library</a> v1.2 (or more recent releases) and make neceesary changes in your Matlab path.

## Conditions of use

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/> .

Whenever you present or publish results that are based on this repository, please cite:

E. Soubies, and W. Bacsa. Fourier Analysis of Interference Scanning Optical Probe Microscopy.  IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.

