# Local extreme value detection for sinusoidal signals corrupted by noise

**A numerical study on (damped) sinusoidal signals**


## Abstract

The analysis of ultrasonic signals is characterised by the analysis of signal data in time domain. In the ultrasonic pulse transmission method, the signal analysis is often aimed at detecting the onset point of the incoming compression or shear wave. This makes it possible to determine the sound propagation time and, if the measuring distance is known, the speed of sound as an essential material parameter.

This starting point is not always easy to localise due to signal noise and possible interference. As an aid and to stabilise the signal analysis, it is useful to estimate the location of the first local extreme value (minimum or maximum) directly behind the onset point of the compression or shear wave. This can significantly reduce the search interval for the onset point and increase the robustness of the entire subsequent analysis.

Experience shows that the first local minimum or maximum is not always identical to the global minimum or maximum. This circumstance very easily leads to false detection. In this work, a search method is presented that enables the robust detection of the first local extreme value's location and comes with a moderate computational effort.

> [!NOTE]
> The entire content of this repository was conceived, implemented and tested by Jakob Harden using the scientific numerical programming language of GNU Octave 6.2.0 and TeXlive/TeXStudio.


## Table of contents

- [Licence](#licence)
- [Prerequisites](#prerequisites)
- [Directory and file structure](#directory-and-file-structure)
- [Installation instructions](#installation-instructions)
- [Usage instructions](#usage-instructions)
- [Help and Documentation](#help-and-documentation)
- [Related work](#related-work)
- [Revision and release history](#revision-and-release-history)


## Licence

Copyright 2025 Jakob Harden (jakob.harden@tugraz.at, Graz University of Technology, Graz, Austria)

This file is part of the PhD thesis of Jakob Harden.

All GNU Octave function files (\*.m) are licensed under the *GNU Affero General Public Licence v3.0*. See also licence information file "LICENSE".

All \*.pdf, \*.png and \*.mat files in the directory **./octave/scicol\_crameri** are licensed under the MIT licence. See also **./octave/scicol_crameri/README.md**. 

All other files are licensed under the *Creative Commons Attribution 4.0 International* licence. See also licence information: [Licence deed](https://creativecommons.org/licenses/by/4.0/deed.en)


## Prerequisites

To be able to use the scripts, GNU Octave 6.2.0 (or a higher version) need to to be installed.
GNU Octave is available via the package management system on many Linux distributions. Windows users have to download the Windows version of GNU Octave and to install the software manually.

[GNU Octave download](https://octave.org/download)


To compile the LaTeX documents in this repository a TeXlive installation is required. The TeXlive package is available via the package management system on many Linux distributions. Windows users have to download the Windows version of TeXlive and to install the software manually.

[TeXlive download](https://www.tug.org/texlive/windows.html)

> [!TIP]
> **TeXStudio** is a convinient and powerful solution to edit TeX files!


## Directory and file structure

GNU Octave script files (\*.m) are written in the scientific programming language of GNU Octave 6.2.0. LaTeX files (\*.tex) are written with compliance to TeXlive version 2020.20210202-3. Text files are generally encoded in UTF-8.

```
peakholdmax   
├── latex   
│   ├── adaptthemePresRIP.tex   
│   ├── biblio.bib   
│   └── main.tex   
├── octave   
│   ├── results   
│   │   └── test_peakholdmax   
│   │       ├── fig   
│   │       ├── oct   
│   │       ├── png   
│   │       └── noise_Nsmp8000_Nmc500.oct   
│   ├── scicol_crameri   
│   │   ├── README.md   
│   │   ├── colormap_crameri.m   
│   │   ├── *.mat   
│   │   ├── *.pdf   
│   │   └── *.png   
│   ├── tools   
│   │   ├── tool_det_localmax.m   
│   │   ├── tool_gen_noise.m   
│   │   ├── tool_gen_signal.m   
│   │   └── tool_scale_noise2snr.m   
│   ├── init.m   
│   └── test_peakholdmax.m   
├── published   
├── LICENSE   
├── README.html   
└── README.md   
```

- **peakholdmax** ... main program directory
  - *LICENSE* ... AGPLv3 licence information file
  - *README.md* ... this file, information about the program (CC BY-4.0)
  - *README.html* ... html version of this file (CC BY-4.0)
- **peakholdmax/latex** ... directory, LaTeX documents, presentation slides (CC BY-4.0)
  - *adaptthemePresRIP.tex* ... file, LaTeX beamer class configuration
  - *biblio.bib* ... file, bibtex bibliography
  - *main.tex* ... file, main LaTeX document, presentation source code
- **peakholdmax/octave** ... directory, GNU Octave script files and analysis results (AGPLv3)
  - *init.m* ... function file, initialization script, load additional packages, add subdirectories to path environment variable
  - *test\_peakholdmax.m* ... GNU Octave function file, main script, algorithm test
- **peakholdmax/octave/results/test\_peakholdmax** ... directory, analysis results (CC BY-4.0)
  - *noise\_Nsmp8000\_Nmc500.oct* ... GNU Octave binary file, standard noise data matrix (to reproduce analysis results), see also function file *tool\_gen\_noise.m*
  - **/oct** ... directory, GNU Octave binary files, analysis results
  - **/fig** ... directory, GNU Octave binary files, figures
  - **/png** ... directory, portable network graphics, figures
- **peakholdmax/octave/scicol\_crameri** ... scientific colour schemes, (AGPLv3, CC BY-4.0, MIT)
  - *README.md* ... information about the scientific colour schemes by Fabio Crameri et. al. and additional licence information. (CC BY-4.0)
  - *colormap\_crameri.m* ... GNU Octave function file, colour map loader (AGPLv3)
  - *\*.mat* ... MatLAB binary files, colour definitions (MIT)
  - *\*.png* ... colour scheme documentation (MIT)
  - *\*.pdf* ... colour scheme documentation (MIT)
- **peakholdmax/octave/tools** ... directory, GNU Octave tool scripts (AGPLv3)
  - *tool\_det\_localmax.m* ... GNU Octave function file, detect the first local maximum of a low-pass signal corrupted by noise
  - *tool\_gen\_noise.m* ... GNU Octave function file, generate standard noise data for reproducible analysis results
  - *tool\_gen\_signal.m* ... GNU Octave function file, generate (damped) sinusoidal signal
  - *tool\_scale\_noise2snr.m* ... GNU Octave function file, scale noise data w.r.t. signal power and signal-to-noise ratio
- **peakholdmax/published** ... directory, published documents, archives (AGPLv3, CC BY-4.0)


## Installation instructions

1. Copy the program directory **peakholdmax** to a location of your choice. e.g. **/home/acme/science/peakholdmax**.   
2. Open GNU Octave.   
3. Make the program directory **/home/acme/science/peakholdmax** the working directory.   


## Usage instructions

1. Open GNU Octave.   
2. Initialize program.   
3. Run script files.   


### Initialize program (command line interface)

The *init* command initializes the program. The initialization must be run once before executing all the other functions. The command is adding the subdirectories included in the main program directory to the 'path' environment variable. Furthermore, *init* is loading additional GNU Octave packages required for the program execution.

```
    octave: >> init;   
```


### Execute function file (command line interface)

```
    octave: >> test_peakholdmax('compile'); # compile test signals and the paramater variation (Monte-Carlo simulation)   
    octave: >> test_peakholdmax('plot'); # plot test signals and the paramater variation results      
```

> [!NOTE]
> To reproduce all analysis results shown in the presentation in **peakholdmax/latex**, run all above-mentioned commands.


## Help and Documentation

All function files contain an adequate function description and instructions on how to use the functions. The documentation can be displayed in the GNU Octave command line interface by entering the following command:

```
    octave: >> help function_file_name;   
```


## Related work

The local extreme value detection method presented here is used in other publications of the author of this work. More information about this publications and preliminary results are available online.

- Harden, J. (2025) "Ultrasonic Pulse Transmission Method: Enhanced threshold detection method enabling fast and robust primary wave time range estimation - Website". [Website](https://jakobharden.at/wordpress/planned-publication-paper-2/)   
- Harden, J. (2025) "Ultrasonic Pulse Transmission Method: A robust method to detect the onset point of shear waves using cross-correlation analysis - Website". [Website](https://jakobharden.at/wordpress/planned-publication-paper-3/)   
- Harden, J. (2025) "Ultrasonic pulse transmission method: Fully automated primary wave onset point detection - preliminary results". [DOI: 10.3217/z5dcs-hf845](https://doi.org/10.3217/z5dcs-hf845)   
- Harden, J. (2025) "Ultrasonic Pulse Transmission Method: A robust method to detect the onset point of shear waves using cross-correlation analysis - preliminary results". [DOI: 10.13140/RG.2.2.29301.90084](https://doi.org/10.13140/RG.2.2.29301.90084)   


## Revision and release history

### 2025-07-26, version 1.0.0

- published/released version 1.0.0, by Jakob Harden   
- URL: [Repository of Graz University of Technology](https://repository.tugraz.at/)   
- Presentation, DOI: [10.3217/yj25h-0h842](https://doi.org/10.3217/yj25h-0h842)   
- GNU Octave code, DOI: [10.3217/m8cyc-2pd96](https://doi.org/10.3217/m8cyc-2pd96)   
- LaTeX code, DOI: [10.3217/tvap0-9y787](https://doi.org/10.3217/tvap0-9y787)   

