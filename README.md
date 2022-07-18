[![arXiv](https://img.shields.io/badge/arXiv-2207.01651-B31B1B.svg)](http://arxiv.org/abs/2207.01651) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Language](https://img.shields.io/badge/language-C++-green.svg)](https://www.cplusplus.com)



# AMIQS

<img src="logo_transparent.png" width="400"> 

##### *`Averaging Method to Integrate Quantum Systems`*<br/> *by [Pilar Hernández](https://inspirehep.net/authors/1006155?ui-citation-summary=true), [Jacobo López-Pávon](https://inspirehep.net/authors/1050355?ui-citation-summary=true), [Nuria Rius](https://inspirehep.net/authors/991635?ui-citation-summary=true) and [Stefan Sandner](https://inspirehep.net/authors/1741540?ui-citation-summary=true)*

Calculates the baryon asymmetry generated via quantum oscillations of right handed neutrinos in the early Universe. 


## Requirements

| Package | Description |
| ------ | ----------- |
| [CMake](https://cmake.org)   | generate the Makefile |
| [GSL](https://www.gnu.org/software/gsl/) | numerical library for C++ |


## Installation

The Makefile can be generate within the AMIQS folder by executing the following:

    mkdir build
    cd build
    cmake ..
    make

This should create the executable amiqs.exe!

## Testing
To check if everything is correctly installed please run

    ./amiqs.exe test

This routine may take roughly 30 seconds.

## Exemplary Usage

    ./amiqs.exe ../example.ini

where the example.ini file contains all initial conditions and its structure should not be changed.
For more explanation please see the file 'example.ini'. 



## Citation

If you use the code, please link this repository, and cite [arXiv:2207.01651](http://arxiv.org/abs/2207.01651).

## Contact

For any type of comments, questions etc. feel free to contact us at <stefan.sandner@ific.uv.es> :otter:.

## Acknowledgement 

The code makes use of the ini file decoder [inifile-cpp](https://github.com/Rookfighter/inifile-cpp) and the numerical integrator LSODA 
in the C++11 version provided by [libsoda-cxx](https://github.com/dilawar/libsoda-cxx).



