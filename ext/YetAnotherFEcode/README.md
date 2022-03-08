# YetAnotherFEcode
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4011281.svg)](https://doi.org/10.5281/zenodo.4011281)

A simple MATLAB-based code for implementing the Finite Element method in an object oriented fashion.

The main idea behind this package is to enable rapid prototyping and reproducible research related to finite element applications 
and/or reduced-order modeling in a user-friendly MATLAB environment. 

On one hand, commercial packages lack the flexibility needed for testing new ideas essential for research, especially in the context of 
reduced-order modeling, where FE problems are indeed just applications but still require mild intrusion/access to the functionality. 
Open source packages, on the other hand, allow endless access to the implementation but tend to be very cumbersome to hack, and 
require significant time and training to be able to test even the simplest of ideas. This code is particularly aimed towards 
users/researchers who are interested in intrusive finite-element modeling without getting lost in gory details of open source FE packages.  

A distinguishing aspect of this package is that apart from using exisiting elements in our library, one can program new elements with
relative ease and flexibility. These elements may also arise from multi-physics problems, e.g., thermo-mechanical 
problems which involve the numerical solution of different partial differential equations governing heat and momentum balance 
on the same physical domain. 

Without worrying about the cumbersome details of finite-element assembly, a researcher can simply focus on the 
element-level implementation to quickly obtain results. At the same time, developers are also encouraged to contribute new and 
alternative ideas to improve this environment and potentially publish them, allowing future users to access and build upon their work. 
This allows for rapid developement and testing of ideas, especially valuable in research efforts.

To use the code, simply add the main folder and its contents to the MATLAB path. Feel free to play with examples in the examples directory.
Further usage and development instructions to follow.  

To showcase the relevance, please cite the following reference if you use this package in your work

Shobhit Jain, Jacopo Marconi & Paolo Tiso (2020). YetAnotherFEcode. Zenodo. http://doi.org/10.5281/zenodo.4011281

Please report any issues/bugs to Shobhit Jain <shjain@ethz.ch>
