% This file is part of the package hdsort that provides a spike sorting
% algorithms. 
%
% Author: Felix Franke
% Affiliation:
%   Bernstein Center for Computational Neuroscience (BCCN) Berlin
%   and
%   Neural Information Processing Group
%   School for Electrical Engineering and Computer Science
%   Berlin Institute of Technology
%   FR 2-1, Franklinstrasse 28/29, 10587 Berlin, Germany
%   Tel: +49-30-314 26756
%
% Date: 8/24/2010
% Copyright ? Felix Franke & Technische Universit?t Berlin 2010
% Acknowledgement: This work was supported by grant DFG GRK 1589/1.
%__________________________________________________________________________
% This verion is for Robin Ince, University of Manchester, please treat
% confidential since it may contain unpublished algorithms. Furthermore,
% treat with care, it is not fully documented or tested and not in a
% state where it should be released.
%__________________________________________________________________________
%
% This program is distributed without any warranty,
% without even the implied warranty of fitness for a particular purpose.
%@_________________________________________________________________________
epochs = [ 1 10
          101 200
          201 220
          300 401];
      
hdsort.epoch.removeShort(epochs, 20)


hdsort.epoch.removeShort(epochs, 99)

hdsort.epoch.removeShort(epochs, 101)
