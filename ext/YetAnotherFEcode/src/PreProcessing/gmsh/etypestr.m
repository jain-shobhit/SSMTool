function [s,nne]=etypestr(etype)

% function [ETYPE,NNE]=ETYPESTR(E)
%
% Converts an element type id as from GMSH to an actual element type string
% NNE - returns teh number of nodes per element
% NOTE: ONLY A FEW TYPE OF ELEMENTS HAVE BEEN CONSIDERED HERE

switch etype
    case 1
        s='Line2';
        nne=2;
    case 2
        s='Tria3';
        nne=3;
    case 3
        s='Quad4';
        nne=4;
    case 4
        s='Tetra4';
        nne=4;
    case 5
        s='Hexa8';
        nne=8;
    case 6
        s='Prism6';
        nne=6;
    case 7
        s='Pyramid5';
        nne=5;
    case 8
        s='Line3';
        nne=3;
    case 9
        s='Tria6';
        nne=6;
    case 10
        s='Quad9';
        nne=9;
     otherwise
        s='Unknown';
        nne=0;
end
end
