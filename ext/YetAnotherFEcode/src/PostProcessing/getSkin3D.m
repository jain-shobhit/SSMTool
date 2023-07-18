% getSkin3D(elements): computes the external faces of the mesh. 
% INPUTS:
%   elements: table of elements (Nelements X Ndofs)
% OUTPUTS:
%   skin: table of nodes labels of external faces. Each column contains the
%         nodes of one face.
%   allfaces: as skin, but with all the faces.
%
% Note: quadratic elements are plotted as linear ones, discarding mid-edge
% nodes to improve plot performances
%
% Supported elements: TET4, TET10, HEX8, HEX20
% Last modified: 12/09/2019, Jacopo Marconi, PoliMi

function [skin,allfaces] = getSkin3D(elements)

nnel  = size(elements,2); % number of nodes per element

% select faces order
switch nnel
    case {4,10} % tetrahedra
        faces = [1 2 4; 2 3 4; 1 3 4; 1 2 3];
    case 15 % wedges
        faces = [4 10 5 11 6 12; 1 7 2 8 3 9; 1 7 2 5 10 4; 2 8 3 6 11 5; 1 9 3 6 12 4];
    case {8,20} % hexahedra
        faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
end

% build face matrix
N = size(elements,1);
FACES = zeros(size(faces,2),N*size(faces,1));
count = 1;
for ii = 1:N
    vertici = elements(ii,:);
    for jj = 1:size(faces,1)
        FACES(:,count) = vertici(faces(jj,:))';
        count = count+1;
    end
end

% function [x1,icc] = OnlyNotRepeatedColumns(x)____________________________
x = FACES;
xs = sort(x,1);
[~,ia,ic] = unique(xs','rows','stable');
ia = ia';
ic = ic';
ic = ia(ic); % a unique index is given to each unique column
icc = hist(ic,1:length(ic));	% counts the occurences 
                                % (occ. after the first appear as zeros)
icc = icc == 1;                 % take only terms occurring ONCE
% x1 = xs(:,icc);__________________________________________________________

% remove all the repeated faces (internal), so to plot only external ones
skin = FACES(:,icc);
allfaces = FACES;


