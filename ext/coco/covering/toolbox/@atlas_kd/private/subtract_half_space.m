function [P, rem] = subtract_half_space(P, index, nrm, on)
%SUBTRACT_HALF_SPACE   Removes the halfspace x.nrm-on>0 from a polyhedron

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

% Mark vertices as to side of the half space
t = max([abs(nrm); abs(on)]);
d = zeros(1, P.n);
for i=1:P.n
  d(i) = (P.v{i}'*nrm - on)/t;
%   if abs(d(i))<eps % Harry: when does this happen?
%     P.nIndices(i) = P.nIndices(i) + 1;
%     P.indices{i}  = [P.indices{i} index];
%   end
end

% Insert vertex where edge crosses halfspace
nold = P.n;
for i=1:nold
  for j=i+1:nold
    if (d(i)<0 && d(j)>0) || (d(i)>0 && d(j)<0)
      inter = intersectSets(P.indices{i}, P.indices{j}); % find common faces
      if numel(inter) >= P.k-1
        newV = P.n + 1;
        t    = -d(i)/(d(j) - d(i));
        
        P.v{newV}        = (1-t)*P.v{i} + t*P.v{j};
        P.nIndices(newV) = numel(inter) + 1;
        P.mark(newV)     = 0;
        P.indices{newV}  = [inter index];
        d(newV)          = -100.;
        P.n              = P.n + 1;
      end
    end
  end
end

% Remove wrong sided vertices
idx = find(d>0);
P.v(idx)        = [];
P.nIndices(idx) = [];
P.indices(idx)  = [];
P.mark(idx)     = [];
P.n             = P.n - numel(idx);

% Add new face to list.
P.nFaces          = P.nFaces+1;
P.face (P.nFaces) = index;
P.faceN{P.nFaces} = nrm;
P.faceO(P.nFaces) = on;

% Remove redundant faces.
P.nFaceV = zeros(1,P.nFaces);
P.faceV  = cell(1,P.nFaces);
for i=1:P.n
  for j=1:P.nIndices(i)
    idx = find(P.face==P.indices{i}(j));
    if numel(idx)>0
      P.nFaceV(idx) = P.nFaceV(idx) + 1;
      P.faceV{idx}  = [P.faceV{idx} i];
    end
  end
end

rem = find(P.nFaceV<P.k);
P.face(rem)   = [];
P.nFaceV(rem) = [];
P.faceN(rem)  = [];
P.faceO(rem)  = [];
P.faceV(rem)  = [];
P.nFaces      = P.nFaces - numel(rem);

end

function C = intersectSets(A, B)

%  Returns the intersection of two SORTED index sets

nA = numel(A);
nB = numel(B);
C  = [];

iA = 1;
iB = 1;
while iA<=nA && iB<=nB
  if A(iA)==B(iB)
    C  = [C A(iA)]; %#ok<AGROW>
    iA = iA+1;
    iB = iB+1;
  elseif A(iA)<B(iB)
    iA = iA+1;
  else
    iB = iB+1;
  end
end

end
