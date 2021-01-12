function P = create_cube(dim, R)
%CREATE_CUBE   Construct a hypercube
%
% P = CREATE_CUBE(DIM, R)
%
% P   - Polyhedron (struct).
% DIM - Manifold dimension (integer).
% R   - Radius (real).
%
% Construct a polyhedral object representing a hypercube.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

P.k = dim;                       % Manifold dimension
P.n = 2^dim;                     % Number of vertices
P.v = {};                        % Cell array of vertices

P.nIndices = dim*ones(1,P.n);    % Array of index cardinalities
P.indices  = cell(1,P.n);        % Cell array of vertex indices

P.mark     = zeros(1,P.n);       % Array of boolean markers
c          = -R*ones(dim,1);     % Lower left corner

% Populate vertices and vertex indices
for i=1:P.n
  P.v{i} = c;
  
  for j=1:dim
    if c(j)<0
      P.indices{i}(j) = 2*j-1;
    else
      P.indices{i}(j) = 2*j;
    end
  end
  
  carry = true;
  j = 1;
  while carry && j<=dim
    if c(j)<0
      c(j)  = R;
      carry = false;
    else
      c(j)  = -R;
      j     = j+1;
    end
  end
end

P.nFaces = 2*dim;             % Number of faces
P.face   = 1:P.nFaces;        % Array of face indices

P.nFaceV = zeros(1,P.nFaces); % Array of vertex cardinalities per face
P.faceN  = cell(1,P.nFaces);  % Cell array of face normals
P.faceV  = cell(1,P.nFaces);  % Cell array of vertex index arrays for each face
P.faceO  = zeros(1,P.nFaces); % Cell array of distances to face centers

% Populate faces and face indices
for i=1:P.nFaces
  P.faceN{i} = zeros(dim,1);
  if mod(i,2)==0
    P.faceN{i}(i/2) = -1;
  else
    P.faceN{i}((i+1)/2) = 1;
  end
  P.faceO(i) = R;
end

for i=1:P.n
  for j=1:P.nIndices(i)
    P.nFaceV(P.indices{i}(j)) = P.nFaceV(P.indices{i}(j))+1;
    P.faceV{P.indices{i}(j)}  = [P.faceV{P.indices{i}(j)} i];
  end
end

P.R = R; % Distance from center to each face

end
