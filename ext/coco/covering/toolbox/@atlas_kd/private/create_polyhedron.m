function chart = create_polyhedron(chart, marg)
% CREATE_POLYHEDRON   Initialize a k-dimensional polyhedron

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

chart.P     = create_cube(chart.dim, marg*chart.R); % Initial dim-dimensional cube;
chart.s     = cell2mat(chart.P.v)/(marg*chart.R*sqrt(chart.dim)); % Vertex directions
chart.bv    = 1:chart.P.n; % Indices of available vertices
chart.nrm   = cell2mat(chart.P.faceN); % stores normal vectors associated with each face of polyhedron
chart.onrm  = chart.P.faceO; % stores distance from center to each face of polyhedron
chart.neigh = zeros(1,chart.P.nFaces);

end
