function [data, y] = int(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T    = u(pr.T_idx);
x    = u(pr.xbp_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

y = (0.5*T/maps.NTST)*mesh.gwt*gcn';

end
