function [data, y] = adjt_int(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(pr.T_idx);
x = u(pr.xbp_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;

J_xbp = (0.5*T/maps.NTST)*gdxcn(:)';
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';

y = [J_xbp J_T];

end
