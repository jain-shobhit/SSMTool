function mesh = coll_mesh(int, maps, tmi)

dim  = int.dim;
NCOL = int.NCOL;
pdim = maps.pdim;
NTST = maps.NTST;

bpnum  = NCOL+1;
cndim  = dim*NCOL;
xcnnum = NCOL*NTST;
xcndim = dim*NCOL*NTST;

mesh.tmi = tmi;
ka       = diff(tmi);
mesh.ka  = ka;

fka          = kron(ka', ones(cndim,1));
fdxka        = kron(ka', ones(dim*cndim,1));
fdpka        = kron(ka', ones(pdim*cndim,1));
fdxdxka      = kron(ka', ones(dim^2*cndim,1));
fdxdpka      = kron(ka', ones(dim*pdim*cndim,1));
fdpdpka      = kron(ka', ones(pdim^2*cndim,1));
mesh.fka     = reshape(fka, [dim xcnnum]);
mesh.fdxka   = reshape(fdxka, [dim dim xcnnum]);
mesh.fdpka   = reshape(fdpka, [dim pdim xcnnum]);
mesh.fdxdxka = reshape(fdxdxka, [dim dim dim xcnnum]);
mesh.fdxdpka = reshape(fdxdpka, [dim dim pdim xcnnum]);
mesh.fdpdpka = reshape(fdpdpka, [dim pdim pdim xcnnum]);

mesh.gka     = mesh.fka(1,:);
mesh.gdxka   = mesh.fdxka(1,:,:);
mesh.gdpka   = mesh.fdpka(1,:,:);
mesh.gdxdxka = mesh.fdxdxka(1,:,:,:);
mesh.gdxdpka = mesh.fdxdpka(1,:,:,:);
mesh.gdpdpka = mesh.fdpdpka(1,:,:,:);

mesh.fwt  = repmat(int.wt, [dim NTST]);
mesh.gwt  = mesh.fwt(1,:);
mesh.wts2 = spdiags(mesh.fwt(:), 0, xcndim, xcndim);
mesh.kas2 = spdiags(mesh.fka(:), 0, xcndim, xcndim);

t  = repmat(tmi(1:end-1)/NTST, [bpnum 1]);
tt = repmat((0.5/NTST)*(int.tm+1), [1 NTST]);
tt = t+repmat(ka, [bpnum 1]).*tt;
tnrm = tt(end);
mesh.tbp = tt(:)/tnrm;

t  = repmat(tmi(1:end-1)/NTST, [NCOL 1]);
tt = repmat((0.5/NTST)*(int.tc+1), [1 NTST]);
tt = t+repmat(ka, [NCOL 1]).*tt;
mesh.tcn = tt(:)/tnrm;

end
