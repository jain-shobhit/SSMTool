function data = torus_bc_update(data, T0, T, x0, x1, p)

coeffs = data.F*x0;
coeffs = reshape(coeffs,[data.dim,2*data.N+1]);
f00    = sum(repmat((1:data.N),[data.dim,1]).*coeffs(:,3:2:end),2); % \sum kb_k
data.x00 = x0(1:data.dim);
data.f00 = f00;

if data.autonomous
    data.x0 = x0(1:data.dim);
    data.f0 = data.fhan(T0,data.x0,p);
end

end