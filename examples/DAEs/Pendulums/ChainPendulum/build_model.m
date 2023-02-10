function [A,B,Fnl,Fext] = build_model(m1,m,c1,c,k1,k,len,n)

% m1 = 1;
% m = 1;
% len = 1;
J = m*len^2/12;
g = 9.8;

ndof  = 2+3*(n-1);
nlamd = 2*n-1;
nauxi = 2*n-2;
qidx  = 1:ndof;
vidx  = ndof+1:2*ndof;
lidx  = 2*ndof+1:2*ndof+nlamd;
uidx  = 2*ndof+nlamd+1:2*ndof+nlamd+nauxi;

xidx = qidx([1 3:3:ndof-2]);
yidx = qidx([2 4:3:ndof-1]);
phi_idx  = [0 qidx(5:3:end)]; % add 0 for (d)phi1 (coding purpose)
dphi_idx = [0 vidx(5:3:end)];

nz = uidx(end);
B  = zeros(nz);
A  = zeros(nz);
B(qidx,qidx) = eye(ndof);
B(vidx,vidx) = diag([m1; m1; repmat([m; m; J],[n-1,1])]);
uidx_sin = uidx(1:2:end-1);
uidx_cos = uidx(2:2:end);
B(uidx_sin,uidx_sin) = eye(nauxi/2);
A(qidx,vidx) = eye(ndof);
% ddot{x}_1
A(vidx(1),1) = -k1; A(vidx(1),vidx(1)) = -c1; A(vidx(1),lidx(2)) = 1;
% ddot{y}_1
A(vidx(2),lidx(1)) = -1; A(vidx(2),lidx(3)) = 1;
% (ddot{x}_i,ddot{y}_i) with 2<=i<= n
for i=2:n
    temp_idx = 2+(i-2)*3;
    % x
    A(vidx(temp_idx+1),lidx(2*i-2)) = -1;
    if i<n
        A(vidx(temp_idx+1),lidx(2*i)) = 1;
    end
    % y
    A(vidx(temp_idx+2),lidx(2*i-1)) = -1;
    if i<n
        A(vidx(temp_idx+2),lidx(2*i+1)) = 1;  
    end
    % phi
    if i<n
        A(vidx(temp_idx+3),phi_idx(i))  = -2*k; 
        A(vidx(temp_idx+3),dphi_idx(i)) = -2*c;
    else
        A(vidx(temp_idx+3),phi_idx(i))  = -k; 
        A(vidx(temp_idx+3),dphi_idx(i)) = -c;
    end        
    if i>2
        A(vidx(temp_idx+3),phi_idx(i-1))  = k;
        A(vidx(temp_idx+3),dphi_idx(i-1)) = c;
    end
    if i<n
        A(vidx(temp_idx+3),phi_idx(i+1))  = k;
        A(vidx(temp_idx+3),dphi_idx(i+1)) = c;
        A(vidx(temp_idx+3),lidx(2*i))     = 0.5*len;
    end
    A(vidx(temp_idx+3),lidx(2*i-2)) = 0.5*len;    
    A(vidx(temp_idx+3),uidx(2*i-3)) = -0.5*len*(2*(n-i)+1)*m*g;
end
A(lidx(1),yidx(1)) = 1;
A(lidx(2),[xidx(2),uidx(1),xidx(1)]) = [1 -0.5*len -1];
A(lidx(3),[yidx(2),uidx(2),yidx(1)]) = [1 0.5*len -1];
for i=2:n-1
    temp = [xidx(i:i+1) uidx([2*i-3,2*i-1])];
    A(lidx(2*i),temp) = [-1 1 -0.5*len -0.5*len];
    temp = [yidx(i:i+1) uidx([2*i-2,2*i])];
    A(lidx(2*i+1),temp) = [-1 1 0.5*len 0.5*len];
end
A(uidx_sin, dphi_idx(2:end)) = eye(nauxi/2); % should be dphi_idx
A(uidx_cos, uidx_cos) = -2*eye(nauxi/2);

% nonlinear parts
% EOM for \ddot{\phi}
subs2 = [];
vals2 = [];
for i=2:n-1
    temp_idx = 2+(i-2)*3;
    eq_idx  = temp_idx+3;
    eq_idx  = vidx(eq_idx);
    subi = [eq_idx lidx(2*i) uidx(2*i-2)
        eq_idx lidx(2*i-2) uidx(2*i-2)
        eq_idx lidx(2*i-1) uidx(2*i-3)
        eq_idx lidx(2*i+1) uidx(2*i-3)];
    vali = 0.5*len*[-1 -1 -1 -1]';
    subs2 = [subs2; subi];
    vals2 = [vals2; vali];
end
eq_idx  = vidx(2+(n-1)*3);
subn = [eq_idx lidx(2*n-2) uidx(2*n-2)
    eq_idx lidx(2*n-1) uidx(2*n-3)];
valn = 0.5*len*[-1;-1];
subs2 = [subs2; subn];
vals2 = [vals2; valn];

% EOM for \dot{u}
subu_sin = [uidx_sin' dphi_idx(2:end)' uidx_cos'];
valu_sin = -ones(nauxi/2,1);
subu_cos = [uidx_cos' uidx_sin' uidx_sin'
    uidx_cos' uidx_cos' uidx_cos'];
valu_cos = ones(nauxi,1);
subu = [subu_sin; subu_cos];
valu = [valu_sin; valu_cos];
subs2 = [subs2; subu];
vals2 = [vals2; valu];

F2 = sptensor(subs2, vals2, [nz nz nz]);
Fnl = {F2};
Fext = zeros(nz,1);
Fext(vidx(1)) = 1;

end