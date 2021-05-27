function [mass,damp,stiff,fnl,fext] = build_model(c1,c2,b1,b2,f1,f2)

n = 2;
mass = eye(n,n);
damp = [c1, 0;
     0, c2];
stiff = [1 0;0 4];
subs2 = [1 1 2
    2 1 1];
vals2 = [b1 b2]';
F2 = sptensor(subs2, vals2, [n,n,n]);
fnl = {F2};
fext = [f1;f2];

end