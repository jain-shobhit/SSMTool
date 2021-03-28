function [mass,damp,stiff,fnl,fext] = build_model(c1,c2,c3,K,epsilon)

n = 3;
mass = eye(n,n);
damp = [c1 0 0;
     0 c2 0;
     0 0 c3];
damp = damp*epsilon;
stiff = eye(n,n);
subs1 = [1 1 1
    1 1 2
    1 2 2
    2 2 2];
subs2 = subs1+1;
vals3 = [1 -3 3 -1]';
subs3 = [ones(4,1), subs1;
        2*ones(4,1), subs1;
        2*ones(4,1), subs2;
        3*ones(4,1), subs2];
vals3 = [vals3; -vals3; vals3; -vals3]*K;
f3 = sptensor(subs3, vals3, [n,n,n,n]);
f3 = f3*epsilon;
fnl = {[],f3};
fext = [1;0;0];

end