function [A,B,F] = build_model()

A = [-1 0;0 -sqrt(24)];
B = eye(2);

F2 = sptensor([2 2 2]);
F3 = sptensor([2,2,2,2]);
F4 = sptensor([2,2,2,2,2]);
F5 = sptensor([2,2,2,2,2,2]);

F2(2,1,1) = 1;
F3(2,1,1,1) = 1;
F4(2,1,1,1,1) = 1;
F5(2,1,1,1,1,1) = 1;

F = {F2,F3,F4,F5};
end