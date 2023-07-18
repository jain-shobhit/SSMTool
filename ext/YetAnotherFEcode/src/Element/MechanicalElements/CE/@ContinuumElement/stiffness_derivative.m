function [Kd] = stiffness_derivative(self,x)
% this function computes the element stiffness derivative matrix
% with respect to the amplitude q of an imposed displacement
% field U (directional derivative), that is:
%      dK|         dK(Uq)|
% Kd = --|       = ----- |
%      dq|_{q=0}    dq   |_{q=0}
% and evaluated for q=0.
% x : element DOF values of U

displ = self.extract_element_data(x);
X = self.quadrature.X;
W = self.quadrature.W;
Kd= self.initialization.K;
C = self.initialization.C;
H = self.initialization.H;
Afun = self.initialization.Afun;	% fun: (th) vector -> matrix (A)
L = L_matrix(self);                 % Quadratic strain matrix: A = L.th, eps_quad = A*th
for ii = 1:length(W)
    Xi = X(:, ii);
    we = W(ii); % weights
    [G,detJ] = shape_function_derivatives(self, Xi);
    th  = G*displ;
    A = Afun(th);
    b1 = einsum('ijk,kl->ijl',L,G);
    b2 = einsum('ijk,il->jkl',b1,C*H*th);
    b = G'*b2;
    int_dK = G'*(H'*C*A + A'*C*H)*G + b;
    Kd = Kd + int_dK * detJ * we;
end
end
        