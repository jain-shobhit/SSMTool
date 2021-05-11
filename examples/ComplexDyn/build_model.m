function [A,B,F,V] = build_model(c,k,kappa,eqform)


As = [0 0 1 0;0 0 0 1;-2*k k -2*c c;k -2*k c -2*c];
[V,Lambda] = eig(As);
Lambda = diag(Lambda);
[~,I] = sortrows([real(Lambda), imag(Lambda)],[1 2],{'descend' 'descend'});

V = V(:,I); % arrange the order based on the real parts and then imag parts

F2 = sptensor([4 4 4]);
F3 = sptensor([4,4,4,4]);
switch eqform
    case 'non-diagonal'
        A = As*V;
        B = V;
        for i=1:4
            for j=1:4
                for k=1:4
                    F3(3,i,j,k)=-kappa*V(1,i)*V(1,j)*V(1,k);
                end
            end
        end
    case 'diagonal'
        Vinv = inv(V);
        A = Vinv*As*V;
        B = eye(4);
        for i=1:4
            for j=1:4
                for k=1:4
                    F3([(1:4)',ones(4,1)*i,ones(4,1)*j,ones(4,1)*k])=...
                        -Vinv(:,3)*kappa*V(1,i)*V(1,j)*V(1,k);
                end
            end
        end        
        
    otherwise
        error('please select from {diagonal, non-diagonal}');
end

F = {sptensor(A),F2,F3};
end