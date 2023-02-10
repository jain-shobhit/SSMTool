function [B,A,Fnl,Fext] = build_model(om1,om2,om3,zeta1,zeta2,zeta3,f1,type)

n   = 3;
mass = eye(n,n);
damp = [2*zeta1*om1, 0, 0;
    0, 2*zeta2*om2, 0;
    0, 0, 2*zeta3*om3];
stiff = [om1^2 0 0;
    0 om2^2 0;
    0 0 om3^2];

switch type
    case 'cubic'
        subs2 = [1 1 1
            1 2 2
            1 3 3
            1 1 2
            1 1 3    
            2 2 2
            2 1 1
            2 3 3
            2 1 2
            2 2 3    
            3 3 3
            3 1 1
            3 2 2
            3 1 3
            3 2 3];

        vals2 = [1.5*om1^2; 
                 0.5*om1^2; 
                 0.5*om1^2; 
                 om2^2; 
                 om3^2; 
                 1.5*om2^2; 
                 0.5*om2^2; 
                 0.5*om2^2; 
                 om1^2; 
                 om3^2; 
                 1.5*om3^2;
                 0.5*om3^2; 
                 0.5*om3^2;
                 om1^2; 
                 om2^2];

        F2 = sptensor(subs2, -vals2, [2*n+1,2*n+1,2*n+1]);

        subs3 = [1 1 1 1
                1 1 2 2
                1 1 3 3
                2 2 1 1
                2 2 2 2
                2 2 3 3
                3 3 3 3
                3 3 1 1 
                3 3 2 2
                1 1 1 7
                2 2 2 7
                7 1 1 1
                7 2 2 2];

        vals3 = [0.5*(om1^2+om2^2+om3^2)*[1 1 1 1 1 1 1 1 1]';-3;-3;1;1];

        F3  = sptensor(subs3, -vals3, [2*n+1,2*n+1,2*n+1,2*n+1]);        
        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        B = [B zeros(2*n,1); zeros(1,2*n+1)];
        A = [A [0;0;-1;0;0;0]; [0 0 1 0 0 0 0]];
        Fext = [f1;0;0;0;0;0;0];
        
    case 'sphere'
        subs2 = [1 1 1
            1 2 2
            1 3 3
            1 1 2
            1 1 3    
            2 2 2
            2 1 1
            2 3 3
            2 1 2
            2 2 3    
            3 3 3
            3 1 1
            3 2 2
            3 1 3
            3 2 3
            1 1 7
            2 2 7
            3 3 7
            7 1 1
            7 2 2
            7 3 3];

        vals2 = [1.5*om1^2; 
                 0.5*om1^2; 
                 0.5*om1^2; 
                 om2^2; 
                 om3^2; 
                 1.5*om2^2; 
                 0.5*om2^2; 
                 0.5*om2^2; 
                 om1^2; 
                 om3^2; 
                 1.5*om3^2;
                 0.5*om3^2; 
                 0.5*om3^2;
                 om1^2; 
                 om2^2;
                 2;
                 2;
                 2;
                 -1;
                 -1;
                 -1];

        F2 = sptensor(subs2, -vals2, [2*n+1,2*n+1,2*n+1]);

        subs3 = [1 1 1 1
                1 1 2 2
                1 1 3 3
                2 2 1 1
                2 2 2 2
                2 2 3 3
                3 3 3 3
                3 3 1 1 
                3 3 2 2];

        vals3 = 0.5*(om1^2+om2^2+om3^2)*[1 1 1 1 1 1 1 1 1]';
        
        F3  = sptensor(subs3, -vals3, [2*n+1,2*n+1,2*n+1,2*n+1]);  
        
        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        B = [B zeros(2*n,1); zeros(1,2*n+1)];
        A = [A [0;0;2;0;0;0]; [0 0 -2 0 0 0 0]];     
        Fext = [f1;0;0;0;0;0;0];
        
    case 'none'
        subs2 = [1 1 1
            1 2 2
            1 3 3
            1 1 2
            1 1 3    
            2 2 2
            2 1 1
            2 3 3
            2 1 2
            2 2 3    
            3 3 3
            3 1 1
            3 2 2
            3 1 3
            3 2 3];

        vals2 = [1.5*om1^2; 
                 0.5*om1^2; 
                 0.5*om1^2; 
                 om2^2; 
                 om3^2; 
                 1.5*om2^2; 
                 0.5*om2^2; 
                 0.5*om2^2; 
                 om1^2; 
                 om3^2; 
                 1.5*om3^2;
                 0.5*om3^2; 
                 0.5*om3^2;
                 om1^2; 
                 om2^2];

        F2 = sptensor(subs2, -vals2, [2*n,2*n,2*n]);

        subs3 = [1 1 1 1
                1 1 2 2
                1 1 3 3
                2 2 1 1
                2 2 2 2
                2 2 3 3
                3 3 3 3
                3 3 1 1 
                3 3 2 2];

        vals3 = 0.5*(om1^2+om2^2+om3^2)*[1 1 1 1 1 1 1 1 1]';

        F3  = sptensor(subs3, -vals3, [2*n,2*n,2*n,2*n]);

        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        Fext = [f1;0;0;0;0;0];      
        
    otherwise
        error('type should be cubic/circle/none');
end

Fnl  = {F2,F3};

end