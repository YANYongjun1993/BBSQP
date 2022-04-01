function [size_C]=generate_functions(size_z)
%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization
%  problems to evaluate values of function and costraints values 
%  and to compute derivatives.
%
% Last update: Feb/26/2021 by Yongjun Yan
%-------------------------------------------------------------------------
% Reorginazed the y and x variables as z
% [z1,z2,z3,z4,z5,z6] = [y1,y2,y3,x1,x2,x6]
    syms z
    for n = 1:1:size_z
        syms (sprintf('z%01d',n));
        z(n) = sym(sprintf('z%01d',n));
    end
%-------------------------------------------------------------------------
% Users should provide following 5 functions in the symbolic form.
% Phi:    terminal cost function
% L:      incremental cost function
% f:      equations of system dynamics 
% C:      input and state mixed constraints
% Cbar:   state only constraints
%--------------------------------------------------------------------------
% Example:
% if your cost function 
%  J=0.5*x(N)'*S*x(N) + 0.5*\Sum_(k=0)^(k=N-1)(x(k)'*Q*x(k) + u(k)'*R*u(k),
% and your eqautaion of system dynamics x(1,k+1) = x(1,k)+Ts*u(1,k), 
%                                       x(2,k+1) = x(2,k)+Ts*x(1,k)*u(2,k),
%                       and constraints x(1,k) + u(1,k) <= -1,       
%                                       x(2,k) + u(2,k) >= 5,
%                                       -2 <= x(1,k) <= 2, 
%                                       x(1,k)^2 + x(2,k)^2 <= 4.
%
% Then, first define all parameters and define the functions based on
% the rules of Symbolic Math Toolbox
%  Phi, L, f, C, and Cbar will be
%
% % define parameters
% s = 1e4;  q = 1,  r= 1e3; 
%
% % Terminal cost
% Phi = 0.5*s*(x1^2+x2^2);
%
% % Incremental cost
% L = 0.5*q*(x1^2+x2^2) + 0.5*r*(u1^2+u2^2);
%
% % Equations of system dynamics
% f = [x1+Ts*u1;
%      x2+Ts*x1*u2];
% 
% % Input and state mixed constriants
% C = [ x1+u1+1;
%       -x2-u2+5];
%        
% % State only constraints 
% Cbar = [ -x1-2;
%           x1-2;
%           x1^2+x2^2-4];
%--------------------------------------------------------------------------     
 
%% ------------------------------------------------------------------------
%  Users' functions: Please provide your functions
%  Phi, L, f, C, Cbar
% -------------------------------------------------------------------------
% Terminal cost


% Incremental cost
L = 5*z1 + 6*z2 + 8*z3 + 10*z4 - 7*z6 - 18*log(z5+1) - 19.2*log(z4-z5+1) + 10;

% Equations of system dynamics


% Input and state mixed constriants
C =  [-(0.8*log(z5+1) + 0.96*log(z4-z5+1) - 0.8*z6);
      z5 - z4;
      z5 - 2*z1;
      z4 - z5 - 2*z2;
      -(log(z5+1) + 1.2*log(z4-z5+1) - z6 - 2*z3 + 2);
      z1 + z2 - 1;
      -z4;
      -z5;
      -z6;
      z4 - 2;
      z5 - 2;
      z6 - 1;
      -z1;
      -z2;
      -z3;
      z1-1;
      z2-1;
      z3-1];

% State only constraints 
 
   
%% ------------------------------------------------------------------------
% Do not modify below this line 
%-------------------------------------------------------------------------

if isempty(L)
    L = zeros(0,1);
    matlabFunction(L, 'file', 'valL','vars',[z]); 
else
    matlabFunction(L, 'file', 'valL','vars',{z}); 
end

if isempty(C)
    C = zeros(0,1);
    matlabFunction(C, 'file', 'valC','vars',[z]);   
else
    matlabFunction(C, 'file', 'valC','vars',{z});
end

%%

if isempty(L)
    Lz = zeros(0,1);
    Lzz = zeros(0,1);
    
    matlabFunction(Lz, 'file', 'DLDz','vars',[z]);
    matlabFunction(Lzz, 'file', 'D2LDzz','vars',[z]);
else
    for j=1:size_z
        Lz(1,j) = diff(L,(sprintf('z%1d',j)));
    end
    matlabFunction(Lz, 'file', 'DLDz','vars',{z});   

    for j = 1:size_z
        for k = 1:size_z
            Lzz(j,k) = diff(Lz(1,k),(sprintf('z%1d',j)));
        end
    end
    matlabFunction(Lzz, 'file', 'D2LDzz','vars',{z}); 

end

%%
size_C = size(C,1);
save size_C;
if size_C == 0
    Cz = zeros(0,1);
    Czz = zeros(0,1);
    
    matlabFunction(Cz, 'file', 'DCDz','vars',[z]);
    matlabFunction(Czz, 'file', 'D2CDzz','vars',[z]);
else
    for j=1:size_C
        for k = 1:size_z
            Cz(j,k) = diff(C(j),(sprintf('z%1d',k)));
        end
    end
    matlabFunction(Cz, 'file', 'DCDz','vars',{z});   

    for n=1:size_C
        for j = 1:size_z
            for k = 1:size_z
                Czz(j,k,n) = diff(Cz(n,k),(sprintf('z%1d',j)));
            end
        end
    end
    matlabFunction(Czz, 'file', 'D2CDzz','vars',{z}); 
end

return

