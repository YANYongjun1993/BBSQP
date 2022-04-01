function [size_C, size_Cbar]=generate_functions(size_x,size_u)
%%-------------------------------------------------------------------------        
% Get users' functions for generating sub-functions of users' optimization
%  problems to evaluate values of function andcostraints values 
%  and to compute derivatives.
%
% Last update: Feb/26/2021 by Yongjun Yan
    syms x u
    for n = 1:1:size_x
        syms (sprintf('x%01d',n));
        x(n) = sym(sprintf('x%01d',n));
    end
    for m = 1:1:size_u
        syms (sprintf('u%01d',m));
        u(m) = sym(sprintf('u%01d',m));
    end
%     for h = 1:1:size_x^2
%         syms (sprintf('pf%01d',h));
%         pf(h) = sym(sprintf('pf%01d',h));
%     end
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
ts = 2;
% load qf
% pf = qf;
pf = [0.7474, 0, 0, 0.0534, 0 ,0; 0, 0.7497, 0, -0, 0.0535, 0; 0, 0, 1.3037, 0, 0, 0.0772; 0.0534, -0, 0, 0.008, 0, 0; 0, 0.0535, 0, 0, 0.008, 0; 0, 0, 0.0772, 0, 0, 0.0094];
% define the model parameters
J = diag([918 920 1365]);
Jp = inv(J);

% Penalty Costs
q = 10*diag([10,10,10,1,1,1]);
r  = diag(0.1*ones(3,1));

% drive the system to the origin
xtrg = zeros(6,1);
% Terminal cost
Phi = 0.5*(   (x1-xtrg(1))*pf(1,1)*(x1-xtrg(1)) + (x1-xtrg(1))*pf(1,2)*(x2-xtrg(2)) + (x1-xtrg(1))*pf(1,3)*(x3-xtrg(3)) + (x1-xtrg(1))*pf(1,4)*(x4-xtrg(4)) + (x1-xtrg(1))*pf(1,5)*(x5-xtrg(5)) + (x1-xtrg(1))*pf(1,6)*(x6-xtrg(6)) ...
          + (x2-xtrg(2))*pf(2,1)*(x1-xtrg(1)) + (x2-xtrg(2))*pf(2,2)*(x2-xtrg(2)) + (x2-xtrg(2))*pf(2,3)*(x3-xtrg(3)) + (x2-xtrg(2))*pf(2,4)*(x4-xtrg(4)) + (x2-xtrg(2))*pf(2,5)*(x5-xtrg(5)) + (x2-xtrg(2))*pf(2,6)*(x6-xtrg(6)) ...
          + (x3-xtrg(3))*pf(3,1)*(x1-xtrg(1)) + (x3-xtrg(3))*pf(3,2)*(x2-xtrg(2)) + (x3-xtrg(3))*pf(3,3)*(x3-xtrg(3)) + (x3-xtrg(3))*pf(3,4)*(x4-xtrg(4)) + (x3-xtrg(3))*pf(3,5)*(x5-xtrg(5)) + (x3-xtrg(3))*pf(3,6)*(x6-xtrg(6)) ...
          + (x4-xtrg(4))*pf(4,1)*(x1-xtrg(1)) + (x4-xtrg(4))*pf(4,2)*(x2-xtrg(2)) + (x4-xtrg(4))*pf(4,3)*(x3-xtrg(3)) + (x4-xtrg(4))*pf(4,4)*(x4-xtrg(4)) + (x4-xtrg(4))*pf(4,5)*(x5-xtrg(5)) + (x4-xtrg(4))*pf(4,6)*(x6-xtrg(6)) ...
          + (x5-xtrg(5))*pf(5,1)*(x1-xtrg(1)) + (x5-xtrg(5))*pf(5,2)*(x2-xtrg(2)) + (x5-xtrg(5))*pf(5,3)*(x3-xtrg(3)) + (x5-xtrg(5))*pf(5,4)*(x4-xtrg(4)) + (x5-xtrg(5))*pf(5,5)*(x5-xtrg(5)) + (x5-xtrg(5))*pf(5,6)*(x6-xtrg(6)) ...
          + (x6-xtrg(6))*pf(6,1)*(x1-xtrg(1)) + (x6-xtrg(6))*pf(6,2)*(x2-xtrg(2)) + (x6-xtrg(6))*pf(6,3)*(x3-xtrg(3)) + (x6-xtrg(6))*pf(6,4)*(x4-xtrg(4)) + (x6-xtrg(6))*pf(6,5)*(x5-xtrg(5)) + (x6-xtrg(6))*pf(6,6)*(x6-xtrg(6)) );
      
% Phi = 0.5*(   (x1-xtrg(1))*pf(1)*(x1-xtrg(1)) + (x1-xtrg(1))*pf(2)*(x2-xtrg(2)) + (x1-xtrg(1))*pf(3)*(x3-xtrg(3)) + (x1-xtrg(1))*pf(4)*(x4-xtrg(4)) + (x1-xtrg(1))*pf(5)*(x5-xtrg(5)) + (x1-xtrg(1))*pf(6)*(x6-xtrg(6)) ...
%   + (x2-xtrg(2))*pf(7)*(x1-xtrg(1)) + (x2-xtrg(2))*pf(8)*(x2-xtrg(2)) + (x2-xtrg(2))*pf(9)*(x3-xtrg(3)) + (x2-xtrg(2))*pf(10)*(x4-xtrg(4)) + (x2-xtrg(2))*pf(11)*(x5-xtrg(5)) + (x2-xtrg(2))*pf(12)*(x6-xtrg(6)) ...
%   + (x3-xtrg(3))*pf(13)*(x1-xtrg(1)) + (x3-xtrg(3))*pf(14)*(x2-xtrg(2)) + (x3-xtrg(3))*pf(15)*(x3-xtrg(3)) + (x3-xtrg(3))*pf(16)*(x4-xtrg(4)) + (x3-xtrg(3))*pf(17)*(x5-xtrg(5)) + (x3-xtrg(3))*pf(18)*(x6-xtrg(6)) ...
%   + (x4-xtrg(4))*pf(19)*(x1-xtrg(1)) + (x4-xtrg(4))*pf(20)*(x2-xtrg(2)) + (x4-xtrg(4))*pf(21)*(x3-xtrg(3)) + (x4-xtrg(4))*pf(22)*(x4-xtrg(4)) + (x4-xtrg(4))*pf(23)*(x5-xtrg(5)) + (x4-xtrg(4))*pf(24)*(x6-xtrg(6)) ...
%   + (x5-xtrg(5))*pf(25)*(x1-xtrg(1)) + (x5-xtrg(5))*pf(26)*(x2-xtrg(2)) + (x5-xtrg(5))*pf(27)*(x3-xtrg(3)) + (x5-xtrg(5))*pf(28)*(x4-xtrg(4)) + (x5-xtrg(5))*pf(29)*(x5-xtrg(5)) + (x5-xtrg(5))*pf(30)*(x6-xtrg(6)) ...
%   + (x6-xtrg(6))*pf(31)*(x1-xtrg(1)) + (x6-xtrg(6))*pf(32)*(x2-xtrg(2)) + (x6-xtrg(6))*pf(33)*(x3-xtrg(3)) + (x6-xtrg(6))*pf(34)*(x4-xtrg(4)) + (x6-xtrg(6))*pf(35)*(x5-xtrg(5)) + (x6-xtrg(6))*pf(36)*(x6-xtrg(6)) );

% Incremental cost
L = 0.5*(   (x1-xtrg(1))*q(1,1)*(x1-xtrg(1)) + (x1-xtrg(1))*q(1,2)*(x2-xtrg(2)) + (x1-xtrg(1))*q(1,3)*(x3-xtrg(3)) + (x1-xtrg(1))*q(1,4)*(x4-xtrg(4)) + (x1-xtrg(1))*q(1,5)*(x5-xtrg(5)) + (x1-xtrg(1))*q(1,6)*(x6-xtrg(6)) ...
          + (x2-xtrg(2))*q(2,1)*(x1-xtrg(1)) + (x2-xtrg(2))*q(2,2)*(x2-xtrg(2)) + (x2-xtrg(2))*q(2,3)*(x3-xtrg(3)) + (x2-xtrg(2))*q(2,4)*(x4-xtrg(4)) + (x2-xtrg(2))*q(2,5)*(x5-xtrg(5)) + (x2-xtrg(2))*q(2,6)*(x6-xtrg(6)) ...
          + (x3-xtrg(3))*q(3,1)*(x1-xtrg(1)) + (x3-xtrg(3))*q(3,2)*(x2-xtrg(2)) + (x3-xtrg(3))*q(3,3)*(x3-xtrg(3)) + (x3-xtrg(3))*q(3,4)*(x4-xtrg(4)) + (x3-xtrg(3))*q(3,5)*(x5-xtrg(5)) + (x3-xtrg(3))*q(3,6)*(x6-xtrg(6)) ...
          + (x4-xtrg(4))*q(4,1)*(x1-xtrg(1)) + (x4-xtrg(4))*q(4,2)*(x2-xtrg(2)) + (x4-xtrg(4))*q(4,3)*(x3-xtrg(3)) + (x4-xtrg(4))*q(4,4)*(x4-xtrg(4)) + (x4-xtrg(4))*q(4,5)*(x5-xtrg(5)) + (x4-xtrg(4))*q(4,6)*(x6-xtrg(6)) ...
          + (x5-xtrg(5))*q(5,1)*(x1-xtrg(1)) + (x5-xtrg(5))*q(5,2)*(x2-xtrg(2)) + (x5-xtrg(5))*q(5,3)*(x3-xtrg(3)) + (x5-xtrg(5))*q(5,4)*(x4-xtrg(4)) + (x5-xtrg(5))*q(5,5)*(x5-xtrg(5)) + (x5-xtrg(5))*q(5,6)*(x6-xtrg(6)) ...
          + (x6-xtrg(6))*q(6,1)*(x1-xtrg(1)) + (x6-xtrg(6))*q(6,2)*(x2-xtrg(2)) + (x6-xtrg(6))*q(6,3)*(x3-xtrg(3)) + (x6-xtrg(6))*q(6,4)*(x4-xtrg(4)) + (x6-xtrg(6))*q(6,5)*(x5-xtrg(5)) + (x6-xtrg(6))*q(6,6)*(x6-xtrg(6)) ...
          + u1*r(1,1)*u1 + u1*r(1,2)*u2 + u1*r(1,3)*u3 ...
          + u2*r(2,1)*u1 + u2*r(2,2)*u2 + u2*r(2,3)*u3 ...
          + u3*r(3,1)*u1 + u3*r(3,2)*u2 + u3*r(3,3)*u3 );

% Equations of system dynamics
% f = [ x1 - ts*(-1/J(1,1))* ((J(3,3)-J(2,2))*x2*x3 + u1);
%       x2 - ts*(-1/J(2,2))* ((J(1,1)-J(3,3))*x1*x3 + u2);
%       x3 - ts*(-1/J(3,3))* ((J(2,2)-J(1,1))*x1*x2 + u3);
%       x4 - ts * (x1 + x2*sin(x4)*tan(x5) + x3*cos(x4)*tan(x5));
%       x5 - ts * (x2*cos(x4) - x3*sin(x4));
%       x6 - ts * (x2*sin(x4)*1/cos(x5) + x3*cos(x4)*1/cos(x5))];
 
f = [   x1 + ts * (-1/J(1,1))* ((J(3,3)-J(2,2))*x2*x3 - u1);
        x2 + ts * (-1/J(2,2))* ((J(1,1)-J(3,3))*x1*x3 - u2);
        x3 + ts * (-1/J(3,3))* ((J(2,2)-J(1,1))*x1*x2 - u3);
        x4 + ts * (x1 + x2*sind(x4)*tand(x5) + x3*cosd(x4)*tand(x5));
        x5 + ts * (x2*cosd(x4) - x3*sind(x4));
        x6 + ts * (x2*sind(x4)*1/cosd(x5) + x3*cosd(x4)*1/cosd(x5))];

% Input and state mixed constriants
C =  [u1-1;
      u2-1;
      u3-1;
      -u1-1;
      -u2-1;
      -u3-1];

% State only constraints 
Cbar = [ x1-1;
        x2-1;
        x3-1;
        x4-1;
        x5-1;
        x6-1;
        -x1-10;
        -x2-10;
        -x3-10;
        -x4-0.1;
        -x5-0.1;
        -x6-0.1;
];       
 

%% ------------------------------------------------------------------------
% Do not modify below this line 
%-------------------------------------------------------------------------

if isempty(Phi)
    Phi = zeros(0,1);
    matlabFunction(Phi, 'file', 'valPhi','vars',x);
else
    matlabFunction(Phi, 'file', 'valPhi','vars',{x}); 
end

if isempty(L)
    L = zeros(0,1);
    matlabFunction(L, 'file', 'valL','vars',[x, u]); 
else
    matlabFunction(L, 'file', 'valL','vars',{x, u}); 
end

if isempty(f)
    f = zeros(0,1);
    matlabFunction(f, 'file', 'valf','vars',[x, u]); 
else
    matlabFunction(f, 'file', 'valf','vars',{x, u}); 
end

if isempty(C)
    C = zeros(0,1);
    matlabFunction(C, 'file', 'valC','vars',[x, u]);   
else
    matlabFunction(C, 'file', 'valC','vars',{x, u});
end

if isempty(Cbar)
    Cbar = zeros(0,1);
    matlabFunction(Cbar, 'file', 'valCbar','vars',x); 
else
    matlabFunction(Cbar, 'file', 'valCbar','vars',{x}); 
end


%%
if isempty(Phi)
    Phix = zeros(0,1);
    Phixx = zeros(0,1);
    
    matlabFunction(Phix, 'file', 'DPhiDx','vars',x);
    matlabFunction(Phixx, 'file', 'D2PhiDxx','vars',x); 
else
    for j=1:size_x
        Phix(1,j) = diff(Phi,(sprintf('x%1d',j)));
    end
    matlabFunction(Phix, 'file', 'DPhiDx','vars',{x});   

    for j = 1:size_x
        for k = 1:size_x
            Phixx(j,k) = diff(Phix(1,k),(sprintf('x%1d',j)));
        end
    end
    matlabFunction(Phixx, 'file', 'D2PhiDxx','vars',{x}); 
end

%%

if isempty(L)
    Lx = zeros(0,1);
    Lu = zeros(0,1);
    Lxx = zeros(0,1);
    Lxu = zeros(0,1);
    Luu = zeros(0,1);
    
    matlabFunction(Lx, 'file', 'DLDx','vars',[x, u]);
    matlabFunction(Lu, 'file', 'DLDu','vars',[x, u]);
    matlabFunction(Lxx, 'file', 'D2LDxx','vars',[x, u]);
    matlabFunction(Lxu, 'file', 'D2LDxu','vars',[x, u]); 
    matlabFunction(Luu, 'file', 'D2LDuu','vars',[x, u]);
else
    for j=1:size_x
        Lx(1,j) = diff(L,(sprintf('x%1d',j)));
    end
    matlabFunction(Lx, 'file', 'DLDx','vars',{x, u});   

    for j=1:size_u
        Lu(1,j) = diff(L,(sprintf('u%1d',j)));
    end
    matlabFunction(Lu, 'file', 'DLDu','vars',{x, u});   

    for j = 1:size_x
        for k = 1:size_x
            Lxx(j,k) = diff(Lx(1,k),(sprintf('x%1d',j)));
        end
    end
    matlabFunction(Lxx, 'file', 'D2LDxx','vars',{x, u}); 

    for j = 1:size_u
        for k = 1:size_x
            Lxu(j,k) = diff(Lx(1,k),(sprintf('u%1d',j)));
        end
    end
    matlabFunction(Lxu, 'file', 'D2LDxu','vars',{x, u}); 

    for j = 1:size_u
        for k = 1:size_u
            Luu(j,k) = diff(Lu(1,k),(sprintf('u%1d',j)));
        end
    end
    matlabFunction(Luu, 'file', 'D2LDuu','vars',{x, u}); 

end


%%
size_f = size(f,1);
save size_f;
if isempty(f)
    fx = zeros(0,1);
    fu = zeros(0,1);
    fxx = zeros(0,1);
    fxu = zeros(0,1);
    fuu = zeros(0,1);
    
    matlabFunction(fx, 'file', 'DfDx','vars',[x, u]);
    matlabFunction(fu, 'file', 'DfDu','vars',[x, u]);
    matlabFunction(fxx, 'file', 'D2fDxx','vars',[x, u]);
    matlabFunction(fxu, 'file', 'D2fDxu','vars',[x, u]); 
    matlabFunction(fuu, 'file', 'D2fDuu','vars',[x, u]);
else
    for j=1:size_x
        for k = 1:size_x
            fx(j,k) = diff(f(j),(sprintf('x%1d',k)));
        end 
    end
    matlabFunction(fx, 'file', 'DfDx','vars',{x, u});   

    for j=1:size_x
        for k = 1:size_u
            fu(j,k) = diff(f(j),(sprintf('u%1d',k)));
        end 
    end
    matlabFunction(fu, 'file', 'DfDu','vars',{x, u});   

    for n=1:size_x
        for j = 1:size_x
            for k = 1:size_x
                fxx(j,k,n) = diff(fx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(fxx, 'file', 'D2fDxx','vars',{x, u}); 

    for n=1:size_x
        for j = 1:size_u
            for k = 1:size_x
                fxu(j,k,n) = diff(fx(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(fxu, 'file', 'D2fDxu','vars',{x, u}); 

    for n=1:size_x
        for j = 1:size_u
            for k = 1:size_u
                fuu(j,k,n) = diff(fu(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(fuu, 'file', 'D2fDuu','vars',{x, u}); 
end


%%
size_C = size(C,1);
save size_C;
if size_C == 0
    Cx = zeros(0,1);
    Cu = zeros(0,1);
    Cxx = zeros(0,1);
    Cxu = zeros(0,1);
    Cuu = zeros(0,1);
    
    matlabFunction(Cx, 'file', 'DCDx','vars',[x, u]);
    matlabFunction(Cu, 'file', 'DCDu','vars',[x, u]);
    matlabFunction(Cxx, 'file', 'D2CDxx','vars',[x, u]);
    matlabFunction(Cxu, 'file', 'D2CDxu','vars',[x, u]); 
    matlabFunction(Cuu, 'file', 'D2CDuu','vars',[x, u]);
else
    for j=1:size_C
        for k = 1:size_x
            Cx(j,k) = diff(C(j),(sprintf('x%1d',k)));
        end
    end
    matlabFunction(Cx, 'file', 'DCDx','vars',{x, u});   

    for j=1:size_C
        for k = 1:size_u
            Cu(j,k) = diff(C(j),(sprintf('u%1d',k)));
        end
    end
    matlabFunction(Cu, 'file', 'DCDu','vars',{x, u});   

    for n=1:size_C
        for j = 1:size_x
            for k = 1:size_x
                Cxx(j,k,n) = diff(Cx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(Cxx, 'file', 'D2CDxx','vars',{x, u}); 

    for n=1:size_C
        for j = 1:size_u
            for k = 1:size_x
                Cxu(j,k,n) = diff(Cx(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(Cxu, 'file', 'D2CDxu','vars',{x, u}); 

    for n=1:size_C
        for j = 1:size_u
            for k = 1:size_u
                Cuu(j,k,n) = diff(Cu(n,k),(sprintf('u%1d',j)));
            end
        end
    end
    matlabFunction(Cuu, 'file', 'D2CDuu','vars',{x, u}); 
end

%%
size_Cbar = size(Cbar,1);
save size_Cbar;
if size_Cbar == 0
   Cbarx = zeros(0,1);
   Cbarxx = zeros(0,1);
   matlabFunction(Cbarx, 'file', 'DCbarDx','vars',x);
   matlabFunction(Cbarxx, 'file', 'D2CbarDxx','vars',x); 
else
    for j=1:size_Cbar
        for k = 1:size_x
            Cbarx(j,k) = diff(Cbar(j),(sprintf('x%1d',k)));
        end
    end
    matlabFunction(Cbarx, 'file', 'DCbarDx','vars',{x});   

    for n=1:size_Cbar
        for j = 1:size_x
            for k = 1:size_x
                Cbarxx(j,k,n) = diff(Cbarx(n,k),(sprintf('x%1d',j)));
            end
        end
    end
    matlabFunction(Cbarxx, 'file', 'D2CbarDxx','vars',{x}); 
end

return

