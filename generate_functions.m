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

% define sampling distance
ts       = 1;
beta     = 0.0001;
a        = 0;% 车辆固有参数
cd       = 0.373;%风阻系数
Af       = 2.58;%有效迎风面积
ro       = 1.205;%空气密度
M        = 1870;%整车质量
f        = 0.011;%滚动阻力系数
g        = 9.8;%重力加速度
itat     = 0.94;%机械效率
rw       = 0.364;%车轮半径
Ig       = 4.66;%4挡位传动比
% w       = 30*Ig/(3.14*rw*1000); %为什么要除1000
w       = 30*Ig/(3.14*rw); 
theta1   = 0;
% Terminal cost
Phi = [];

% Incremental cost
L = beta*(p1*w*u3 + p2*w^2*x + p3*w*u1^2) + (1-beta)*(1/x);

% Equations of system dynamics
f = x + ts*( (itat*Ig*u1*u3-u2)/(rw*M*x) - cd*Af*ro*x/2/M - g*(f + sind(theta1))/x );

% Input and state mixed constriants
C =  [-u1;
      -u2;
      -u3;
      u1-120;
      u2-500;
      u3-1];

% State only constraints 
Cbar = [-x+15;
        x-25];       
   
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
    Phix = zeros(1,1);
    Phixx = zeros(1,1);
    
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

