function [X, C, I, signals] = Model(inp,par)
%% vehicle parameters
cd       = 0.373;%风阻系数
Af       = 2.58;%有效迎风面积
ro       = 1.205;%空气密度
M        = 1870;%整车质量
f        = 0.011;%滚动阻力系数
g        = 9.8;%重力加速度
itat     = 0.94;%机械效率
rw       = 0.364;%车轮半径
p1       = 0.00385165576001348;
p2       = 1.19220483871151E-6;
p3       = 2.95706278924791E-10;
p4       = 5.45299140436618E-5;
p5       = -1.59515599116976E-9;
p6       = 1.74678055479864E-7; %z = p1*x+p2*x^2+p3*x^3+p4*x*y+p5*x^2*y+p6*x*y^2

%% state update

if inp.U{2}== 4
    Ig = 4.66;
elseif inp.U{2}== 5
    Ig = 3.66;
else
    Ig = 2.63;
end
func      = itat*Ig.*inp.U{1}./(rw.*M.*inp.X{1}) - cd.*Af.*ro.*inp.X{1}./2./M - g.*(f.*cosd(inp.W{1}) + sind(inp.W{1}))./inp.X{1};
X{1}      = inp.Ts.*func + inp.X{1};

%% cost
% T         = (rw.*(cd.*Af.*ro.*inp.X{1}.^2*0.5 + M.*g.*(cos(inp.W{1}) + sin(inp.W{1})))./(itat.*Ig));
% Tmax      = (rw.*(cd.*Af.*ro.*20.^2.*0.5 + M.*g.*(cos(inp.W{1}) + sin(inp.W{1})))./(itat.*Ig));
% Jfuel_max = (p1.*kk + p2.*kk.^2.*20 + p3.*kk.^3.*20.^2 + p4.*kk.*Tmax + p5.*kk.^2.*20.*Tmax + p6.*kk.*Tmax.^2);
% Tmin      = (rw.*(cd.*Af.*ro.*10.^2.*0.5 + M.*g.*(cos(inp.W{1}) + sin(inp.W{1})))./(itat.*Ig));
% Jfuel_min = (p1.*kk + p2.*kk.^2.*10 + p3.*kk.^3.*10.^2 + p4.*kk.*Tmin + p5.*kk.^2.*10.*Tmin + p6.*kk.*Tmin.^2);
kk        = 30*Ig/(3.14*rw);
Jfuel     = (p1.*kk + p2.*kk.^2.*inp.X{1} + p3.*kk.^3.*inp.X{1}.^2 + p4.*kk.*inp.U{1} + p5.*kk.^2.*inp.X{1}.*inp.U{1} + p6.*kk.*inp.U{1}.^2);
Jtime     = 1./inp.X{1};
alpha     = 0.05;
% C{1}      = inp.Ts.*sqrt(((Jfuel-inp.W{2})/(Jfuel_max-inp.W{2})).^2 + ((5./inp.X{1}-0.25)./0.25).^2);
C{1}      = inp.Ts.*(alpha*Jfuel + (1-alpha)*(Jtime));
 
%% feasibility
I  =  0;

%% Output Signals
signals.u1= inp.U{1};
signals.u2= inp.U{2};
signals.x1= inp.X{1};
signals.w1= inp.W{1};
signals.c1= Jfuel;
signals.c2= Jtime;

% if numel(find(I==0))==0
%     keyboard
% end
