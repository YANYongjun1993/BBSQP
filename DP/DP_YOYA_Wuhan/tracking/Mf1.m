%Mf.m
% È¼ÓÍÏûºÄÂÊ
function mf=Mf1(n,p,z,v0,v_ref)
h = 10;
v0 = v0/3.6;
v_ref = v_ref/3.6;
beta = 0.7;
p1 =   0.2159;% µ¡ËÙÓÍºÄm/s
p2 =   0.005676;
p3 =   0.0004349; 
p4 =   8.899e-07; % f(x,y) = p1 + p2*x*y + p3*x^2*y + p4*x*y^2
if p >= 0
    mf = (beta*(p1 + p2*n*p + p3*n^2*p + p4*n*p.^2)*(1/v0) + (1-beta)*(v0-v_ref)^2)*h;
elseif p == -30
    mf = (1-beta)*(v0-v_ref)^2*h;
else
    mf = (beta*p1*(1/v0) + (1-beta)*(v0-v_ref)^2)*h;
end