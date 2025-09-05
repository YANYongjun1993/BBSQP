%Mf.m
% È¼ÓÍÏûºÄÂÊ
function mf=Mf2(n,p,z,v0)
h = 10;
v0 = v0/3.6;
p1 =   0.2159;% µ¡ËÙÓÍºÄm/s
p2 =   0.005676;
p3 =   0.0004349; 
p4 =   8.899e-07; % f(x,y) = p1 + p2*x*y + p3*x^2*y + p4*x*y^2
if z == 0
    mf = 0;
else
    mf = (p1 + p2*n*p + p3*n^2*p + p4*n*p.^2)*(1/v0)*h;
end