clc
clear all

%% https://optimization.mccormick.northwestern.edu/index.php/Branch_and_bound_(BB)
tfinal = 200;
ts = 2;
nsim = ceil(tfinal/ts);

%% model parameter
size_x = 6;                   % number of states
size_u = 3;                   % number of control input
% x_init = 20;                  % initial states
N = 3;                        % prediction horizon
[size_C, size_Cbar]=generate_functions(size_x,size_u);
load size_C;
load size_Cbar;
uu = zeros(size_u,N);
xx = zeros(size_x,N);
s = zeros(N,1);
X2 = zeros(6,nsim);
t2 = zeros(nsim,1);
U2 = zeros(3,nsim);
xsim(:,1) = x0;
z = [reshape(uu,[size_u*N,1]);reshape(xx,[size_x*N,1]);s];
lambda = zeros(size_x*N,1);   % initial guess of dual variables
nu = zeros(N*(size_C+size_Cbar+1),1);   % initial guess of dual variables
km.z = z;
km.l = lambda;
km.v = nu;
itmax = 5;
tol = 10e-4;
t2(1) = 0;
for i = 1:nsim-1
%     x0 = xsim(:,i);
    [u_sqp,xtra,fval,it,res1(i),time(i)] = seq_quad_prog(xsim(:,i),km,itmax,tol,N,size_x,size_u,gamma);
    xsim(:,i+1) = valf(xsim(:,i)',u_sqp');
    km.z = xtra.z;
    km.l = xtra.l;
    km.v = xtra.v;
    U2(:,i) = u_sqp;
    t2(i) = (i-1)*ts;
end
U2(:,nsim) = U2(:,nsim-1);
W2 = xsim(1:3,:);
T2 = xsim(4:6,:);
res1(nsim) = 0;
time(nsim) = 0;
%% plot
figure();
ax(1) = subplot(2,2,1);
hold on
plot(t2,T2);
plot(t2,-0*ones(size(t2)),'k-.');
% plot(t1,T1,'--');

legend('\theta_1','\theta_2','\theta_3','\theta_{lb}');
ylabel('Euler Angles [rad]');

ax(2) = subplot(2,2,2);
hold on
plot(t2,U2);
plot(t2,ones(size(t2)),'k.-',t2,-ones(size(t2)),'-.k')
% plot(t1,U1,'--');
ylabel('Torque [Nm]');

legend('u_1','u_2','u_3','u_{ub}','u_{lb}');

ax(3) = subplot(2,2,3);
hold on
plot(t2,W2);
% plot(t1,W1,'--');
legend('\omega_1','\omega_2','\omega_3');
ylabel('Angular Velocity [rad/s]');
xlabel('Time [s]')

ax(4) = subplot(2,2,4);
% semilogy(t1,res1,'--',t2,res2);
semilogy(t2,res1,'--');
legend('Unconstrained','Constrained');
xlabel('Time [s]');
ylabel('||F_{NR}||');

% load time_FBsstab
% figure
% plot(t2(1:99),time(1:99),'k');
% hold on
% load time_quadprog
% plot(t2(1:99),time(1:99),'r');
% xlabel('Time [s]');
% ylabel('Computational time [s]');