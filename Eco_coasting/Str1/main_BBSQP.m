clc
clear all
%% https://optimization.mccormick.northwestern.edu/index.php/Branch_and_bound_(BB)
tfinal = 200;
Horizon = 50;
Samp_dis = 10;
Num = Horizon/Samp_dis;
nsim = ceil(tfinal/Samp_dis);
size_x = 1;                   % number of states
size_u = 3;                   % number of control input
%% obtain first optimal solution using SQP
[size_C, size_Cbar]=generate_functions(size_x,size_u);
load size_C;
load size_Cbar;
opt.itmax = 100; % maximum iteration of the SQP
opt.tol = 10e-4; % convergence tolerance of SQP
opt.gamma = 10;   % weight of slack variable 
opt.maxNodes = 100000; % maximum modes
opt.branchMethod = 3;
opt.branchCriteria = 1;
opt.intTol = 1e-4;
opt.display = 'off';
opt.iterplot = false;
opt.warmStart = false;
opt.gap = false;
if opt.gap
    opt.Delta = 1e-8; %Stopping tolerance of the gap (f_integer-f_lp)/(f_integer+f_lp)
else
    opt.Delta = -inf;
end
xsim(:,1) = 20; % Initial state (velocity)
%% Initialization
uu = [10,10,0.5];
uu = repmat(uu,1,Num);
xx = xsim(:,1);
xx = repmat(xx,1,Num);
s = zeros(Num,1);
X2 = zeros(size_x,nsim);
t2 = zeros(nsim,1);
iter = zeros(nsim,1);
U2 = zeros(size_u,nsim);
z = [reshape(uu,[size_u*Num,1]);reshape(xx,[size_x*Num,1]);s];
lambda = zeros(size_x*Num,1);   % initial guess of dual variables
nu = zeros(Num*(size_C+size_Cbar+1),1);   % initial guess of dual variables
km.z = z;
km.l = lambda;
km.v = nu;
yidx = false(2,Num);
yidx = [yidx;true(1,Num)];
yidx = [reshape(yidx,[size_u*Num,1]);false(2*Num,1)];
warmpath = nan(Num,4);
limit_num = 5;
z_hist = [1;1;1;1;1]; % initial history signal of fuel cut-off, 1- fuel cut-off is inactivated
%% slope and reference speed resampling
aa = zeros(nsim + Num,1);
load dist100
load slope
for i = 1:nsim + Num
    for jj = 1:1:(size(dist100,2)-1)
        if (dist100(jj) < i*Samp_dis ) && ( i*Samp_dis <= dist100(jj+1))
            aa(i) = slope(jj);
        end
    end
end

load dist.mat
load speed_ref.mat
sped_ref = 21.623*ones(nsim + Num,1);
for k = 1:nsim + Num
    if k == 1
        sped_ref(k) = 19.1305;
    else
        for jj = 2:1:(size(dist,2))
            if (dist(jj-1) < k*Samp_dis) && (k*Samp_dis <= dist(jj))
                sped_ref(k) = speed_ref(jj-1);
            end
        end
    end
end
%% Simulation
for i =1:nsim-1
    a = aa(i:i+(Num-1))'; % Give the slope profile of the current predictive horizon
    ref = sped_ref(i:i+(Num-1))';
    p = [a;ref]; % Parameter comprise slope data and reference speed vector
    tic
    [u_sqp, o, s, d, f, it] = BandB(xsim(:,i),Num,km,opt,yidx,size_x,size_u,p,warmpath);
    time(i) = toc;
    uu = reshape(u_sqp(1:size_u*Num,1),[size_u,Num]);
    z_curr= uu(3,:);
    z_hist(1) = [];
    z_hist(limit_num) = z_curr(1);
    warmpath = [o s d f]; % Optimal path at current sampling interval
    xx = reshape(u_sqp(size_u*Num+1:size_u*Num+size_x*Num,1),[size_x,Num]);
    xsim(:,i+1) = valf(xsim(:,i)',uu(:,1)',a(1));
    U2(:,i) = uu(:,1)';
    t2(i) = (i-1)*Samp_dis;
    iter(i) = it;
end
t2(1) = 0; % fill the first element of time sequence for plotting
save Str123_50.mat
%% Plotting stuff

% play with the colours a bit for prettiness
co =  [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];

% co = flipud(co);
co = circshift(co,1,1);

set(groot,'defaultAxesColorOrder',co);


figure();
load withws100.mat
ax(1) = subplot(3,2,1);
hold on
plot(t2(1:nsim-1),aa(1:nsim-1),'k');
% plot(t2,-0*ones(size(t1)),'k-.');
% plot(t1,T1,'--');
% legend('\theta_1','\theta_2','\theta_3','\theta_{lb}');
ylabel('坡度 [rad]');
box on;

ax(2) = subplot(3,2,2);
hold on
load withws100.mat
plot(t2(1:nsim-1),sped_ref(1:nsim-1),'k--');
hold on
% plot(t1,ones(size(t1)),'k.-',t1,-ones(size(t1)),'-.k')
plot(t2(1:nsim-1),xsim(1:nsim-1)','r--');
hold on
load withoutws100.mat
plot(t2(1:nsim-1),xsim(1:nsim-1)','g');
ylabel('车速 [m/s]');
legend('参考车速','真实车速');
box on;

ax(3) = subplot(3,2,3);
hold on
load withws100.mat
plot(t2(1:nsim-1),U2(1,1:nsim-1),'r--');
hold on
load withoutws100.mat
plot(t2(1:nsim-1),U2(1,1:nsim-1),'g');
% legend('\omega_1','\omega_2','\omega_3');
ylabel('转矩 [Nm]');
box on;

ax(4) = subplot(3,2,4);
hold on
load withws100.mat
% semilogy(t1,res1,'--',t2,res2);
plot(t2(1:nsim-1),U2(2,1:nsim-1),'r--');
hold on
load withoutws100.mat
plot(t2(1:nsim-1),U2(2,1:nsim-1),'g');
% legend('Unconstrained','Constrained');
ylabel('制动 [Nm]');
box on;

ax(5) = subplot(3,2,5);
hold on
load withws100.mat
% semilogy(t1,res1,'--',t2,res2);
stairs(t2(1:nsim-1),U2(3,1:nsim-1),'r--');
hold on
load withoutws100.mat
stairs(t2(1:nsim-1),U2(3,1:nsim-1),'g');
% legend('Unconstrained','Constrained');
xlabel('距离 [m]');
ylabel('断油 [-]');
box on;

ax(6) = subplot(3,2,6);
hold on
load withws100.mat
% semilogy(t1,res1,'--',t2,res2);
stairs(t2(1:nsim-1),iter(1:nsim-1),'r--');
hold on
load withoutws100.mat
stairs(t2(1:nsim-1),iter(1:nsim-1),'g');
% legend('Unconstrained','Constrained');
xlabel('距离 [m]');
ylabel('节点个数[-]');
box on;

linkaxes(ax,'x');