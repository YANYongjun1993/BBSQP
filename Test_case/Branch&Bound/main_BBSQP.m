clc
clear all
%% https://optimization.mccormick.northwestern.edu/index.php/Branch_and_bound_(BB) 
%% refer to: the first test problem in An outer-approximation algorithm for a class of mixed-integer nonlinear programs
tfinal = 200;
ts = 2;
nsim = ceil(tfinal/ts);

Horizon = 10;
Samp_dis = 1;
N = Horizon/Samp_dis;
size_z= 6;                   % number of states
%% obtain first optimal solution using SQP
% [size_C]=generate_functions(size_z);
load size_C;
opt.itmax = 500; % maximum iteration of the SQP
opt.tol = 10e-4; % convergence tolerance of SQP
opt.gamma = 10;   % weight of slack variable 
opt.maxNodes = 100000; % maximum modes
opt.branchMethod = 3;
opt.branchCriteria = 1;
opt.intTol = 1e-6;
opt.display = 'iter';
opt.iterplot = false;
opt.solver = 'linprog';
opt.Delta = 1e-8; %Stopping tolerance of the gap (f_integer-f_lp)/(f_integer+f_lp)
xsim(:,1) = 20; % Initial state (velocity)

z = [1,1,0,0.5,0.5,0.5];
nu = zeros(size_C,1);   % initial guess of dual variables
z0 = z;
v0 = nu;
yidx = [1;1;1;0;0;0];
for i = 1:nsim-1
    u_sqp = BandB(z0,v0,opt,yidx);
    xsim(:,i+1) = valf(xsim(:,i)',u_sqp');
    km.z = xtra.z;
    km.l = xtra.l;
    km.v = xtra.v;
    U2(:,i) = u_sqp;
    t2(i) = (i-1)*ts;
end
t2(1) = 0; % fill the first element of time sequence for plotting
%% B&B
% function [x,out] = BandB(v0,N)
%     gamma = 10;   % weight of slack variable 
%     uu = [10,100,0.5]';
%     uu = repmat(uu,1,N);
%     xx = v0;
%     xx = repmat(xx,1,N);
%     s = zeros(N,1);
%     X2 = zeros(6,nsim);
%     t2 = zeros(nsim,1);
%     U2 = zeros(3,nsim);
%     xsim(:,1) = x0;
%     z = [reshape(uu,[size_u*N,1]);reshape(xx,[size_x*N,1]);s];
%     lambda = zeros(size_x*N,1);   % initial guess of dual variables
%     nu = zeros(N*(size_C+size_Cbar+1),1);   % initial guess of dual variables
%     km.z = z;
%     km.l = lambda;
%     km.v = nu;
%     equ_con = zeros(size_u*N,1); % change the inequality constraint into the equality constraint
%     [u_sqp,xtra,fval,it,res1(i),time(i),flag] = seq_quad_prog(xsim(:,i),km,itmax,tol,N,size_x,size_u,gamma,equ_con); % father node
%     que = cell(0); % create a queue for the branch candidates
%     t = tree('root');
%     opt_fval = fval; % optimal objective value so far
%     for i = 3:3:3*N
%         sum = sum + u_sqp(i);
%     end
%     if (flag == 0) 
%         fprintf('The Problem is infeasiable,fathom this point\n');
%     elseif sum == fix(sum)
%         u_opt = u_sqp;
%         fprintf('The optimal integer value are found,fathom this point\n');
%     else
%         for i = 3:3:3*N % find the first value to branch, it need to be optimated later
%             if 0 < u_sqp(i) < 1
%                 [ t node1 ] = t.addnode(1, 0);
%                 [ t node2 ] = t.addnode(1, 1);
%                 que(:,1) = {i 0}; 
%                 que(:,2) = {i 1}; 
%                 break
%             end
%         end  % initination is finished so far
%     end
%     while( ~isempty(que))
%         km = xtra; % warm start
%         if que{2, 1} == 0 % tell the first value of the queue
%             km.z(i) = 0; % set the initial guess of the control variable as 0
%             equ_con(2*i) = 1;
%         else
%             km.z(i) = 1; % set the initial guess of the control variable as 0
%             equ_con(2*i) = -1;
%         end
%         [u_sqp,xtra,fval,it,res1(i),time(i)] = seq_quad_prog(xsim(:,i),km,itmax,tol,N,size_x,size_u,gamma,equ_con); % father node
%         for i = 3:3:3*N
%             sum = sum + u_sqp(i);
%         end
%         if (flag == 0)           
%             fprintf('The Problem is infeasiable,fathom this point\n');
%         elseif fval < opt_fval
%             fprintf('Worse than incumbent,fathom this point\n');
%         elseif sum == fix(sum)
%             u_opt = u_sqp;
%             fprintf('The optimal integer value are found,fathom this point\n');
%         else % branch
%             for i = 3:3:3*N
%                 if 0 < u_sqp(i) < 1
%                     num = size(que,2);
%                     [ t eval(['z',num2str(i),'=','0']) ] = t.addnode(0);
%                     [ t eval(['z',num2str(i),'=','1']) ] = t.addnode(1);
%                     que(:,num+1) = {i 0}; 
%                     que(:,num+2) = {i 1}; 
%                     break
%                 end
%             end  
%         end
%         que(:,1) = []; % delete the first element of the queue
%     end
% end