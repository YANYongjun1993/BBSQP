function [x_best, i] = BandB(x0, N, km, o, yidx, size_x, size_u, p, z_hist, limit_num)
%BandB
% Branch and bound algorithm for mixed integer nonlinear programs
% x_best = BandB(c,A,b,Aeq,beq,lb,ub,yidx,[options])
%
%   min c'*x s.t. A*x <= b Aeq*x == beq lb <= x <= ub x(yidx) in [0,1,2,...]
% 
%   where yidx is a logical index vector such that
%   x(yidx) are the integer variables
% 
% Input:
% c,A,b,Aeq,beq,lb,ub - problem matrices
% yidx - a logical index vector such that y = x(yidx)
% options - (optional) mipoptions 
% 
% Output:
% x_best - best integer solution found
% 
% See Also:
% mip_gridoptim, gridoptions
% 
% Programmed by Yongjun Yan 2021
%  


%Assume no initial best integer solution
%Add your own heuristic here to find a good incumbent solution, store it in
%f_best,y_best,x_best
f_best = inf;
y_best = [];
x_best = [];

%Variable for holding the objective function variables of the lp-relaxation
%problems
f = inf(o.maxNodes,1);
f(1) = 0;
fs = inf;
numIntSol = double(~isempty(y_best));

%Set of problems
S = nan(sum(yidx),1);
D = zeros(sum(yidx),1);

%The priority in which the problems shall be solved
priority = [1];
%The indices of the problems that have been visited
visited = nan(o.maxNodes,1);
% History binary variable
k = 0;
for i = limit_num:-1:1
    if abs(z_hist(i)) <= o.intTol
        k = k+1;
    else
        break
    end
end
if k ~= 0
    for i = 1:limit_num-k
            S(i) = 0;
            D(i) = 1;
    end
end
if abs(z_hist(limit_num)) <= o.intTol
    s = zeros(sum(yidx),1);
    d = ones(sum(yidx),1);
    S = [S s];
    D = [D d];
    for i =N:-1:limit_num-k+1
        s(i) = 1;
        d(i) = -1;
        S = [S s];
        D = [D d];
    end
else
    s = ones(sum(yidx),1);
    d = -1*ones(sum(yidx),1);
    S = [S s];
    D = [D d];
    for i =N:-1:1
        s(i) = 0;
        d(i) = 1;
        S = [S s];
        D = [D d];
    end
end
for i = 2:size(S,2)
    s = S(:,i);
    d = D(:,i);
    [x,this_f,flag] = seq_quad_prog(x0,km,o,N,size_x,size_u,s,d,yidx,p); 
    if flag==0
        %infeasible, dont branch
        if i==1
            error('BBSQP: Infeasible initial SQP problem. Try another SQP solver.')
        end
        f(i) = inf;
    elseif flag==1        
        f(i) = this_f;  
        if this_f<f_best
            %find the integer variables among the control sequence
            y = x(yidx);
            numIntSol = numIntSol+1;
            f_best = this_f;
            y_best = round(x(yidx));
            x_best = x;
        end
    end
end
end