% clc
% clear all
load('dist_partial.mat');
load('slope_partial.mat');
global theta
Utopia   = zeros(2,size(dist_partial,2));
Knee     = zeros(2,size(dist_partial,2));
% for i = 1:1:size(dist_partial,2)
for i = 1:1
    theta = slope_partial(1,i);
    k = 1;
    [min1,minfn1] = fminbnd(@(x)pickindex(x,k),15,25);
    k = 2;
    [min2,minfn2] = fminbnd(@(x)pickindex(x,k),15,25);
    goal = [minfn1,minfn2];
    Utopia(:,i) = goal;
    nf = 2; % number of objective functions
    N = 500; % number of points for plotting
    onen = 1/N;
    x = zeros(N+1,1);
    f = zeros(N+1,nf);
    fun = @simple_mult;
    x0 = 5;
    options = optimoptions('fgoalattain','Display','off');
    for r = 0:N
        t = onen*r; % 0 through 1
        weight = [t,1-t];
        [x(r+1,:),f(r+1,:)] = fgoalattain(fun,x0,goal,weight,...
            [],[],[],[],[],[],[],options);
    end
    min_dis = 100000000; %knee point
    j = 1;
    for ii = 1:N+1
        dis = sqrt((f(ii,1)-minfn1)^2 + (f(ii,2)-minfn2)^2);
        if dis<min_dis
            min_dis = dis;
            j = ii;
        end
    end
    Knee(:,i) = [f(j,1),f(j,2)];
end
