function [x_best, od_best, s_best, d_best, fod_best,i] = BandB(x0,N,km,o,yidx,size_x,size_u,p,warmpath)
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
Od = nan(sum(yidx),1);
Fod = nan(sum(yidx),1);
Fod_last = nan(sum(yidx),1);
D = zeros(sum(yidx),1);

%The priority in which the problems shall be solved
priority = [1];
%The indices of the problems that have been visited
visited = nan(o.maxNodes,1);

%Plot each iteration?
i=0;
if o.iterplot
    figure;
    hold on;
    title('Bounds')
    xlabel('Iteration')
    ylabel('Obj. fun. val')
end

%WarmStart tree have not been built
Status_WarmStart_construction = 0;
Length_WarmStart = 0;
if  sum(isnan(warmpath(:,1))) == N
    warmpath = [];
end
%% Branch and bound loop
while i==0 || isinf(f_best) || (~isempty(priority) &&  ((f_best-min(fs(priority)))/abs(f_best+min(fs(priority))) > o.Delta) &&  i<o.maxNodes)
    %Is the parent node less expensive than the current best
    if i==0 || fs(priority(1))<f_best
        %Solve the NLP problem
        i=i+1;
        s = S(:,priority(1));
        od1 = Od(:,priority(1));
        d = D(:,priority(1));
        fod_last = Fod_last(:,priority(1));
        [x,this_f,flag] = seq_quad_prog(x0,km,o,N,size_x,size_u,s,d,yidx,p); 
        
        %Visit this node
        visited(i) = priority(1);
        priority(1) = [];
        if flag==0
            %Infeasible, dont branch
            if i==1
%                 error('BBSQP: Infeasible initial SQP problem. Try another SQP solver.')
            disp('BBSQP: Infeasible initial SQP problem. Try another SQP solver.');
            end
            f(i) = inf;
            i = i - 1;
        elseif flag==1
            %convergence
            f(i) = this_f;  
            if i > 1
                fod1 = fod_last;
                fod1 = circshift(fod1,-1);
                fod1(N) = NaN;
                sm = sum(~isnan(s));
                fod1(sm) = this_f;
                Fod = [Fod fod1]; 
            end
            if this_f<f_best
                %find the integer variables among the control sequence
                y = x(yidx);
                %fractional part of the integer variables -> diff
                diff = abs(round(y)-y);
                if all(diff<o.intTol)
                    %all fractions less than the integer tolerance
                    %we have integer solution
                    numIntSol = numIntSol+1;
                    f_best = this_f;
                    y_best = round(x(yidx));
                    x_best = x;
                    od_best = od1;
                    s_best = s;
                    d_best = d;
                    if i > 1
                        fod_best = fod1;
                    else
                        fod_best = fod_last;
                        break
                    end
                else
                    %% WarmStart 
                    if ~isempty(warmpath) && o.warmStart && ~Status_WarmStart_construction
                        %shift
                        od_best_ws = warmpath(:,1)-1;
                        s_best_ws = circshift(warmpath(:,2),-1);
                        s_best_ws(N) = NaN;
                        d_best_ws = circshift(warmpath(:,3),-1);
                        d_best_ws(N) = 0;
                        fod_best_ws = warmpath(:,4);
                        %remove the variable on stage 0 after shifting 
                        for ii = 1:N
                            if od_best_ws(ii) == 0
                                od_best_ws(ii) = [];
                                od_best_ws(N) = NaN;
                                fod_best_ws(ii) = [];
                                fod_best_ws(N) = NaN;
                            end
                        end
                        for ii = N:-1:1
                            if ~isnan(fod_best_ws(ii))
                                Fod_last_ws = fod_best_ws;
                                Fod_last_ws(ii) = NaN;
                                Fod_last_ws = circshift(Fod_last_ws,1);
                                Fod_last_ws(1) = this_f;
                                break
                            end
                        end
                        % store the longest path to the warm start list
                        ii = N;
                        while sum(isnan(od_best_ws))~=N
                            if ~isnan(od_best_ws(ii)) 
                                S = [S s_best_ws];
                                Od = [Od od_best_ws];
                                D = [D d_best_ws];
                                fs = [fs warmpath(ii,4)]; 
                                Fod_last = [Fod_last Fod_last_ws];
                                break
                            end
                        ii = ii -1;
                        end
                        % store the child nodes to the warm start list
                        while sum(isnan(od_best_ws))~=N
                            if s_best_ws(od_best_ws(ii)) == 1
                                s_best_ws(od_best_ws(ii)) = 0;
                                d_best_ws(od_best_ws(ii)) = 1;
                            else
                                s_best_ws(od_best_ws(ii)) = 1;
                                d_best_ws(od_best_ws(ii)) = -1;
                            end
                            S = [S s_best_ws];
                            Od = [Od od_best_ws];
                            D = [D d_best_ws];
                            Fod_last = [Fod_last Fod_last_ws];
                            fs = [fs warmpath(ii,4)]; 
                            s_best_ws(od_best_ws(ii)) = NaN;
                            d_best_ws(od_best_ws(ii)) = 0;
                            od_best_ws(ii) = NaN;
                            Fod_last_ws(ii) = NaN;
                            ii = ii -1;
                        end
                        % WarmStart tree construction is finished
                        Status_WarmStart_construction = 1;
                        Length_WarmStart = size(S,2);
                        nsold = 1;
                    end
                    if i ~= 1 || isempty(warmpath) || ~o.warmStart
                        if o.branchCriteria==1
                            %branch on the most fractional variable
                            [maxdiff,branch_idx] = max(diff,[],1);
                        elseif o.branchCriteria==2
                            %branch on the least fractional variable
                            diff(diff<o.intTol)=inf;
                            [mindiff,branch_idx] = min(diff,[],1);
                        elseif o.branchCriteria==3
                            %branch on the variable with highest cost
                            cy = c(yidx);
                            cy(diff<o.intTol)=-inf;
                            [maxcost,branch_idx] = max(cy,[],1);
                        elseif o.branchCriteria==4
                            %branch on the variable with lowest cost
                            cy = c(yidx);
                            cy(diff<o.intTol)=inf;
                            [mincost,branch_idx] = min(cy,[],1);  
                        else
                            error('BBSQP: Unknown branch criteria.')
                        end
                        %Branch into two subproblems
                        s1 = s;
                        s2 = s;
                        sm = sum(~isnan(s));
                        od1(sm+1) = branch_idx;
                        od2 = od1;
                        fod1_last = fod_last;
                        fod2_last = fod_last;
                        fod1_last(sm+1) = f(i);
                        fod2_last = fod1_last;
                        Fod_last = [Fod_last fod1_last fod2_last];
                        d1 = d;
                        d2 = d;
                        for ii = 1:N
                            if ii == branch_idx
                                s1(ii) = 1;
                            end
                        end
                        d1(branch_idx)=-1; %direction of bound is GE
                        for ii = 1:N
                            if ii == branch_idx
                                s2(ii) = 0;
                            end
                        end
                        d2(branch_idx)=1; %direction of bound is LE
                        nsold = size(S,2);

                        % add subproblems to the problem tree
                        S = [S s1 s2];
                        Od = [Od od1 od2];
                        D = [D d1 d2];
                        fs = [fs f(i) f(i)];  
                        nsnew = nsold+2;
                    else
                        nsnew = size(S,2);
                    end
           
                    if Length_WarmStart == 0
                        if o.branchMethod==1 || (o.branchMethod==13 && numIntSol<6)
                            %depth first, add newly branched problems to the
                            %beginning of the queue
                            priority = [nsold+1:nsnew priority];
                        elseif o.branchMethod==2
                            %breadth first, add newly branched problems to the
                            %end of the queue
                            priority = [priority nsold+1:nsnew];
                        elseif o.branchMethod==3 || (o.branchMethod==13 && numIntSol>=6)
                            %branch on the best lp solution
                            priority = [nsold+1:nsnew priority];
                            [dum,pidx] = sort(fs(priority));
                            priority=priority(pidx);
                        elseif o.branchMethod==4
                            %branch on the worst lp solution
                            priority = [nsold+1:nsnew priority];
                            [dum,pidx] = sort(-fs(priority));
                            priority=priority(pidx);
                        else
                            error('BBSQP: Unknown branch method.')
                        end
                    else
                        priority = [priority nsold+1:nsnew];
                        Length_WarmStart = Length_WarmStart -1;
                    end
                end
            end
            if (strcmp(o.display,'improve') || strcmp(o.display,'iter')) && (f(i)==f_best || i==1) 
                disp(['It. ',num2str(i),'. Best integer solution: ',num2str(f_best),' Delta ',num2str(max([100*(f_best-min(fs(priority)))/abs(f_best+min(fs(priority))) 0 100*isinf(f_best)])),'%']);
                %disp(['Queue: ', num2str(length(priority))]);
            end
        else
            error('BBSQP: Problem neither infeasible nor solved, try another solver or reformulate problem!')
        end
        if strcmp(o.display,'iter')
            disp(['It. ',num2str(i),'. F-val(It): ',num2str(f(i)),' Delta ',num2str(max([100*(f_best-min(fs(priority)))/abs(f_best+min(fs(priority))) 0 100*isinf(f_best)])),'%. Queue len. ', num2str(length(priority))]);
        end
        if o.iterplot
            plot(i,[f_best min(fs(priority))],'x')
            if i==1
                legend('MIP','LP')
            end
            drawnow;
        end
    else %parent node is more expensive than current f-best -> don't evaluate this node
        priority(1) = [];
    end
end

if ~strcmp(o.display,'off')
    disp(['Iteration ', num2str(sum(~isnan(visited))), '. Optimization ended.']);
    if isempty(priority) || f_best<min(fs(priority))
        disp('Found optimal solution!')
    elseif ~isinf(f_best)
        disp(['Ended optimization. Current delta ',num2str(max([100*(f_best-min(fs(priority)))/f_best 0 100*isinf(f_best)])),'%']);
    else
        disp(['Did not find any integer solutions']);
    end
    disp(['Time spent ', num2str(toc), ' seconds']);
    disp(['Objective function value: ',num2str(f_best)]);
end