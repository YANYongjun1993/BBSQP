%  Function: [x,fval,it,X] = seq_quad_prog (f, gradf, hessf,
%                                    A, b, G, r, x0, itmax, tol)
%
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az <= b
%
%  using the sequential quadratic programming method.
%
%  f     : The name of the objective function f : R^n -> R.
%  gradf : The name of the function that returns the gradient of f.
%  hessf : The name of the function that returns the hessian matrix of f.
%  A     : A matrix with dimension m x n for the linear equality contraints.
%  b     : A vector with m elements for the linear equality contraints.
%  G     : A matrix with dimension p x n for the linear inequality constraints.
%  r     : A vector with p elements for the linear inequality constraints.
%  x0    : The start point for the algorithm.
%  itmax : The maximal number of iterations allowed.
%  tol   : The bound needed for the stop criteria.

function [x0,fval,flag] = seq_quad_prog(x0,v0,opt,s,d,yidx)
% %% compute weights of the terminal cost
% q = 10*diag([10,10,10,1,1,1]);
% r  = diag(0.1*ones(3,1));
% [K,pf,~] = dlqr(DfDx(xx(:,1)',uu(:,1)'),DfDu(xx(:,1)',uu(:,1)'),q,r);

%Calculate the extra constraints imposed by s
%if s is all NaN then no new constraints
sidx = ~isnan(s);
ydiag = double(yidx);
if any(sidx)
    ydiag(yidx) = d;
    l = length(yidx);  
    Y = spdiags(ydiag,0,l,l);
    %Remove zero entries
    Y = Y(any(Y,2),:);
end
stop = false;
it = 1;
flag = 1;
% uu = reshape(km.z(1:size_u*N,1),[size_u,N]);
% xx = reshape(km.z(size_u*N+1:size_u*N+size_x*N,1),[size_x,N]);
% s = km.z(size_u*N+size_x*N + 1:size_u*N+size_x*N + N, 1);
tic;
    while( ~stop )
       %% compute gradients of the constraints 
        DzDh = DCDz(x0);
        if any(sidx)
            DzDh = [DzDh;Y];
        end
        if any(sidx) && it == 1
            v0 = [v0;zeros(size(Y,1),1)];
        end
        % compute gradients of the dynamic equation, referring to the page 38 of
        % modele6 slides
        DzDg = [];
        % compute gradients of the cost function, referring to the page 37 of
        % modele6 slides
        DzDJ = DLDz(x0);
        % compute gradient of the Lagrangian, referring to the page 41 of
        % modele6 slides
        DzDL = DzDJ' + DzDh'*v0;
        % natural residual
        H = [];
        H =  valC(x0);
        sidx_1 = [double(sidx);0;0;0];
        sidx_1 = logical(sidx_1)';
        if any(sidx)
            H = [H; -( s(sidx).*d(sidx)) + d(sidx).*x0(sidx_1)' ]; % The extra constraints imposed by s
        end
        F_NR = [DzDL; min(-H,v0)];  % |-H|?
        ocp_res = norm(F_NR);
        % compute Hessian of the cost function, referring to the page 37 of
        % modele6 slides
        DzzDL = [];
        DzzDL = D2LDzz(x0);
        Hessian = [DzzDL];
        %% solve using quadprog
%         qp.H = Hessian; 
%         qp.f = DzDJ';        
%         qp.G = DzDg;
%         qp.h = [];  
%         qp.A = full(DzDh);
%         qp.b = -H;
%         options = optimoptions('quadprog');
%         options = optimoptions(options,'StepTolerance',1e-8);
%         [x,~,~,~,duals] = quadprog(qp.H,qp.f,qp.A,qp.b,[],[],[],[],[],options);
        %% solve using fastab_dense
        qp.H = Hessian; 
        eigenvalues = eig(qp.H);
        if all(eigenvalues(:) >= 0) == false
            fprintf('The Hessian is not positive semidefinite, the sub-QP is non-convex QPs\n');
            break
        end  
        qp.f = DzDJ';  
        qp.A = full(DzDh);
        qp.b = -H;
        opts.display_level = 2;
        opts.max_iters = 5;
        opts.tol = 1e-4;
        [x,v,out] = fbstab_dense(qp,x0',v0,opts);
        if( norm(F_NR) < opt.tol )
			stop = true; % => x is the solution
        else
            x0 = x0 + x';
            v0 = v;
% 			  km.z = km.z + x';
%             km.l = duals.eqlin;
%             km.v = duals.ineqlin;
			it = it + 1;
%             clear der.Cu der.Cbarx
		end
        % If there are too many iterations
		if (it >= opt.itmax)
            flag = 0;
			stop = true;
        end
    end  
    time = toc;
    fval =  valL(x0);
end