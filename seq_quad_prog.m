%  Function: [x,fval,it,X] = seq_quad_prog (f, gradf, hessf,
%                                    A, b, G, r, x0, itmax, tol)
%
% min.  1/2 z'Hz + f'z
% s.t.  Gz = h
%       Az â‰? b
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

function [u_sqp,xtra,fval,it,ocp_res,time] = seq_quad_prog(xhat,km,opt,N,size_x,size_u,s,d)
% %% compute weights of the terminal cost
% q = 10*diag([10,10,10,1,1,1]);
% r  = diag(0.1*ones(3,1));
% [K,pf,~] = dlqr(DfDx(xx(:,1)',uu(:,1)'),DfDu(xx(:,1)',uu(:,1)'),q,r);
stop = false;
it = 1;
uu = reshape(km.z(1:size_u*N,1),[size_u,N]);
xx = reshape(km.z(size_u*N+1:size_u*N+size_x*N,1),[size_x,N]);
s = km.z(size_u*N+size_x*N + 1:size_u*N+size_x*N + N, 1);
load size_Cbar;
load size_f;
tic;
    while( ~stop )
        uu = reshape(km.z(1:size_u*N,1),[size_u,N]);
        xx = reshape(km.z(size_u*N+1:size_u*N+size_x*N,1),[size_x,N]);
        s = km.z(size_u*N+size_x*N + 1:size_u*N+size_x*N + N, 1);
       %% compute gradients of the constraints 
        DuDp = [];
        DxDc = [];
        for j = 1:N
            der.Cu(:,:,j) = DCDu(xx(:,j)',uu(:,j)');
            DuDp = blkdiag(DuDp,der.Cu(:,:,j));
            der.Cbarx(:,:,j) = DCbarDx(xx(:,j)');
            DxDc = blkdiag(DxDc,der.Cbarx(:,:,j));
        end
        Gamma = kron(speye(N), ones(size_Cbar,1));
        DzDh = blkdiag(DuDp, DxDc);
        DzDh = [DzDh,[zeros(size(DuDp,1),N);-Gamma];zeros(N,size(DzDh,2)),-eye(N)];
        % compute gradients of the dynamic equation, referring to the page 38 of
        % modele6 slides
        DuDg = [];
        for j = 1:N
            if j == 1
                der.fu(:,:,j) = DfDu(xhat',uu(:,j)');
                DuDg = blkdiag(DuDg,-der.fu(:,:,j));
            else
                der.fu(:,:,j) = DfDu(xx(:,j-1)',uu(:,j)');
                DuDg = blkdiag(DuDg,-der.fu(:,:,j));
            end
        end
        for j = 1:N-1
            der.fx(:,:,j) = DfDx(xx(:,j)',uu(:,j+1)');                
        end
        DxDg = [eye(size_x),zeros(size_x),zeros(size_x);-der.fx(:,:,1),eye(size_x),zeros(size_x);zeros(size_x),-der.fx(:,:,2),eye(size_x)];
        DzDg = [DuDg,DxDg,zeros(N*size_f,N)];
        % compute gradients of the cost function, referring to the page 37 of
        % modele6 slides
        DxDJ = [];
        DuDJ = [];
        for j = 1:N
            if j  < N
                der.Lx(:,:,j) = DLDx(xx(:,j)',uu(:,j+1)');
                DxDJ = [DxDJ,der.Lx(:,:,j)];
            else
                der.Phix(:,:,j) = DPhiDx(xx(:,j)');
                DxDJ = [DxDJ,der.Phix(:,:,j)];
            end
            if j == 1
                der.Lu(:,:,j) = DLDu(xhat',uu(:,j)');
                DuDJ = [DuDJ,der.Lu(:,:,j)];
            else
                der.Lu(:,:,j) = DLDu(xx(:,j-1)',uu(:,j)');
                DuDJ = [DuDJ,der.Lu(:,:,j)];
            end
        end
        DzDJ = [DuDJ,DxDJ,opt.gamma*ones(1,N)];
        % compute gradient of the Lagrangian, referring to the page 41 of
        % modele6 slides
        DzDL = DzDJ' + DzDh'*km.v + DzDg'*km.l;
        % natural residual
        G = [];
        G = [xx(:,1)- valf(xhat', uu(:,1)')];
        for j=1:N-1                   % generate dynamic constraints
             G = [G;xx(:,j+1)- valf(xx(:,j)', uu(:,j+1)')];
        end
        H = [];
        for j = 1:N
            H = [H; valC(xx(:,j)',uu(:,j)') + equ_con];
        end
        for j = 1:N
            H = [H; valCbar(xx(:,j)')];
        end
        H = [H; -s];
        F_NR = [DzDL; G; min(-H,km.v)];  % |-H|?
        ocp_res = norm(F_NR);
        % compute Hessian of the cost function, referring to the page 37 of
        % modele6 slides
        DuuDL = [];
        DxxDL = [];
        DxuDL = [];
        for j = 1:N
            if j == 1
                der.Luu(:,:,j) = D2LDuu(xhat',uu(:,j)');
                DuuDL = blkdiag(DuuDL,der.Luu(:,:,j));
            else
                der.Luu(:,:,j) = D2LDuu(xx(:,j-1)',uu(:,j)');
                DuuDL = blkdiag(DuuDL,der.Luu(:,:,j));
            end
            if j == 1
                der.Lxu(:,:,j) = D2LDxu(xx(:,j)',uu(:,j)');
                DxuDL = blkdiag(DxuDL,der.Lxu(:,:,j));
            else
                der.Lxu(:,:,j) = D2LDxu(xx(:,j-1)',uu(:,j)');
                DxuDL = blkdiag(DxuDL,der.Lxu(:,:,j));
            end
            if j == 1
                der.Lxx(:,:,j) = D2LDxx(xhat',uu(:,j)');
                DxxDL = blkdiag(DxxDL,der.Lxx(:,:,j));
            elseif j < N
                der.Lxx(:,:,j) = D2LDxx(xhat',uu(:,j)');
                DxxDL = blkdiag(DxxDL,der.Lxx(:,:,j));
            else
%                 der.Phixx(:,:,j) = D2PhiDxx(xx(:,j)',reshape(pf,[1,size_x^2]));
                der.Phixx(:,:,j) = D2PhiDxx(xx(:,j)');
                DxxDL = blkdiag(DxxDL,der.Phixx(:,:,j));
            end
        end
        Hessian = [DuuDL,DxuDL,zeros(size(DuuDL,1),N); DxuDL',DxxDL,zeros(size(DxxDL,1),N); zeros(N,size(DuuDL,1)),zeros(N,size(DxxDL,1)),zeros(N,N)];
        qp.H = Hessian; 
        eigenvalues = eig(qp.H);
        if all(eigenvalues(:) >= 0) == false
            fprintf('The Hessian is not positive semidefinite, the sub-QP is non-convex QPs\n');
            break
        end            
        qp.f = DzDJ';        
        qp.G = DzDg;
        qp.h = -G;  
        qp.A = full(DzDh);
        qp.b = -H;
        opts.display_level = 2;
        opts.max_iters = 5;
        opts.tol = 1e-4;
        [x,out] = fbstab_sparse(km,qp,opts);
        if( norm(F_NR) < opt.tol )
			stop = true; % => x is the solution
        else
            km.z = km.z + x.z;
            km.l = x.l;
            km.v = x.v;
			it = it + 1;		
%             clear der.Cu der.Cbarx
		end
        % If there are too many iterations
		if (it >= opt.itmax)
			stop = true;
        end
    end  
    time = toc;
    u_sqp = km.z(1:size_u*N);
    fval =0;
    for j = 1:N
        if j == 1
            fval = fval + valL(xhat',uu(:,j)');
        elseif j == N
%             fval = fval + valPhi(xx(:,N)',reshape(pf,[1,size_x^2]));
            fval = fval + valPhi(xx(:,N)');
        else
            fval = fval + valL(xx(:,j-1)',uu(:,j)');
        end
    end
    xtra = km;
end