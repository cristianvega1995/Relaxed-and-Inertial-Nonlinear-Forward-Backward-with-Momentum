function [x,time,fo,error,iter] = FRHF_Semi_inertial(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a)
% FRHF_Semi_inertial: Semi-inertial forward reflected half forward splitting for image restoration
%Semi inertial theta>0 other parameters equal to zero alpha=beta=0
% Inputs:
% b: observed noisy and blurred image
% lam1, lam2: regularization parameters for the L1 and L2 norms
% K: forward operator representing the blur process
% KT: transpose of the operator K
% L1, L2: operators modeling horizontal and vertical differences
% L1T, L2T: transposes of L1 and L2
% mu, zeta: cocoercive constant and Lipschitz constant
% tol: tolerance for stopping criteria
% maxiter: maximum number of iterations
% gam: step size parameter
% a: inertial parameter

% Outputs:
% x: restored image after the algorithm has converged
% time: execution time taken for the algorithm to run
% fo: objective function value
% error: error between iterations for stopping condition
% iter: total number of iterations

% Precompute frequently used operators
% KTK computes K^T * K
KTK = @(X) KT(K(X));
% Kb computes K^T * b
Kb = KT(b);
%Operator C cocoercive
C = @(X)  lam1*(KTK(X)-Kb);
%Initial error
error = 1;
% Initialization of variables
x1 = b;
x1_ = x1;
x21 = L1(b);
x22 = L2(b);
x2_1 = x21;
x2_2 = x22;
y1 = x1;
y21 = x21;
y22 = x22;
lams=lam2/gam;% lambda scaled by gamma
%iteration counter
iter=0;
% Start measuring time
tic
% Main loop: iterate until the error is below the tolerance or max iterations are reached
while error > tol & iter<maxiter
         %update iteration
         iter = iter +1;
         %Stored previous variable
         xo1_ = x1_;
         x1_ = x1;
    
         xo2_1 = x2_1;
         xo2_2 = x2_2;
    
         x2_1 = x21;
         x2_2 = x22;
    
        xx1 = 2*x2_1-xo2_1;
        xx2 = 2*x2_2-xo2_2;
        q1 = x1_ +a*(x1_-xo1_)  - gam*(L1T(xx1)+L2T(xx2)+C(x1_));%Primal inertial step
    
        qq = 2*x1_-xo1_;
        %qq1 = L1(qq);
        %qq2 = L2(qq);
        q21 = (x2_1 + a*(x2_1-xo2_1))/gam+ L1(qq);%Dual inertial step
        q22 = (x2_2 + a*(x2_2-xo2_2))/gam+ L2(qq);%Dual inertial step
    
        x1 = max(min(q1,255),0); %Projection
    
        

        xx2 = abs(q21) -  lams; 
        x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %Prox conjugate l1 norm
 
        xx2 = abs(q22) -  lams;      
        x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %Prox conjugate l1 norm
    
        % --- Compute relative primal dual error 
        error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
  
end




% --- Post-processing ---

time = toc;% Record total execution time
x = x1;% Final restored image
% Compute final objective function value (data term + L1 regularization)
fo = lam1*norm(K(x)-b,2)^2/2+lam2*sum(sum(abs(x)));