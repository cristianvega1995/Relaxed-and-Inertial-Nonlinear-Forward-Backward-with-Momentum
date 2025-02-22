function [x,time,fo,error,iter] = FRHF_inertial_restart(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a,re)
% FRHF_double_inertial: double-inertial forward reflected half forward splitting for image restoration
%Semi inertial beta>0 other parameters equal to zero alpha=theeta=0
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
% re: iteration when we start to use the parameter alpha and beta

% Outputs:
% x: restored image after the algorithm has converged
% time: execution time taken for the algorithm to run
% fo: objective function value
% error: error between iterations for stopping condition
% iter: total number of iterations



%Avoid the case when alpha=0
if a == 0
   error('The parameters are incorrect')
end



lam1g= lam1*gam;% lambda multiplied by gamma
% Precompute frequently used operators
% KTK computes K^T * K
KTK = @(X) KT(K(X));
% Kb computes K^T * b
Kb = KT(b);
C = @(X)  lam1g*(KTK(X)-Kb);%Operator C cocoercive
%Initial error
error = 1;
% Initialization of variables
x1 = b;
x1_ = x1;
xo1_ = x1;
x21 = L1(b);
x22 = L2(b);
x2_1 = x21;
x2_2 = x22;
y1 = x1;
y21 = x21;
y22 = x22;
xo2_1 = x21;
xo2_2 = x22;
lams=lam2/gam;% lambda scaled by gamma
%iteration counter
iter=0;
%Parameters 
a1=a+2;
a2=1+2*a;
a3=1+a;
a3g= a3/gam;
ag = a/gam;

a1g = a1*gam;
a2g = a2*gam;
apg = a*gam;

tic % Start measuring time
 % Main loop: iterate until the error is below the tolerance or max iterations are reached
    while error > tol & iter<re%This is the first loop before restarting time
     iter = iter +1;%update iteration
     %Stored previous variable
     xoo1_ = xo1_;
     xo1_ = x1_;
     x1_ = x1;
     
     xoo2_1 = xo2_1;
     xoo2_2 = xo2_2;

     xo2_1 = x2_1;
     xo2_2 = x2_2;

     x2_1 = x21;
     x2_2 = x22;


     y1 = a3*x1_-a*xo1_;%Inertial step with alpha
    
    xx1 = a1g*x2_1 -a2g*xo2_1 + apg*xoo2_1;
    xx2 = a1g*x2_2 -a2g*xo2_2 + apg*xoo2_2;
    q1 = y1 - (L1T(xx1)+L2T(xx2)+C(y1));

    qq=a1*x1_ -a2*xo1_ + a*xoo1_;
    qq1 = L1(qq);
    qq2 = L2(qq);
    q21 = a3g*x2_1 - ag*xo2_1 + qq1;
    q22 = a3g*x2_2 - ag*xo2_2 + qq2;

    x1 = max(min(q1,255),0);%Projection

    

    xx2 = abs(q21) -  lams; 
    x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %Prox conjugate l1 norm

    xx2 = abs(q22) -  lams;      
    x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %Prox conjugate l1 norm

     % --- Compute relative primal dual error 
    error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    
end
while error > tol & iter<maxiter%This is the second  loop before restarting time
     iter = iter +1;%update iteration
     %Stored previous variable
     xo1_ = x1_;
     x1_ = x1;

     xo2_1 = x2_1;
     xo2_2 = x2_2;

     x2_1 = x21;
     x2_2 = x22;

    xx1 = 2*x2_1-xo2_1;
    xx2 = 2*x2_2-xo2_2;
    q1 = x1_ - gam*(L1T(xx1)+L2T(xx2))-C(x1_);

    qq = 2*x1_-xo1_;
    qq1 = L1(qq);
    qq2 = L2(qq);
    q21 = x2_1/gam + qq1;
    q22 = x2_2/gam + qq2;

    x1 = max(min(q1,255),0);%Projection

    

    xx2 = abs(q21) -  lams; 
    x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %Prox conjugate l1 norm

    xx2 = abs(q22) -  lams;      
    x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %Prox conjugate l1 norm

     % --- Compute relative primal dual error 
    error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    end5
% --- Post-processing ---
 time = toc;% Record total execution time
 x = x1;% Final restored image
 fo = lam1*norm(K(x)-b,2)^2/2+lam2*sum(sum(abs(x)));% Compute final objective function value (data term + L1 regularization


