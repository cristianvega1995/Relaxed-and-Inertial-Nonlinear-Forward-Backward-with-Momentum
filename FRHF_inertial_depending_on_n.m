function [x,time,fo,error,iter] = FRHF_inertial_depending_on_n(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,kappa1,kappa2,lam)
% FRHF inertial varying parameter: Inertial forward reflected half forward splitting for image restoration
%Semi inertial beta>0 other parameters equal to zero alpha>=0theta>=0
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
% kappa1,kappa2: tunning parameter
% lam: relaxation parameter

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

C = @(X)  lam1g*(KTK(X)-Kb);%Operator C cocoercive
%Initial error
error=1;
% Initialization of variables
error = 1;
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
p1 = x1;
p21 = x21;
p22 = x22;


kappa3 = kappa2-kappa1;
%iteration counter
iter=0;

tic
while error > tol & iter<maxiter
     iter = iter +1;%update iteration
     ka= kappa1 + kappa3/iter;%update  paramerters
     gam = ka*2/(4*zeta+1/mu);%update  paramerters
     A = 1-zeta*gam-gam/(2*mu);%update  paramerters
     a = 0.9999999*(2*A+1-sqrt((2*A+1)^2-4*(A-1)*(A-zeta*gam)))/(2*(A-1));%update  paramerters
     %(1-a).^2.*(1-zeta*gam-gam/(2*mu))-zeta*gam-(1+a).*a
    %Parameters
    lams=lam2/gam;
    a1=a+2;
    a2=1+2*a;
    a3=1+a;
    a3g= a3/gam;
    ag = a/gam;
    
    a1g = a1*gam;
    a2g = a2*gam;
    apg = a*gam;
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


     y1 = (1+a)*x1_ - a*xo1_;%Inertial step

    xx1 = (2+a)*x2_1 - (1+2*a)*xo2_1 + a*xoo2_1;
    xx2 = (2+a)*x2_2 - (1+2*a)*xo2_2 + a*xoo2_2;
    q1 = y1 - gam*(L1T(xx1)+L2T(xx2)+C(y1));

    qq = (2+a)*x1_ -(1+2*a)*xo1_ + a*xoo1_;
    qq1 = L1(qq);
    qq2 = L2(qq);
    q21 = ((1+a)*x2_1 - a*xo2_1)/gam + qq1;%Inertial step
    q22 = ((1+a)*x2_2 - a*xo2_2)/gam + qq2;%Inertial step

    x1 = max(min(q1,255),0);%Projection



    xx2 = abs(q21) -  lams; 
    x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %Prox conjugate l1 norm

    xx2 = abs(q22) -  lams;      
    x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %Prox conjugate l1 norm
    % % 
    % 
    %  xoo1_ = xo1_;
    %      xo1_ = x1_;
    %      x1_ = x1;
    % 
    %      xoo2_1 = xo2_1;
    %      xoo2_2 = xo2_2;
    % 
    %      xo2_1 = x2_1;
    %      xo2_2 = x2_2;
    % 
    %      x2_1 = x21;
    %      x2_2 = x22;
    % 
    % 
    %      y1 = a3*x1_-a*xo1_;
    % 
    %     xx1 = a1g*x2_1 -a2g*xo2_1 + apg*xoo2_1;
    %     xx2 = a1g*x2_2 -a2g*xo2_2 + apg*xoo2_2;
    %     q1 = y1 - (L1T(xx1)+L2T(xx2)+gam*C(y1));
    % 
    %     qq=a1*x1_ -a2*xo1_ + a*xoo1_;
    %     qq1 = L1(qq);
    %     qq2 = L2(qq);
    %     q21 = a3g*x2_1 - ag*xo2_1 + qq1;
    %     q22 = a3g*x2_2 - ag*xo2_2 + qq2;
    % 
    %     x1 = max(min(q1,255),0);
    % 
    % 
    % 
    %     xx2 = abs(q21) -  lams; 
    %     x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1
    % 
    %     xx2 = abs(q22) -  lams;      
    %     x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1
    % --- Compute relative primal dual error 
    error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));

end
% --- Post-processing ---
time = toc;% Record total execution time
x = x1;% Final restored image% Record total execution time
fo = lam1*norm(K(x)-b,2)^2/2+lam2*sum(sum(abs(x)));% Compute final objective function value (data term + L1 regularization)


