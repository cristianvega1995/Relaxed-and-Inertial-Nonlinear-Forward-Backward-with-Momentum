function [x,time,fo,error,iter] = FRHF_inertial_restart(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a,re)



if a == 0
   error('The parameters are incorrect')
end


%Operator C cocoercive
lam1g= lam1*gam;
KTK = @(X) KT(K(X));
Kb = KT(b);
C = @(X)  lam1g*(KTK(X)-Kb);

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
lams=lam2/gam;

iter=0;
a1=a+2;
a2=1+2*a;
a3=1+a;
a3g= a3/gam;
ag = a/gam;

a1g = a1*gam;
a2g = a2*gam;
apg = a*gam;

tic
    while error > tol & iter<re
     iter = iter +1;
     
     xoo1_ = xo1_;
     xo1_ = x1_;
     x1_ = x1;
     
     xoo2_1 = xo2_1;
     xoo2_2 = xo2_2;

     xo2_1 = x2_1;
     xo2_2 = x2_2;

     x2_1 = x21;
     x2_2 = x22;


     y1 = a3*x1_-a*xo1_;
    
    xx1 = a1g*x2_1 -a2g*xo2_1 + apg*xoo2_1;
    xx2 = a1g*x2_2 -a2g*xo2_2 + apg*xoo2_2;
    q1 = y1 - (L1T(xx1)+L2T(xx2)+C(y1));

    qq=a1*x1_ -a2*xo1_ + a*xoo1_;
    qq1 = L1(qq);
    qq2 = L2(qq);
    q21 = a3g*x2_1 - ag*xo2_1 + qq1;
    q22 = a3g*x2_2 - ag*xo2_2 + qq2;

    x1 = max(min(q1,255),0);

    

    xx2 = abs(q21) -  lams; 
    x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1

    xx2 = abs(q22) -  lams;      
    x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1


    error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    
end
while error > tol & iter<maxiter
     iter = iter +1;
     
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

    x1 = max(min(q1,255),0);

    

    xx2 = abs(q21) -  lams; 
    x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1

    xx2 = abs(q22) -  lams;      
    x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1


    error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
 end
 time = toc;
 x = x1;
 fo = lam1*norm(K(x)-b,2)^2/2+lam2*sum(sum(abs(x)));

