function [x,time,fo,error,iter] = FRHF_Semi_inertial(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a)



%Operator C cocoercive
KTK = @(X) KT(K(X));
Kb = KT(b);
C = @(X)  lam1*(KTK(X)-Kb);

error = 1;
x1 = b;
x1_ = x1;
x21 = L1(b);
x22 = L2(b);
x2_1 = x21;
x2_2 = x22;
y1 = x1;
y21 = x21;
y22 = x22;
lams=lam2/gam;

iter=0;

tic
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
        q1 = x1_ +a*(x1_-xo1_)  - gam*(L1T(xx1)+L2T(xx2)+C(x1_));
    
        qq = 2*x1_-xo1_;
        %qq1 = L1(qq);
        %qq2 = L2(qq);
        q21 = (x2_1 + a*(x2_1-xo2_1))/gam+ L1(qq);
        q22 = (x2_2 + a*(x2_2-xo2_2))/gam+ L2(qq);
    
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