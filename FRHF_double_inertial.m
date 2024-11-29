function [x,time,fo,error,iter] = FRHF_double_inertial(b,lam1,lam2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a,beta,t)



%Operator C cocoercive
lam1g = lam1*gam;
KTK = @(X) KT(K(X));
Kb = KT(b);
C = @(X)  lam1g*(KTK(X)-Kb);

error=1;
x1 = b;
x1_ = x1;
x21 = L1(b);
x22 = L2(b);
x2_1 = x21;
x2_2 = x22;
xo1_ = x1;
xo2_1 = x21;
xo2_2 = x22;
z1 = 2*x1_ - xo1_;
z21 = 2*x2_1 - xo2_1; 
z22 = 2*x2_2 - xo2_2;    
lams=lam2/gam;

iter=0;

a1=a+2;
a2=1+2*a;
at= a+t;
atg = at/gam;
at2 = 1+at;
at2g = at2/gam;
b1= beta+1;

a1g = a1*gam;
a2g = a2*gam;
apg = a*gam;

if beta ~= 1
   % 'beta ~= 1'

    tic
    while error > tol & iter<maxiter
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
    
    
         z1 = b1*x1_ - beta*xo1_;
    
       
        xx1 = a1g*x2_1 -a2g*xo2_1 + apg*xoo2_1;
        xx2 = a1g*x2_2 -a2g*xo2_2 + apg*xoo2_2;
        q1 = at2*x1_ - at*xo1_- (L1T(xx1)+L2T(xx2)+C(z1));
    
        qq = a1*x1_ -a2*xo1_ + a*xoo1_;
        %qq1 = L1(qq);
        %qq2 = L2(qq);
        q21 = at2g*x2_1 - atg*xo2_1 + L1(qq);
        q22 = at2g*x2_2 - atg*xo2_2 + L2(qq);
    
        x1 = max(min(q1,255),0);
    
        
    
        xx2 = abs(q21) -  lams; 
        x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1
    
        xx2 = abs(q22) -  lams;      
        x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1
    
    
        error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    end
    time = toc;
elseif beta == 1 && a ~= 0
    %'beta = 1 y alpha ~= 0'
    C = @(X)  lam1*(KTK(X)-Kb);
    tic
    while error > tol & iter<maxiter
         iter = iter +1;
         
         xo1_ = x1_;
         x1_ = x1;
         
         xo2_1 = x2_1;
         xo2_2 = x2_2;
    
         x2_1 = x21;
         x2_2 = x22;

         z1_ = z1;
         z21_ = z21;
         z22_ = z22;
    
         z1  =  2*x1_ - xo1_;
         z21 = 2*x2_1 - xo2_1; 
         z22 = 2*x2_2 - xo2_2;    
       
        xx1 = z21 + a*(x2_1-z21_);
        xx2 = z22 + a*(x2_2-z22_);
        q1 = at2*x1_ - at*xo1_- gam*(L1T(xx1)+L2T(xx2)+C(z1));
    
        qq = z1 + a*(x1_- z1_);
        %qq1 = L1(qq);
        %qq2 = L2(qq);
        q21 = at2g*x2_1 - atg*xo2_1 + L1(qq);
        q22 = at2g*x2_2 - atg*xo2_2 + L2(qq);
    
        x1 = max(min(q1,255),0);
    
        
    
        xx2 = abs(q21) -  lams; 
        x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1
    
        xx2 = abs(q22) -  lams;      
        x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1
    
    
        error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    end
    time = toc;
elseif beta == 1 && a == 0
       %'beta = 1 y alpha == 0'
       C = @(X)  lam1*(KTK(X)-Kb);
    tic
    while error > tol & iter<maxiter
         iter = iter +1;
         
         xo1_ = x1_;
         x1_ = x1;
         
         xo2_1 = x2_1;
         xo2_2 = x2_2;
    
         x2_1 = x21;
         x2_2 = x22;

         z1_ = z1;
         z21_ = z21;
         z22_ = z22;
    
         z1 = 2*x1_ - xo1_;
         z21 = 2*x2_1 - xo2_1; 
         z22 = 2*x2_2 - xo2_2;    
       
        q1 = at2*x1_ - at*xo1_- gam*(L1T(z21)+L2T(z22)+C(z1));
    
        
        %qq1 = ;
        %qq2 = ;
        q21 = at2g*x2_1 - atg*xo2_1 + L1(z1);
        q22 = at2g*x2_2 - atg*xo2_2 + L2(z1);
    
        x1 = max(min(q1,255),0);
    
        
    
        xx2 = abs(q21) -  lams; 
        x21 = gam*(q21-sign(q21).*((xx2>0).*xx2)); %norma 1
    
        xx2 = abs(q22) -  lams;      
        x22 = gam*(q22-sign(q22).*((xx2>0).*xx2)); %norma 1
    
    
        error = sqrt((norm(x1(:)-x1_(:))^2+norm(x21(:)-x2_1(:))^2+norm(x22(:)-x2_2(:))^2)/(norm(x1_(:))^2+norm(x2_1(:))^2+norm(x2_2(:))^2));
    end
    time = toc;
end



x = x1;
fo = lam1*norm(K(x)-b,2)^2/2+lam2*sum(sum(abs(x)));

