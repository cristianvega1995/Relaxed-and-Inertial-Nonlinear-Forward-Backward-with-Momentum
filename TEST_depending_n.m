clear all;
close all;
clc;

rng('default');%Random seed

%% original image
x_bar = double( imread('bird.tif') );%load the image


% blur operator
psf = fspecial('average', 3);%type and size of blur
%psf = fspecial('motion');


%parameters
reg1 = 1;
reg2 = 5;

% data fidelity
K = @(x) imfilter(x, psf,'symmetric');
KT = @(x) imfilter(x, rot90(psf,2),'symmetric');%K transpose 
mu = reg1*sum(abs(psf(:)));
mu=1/mu;%cocoercivity parameter

%% forward finite differences (with Neumann boundary conditions)
L1= @(x) [x(:,2:end,:)-x(:,1:end-1,:), zeros(size(x,1),1,size(x,3))]; % horizontal finite diference
L2 = @(x) [x(2:end,:,:)-x(1:end-1,:,:); zeros(1,size(x,2),size(x,3))]; % vertical finite diference

% backward finite differences (with Neumann boundary conditions)
L1T = @(x) [-x(:,1,:), x(:,1:end-2,:)-x(:,2:end-1,:), x(:,end-1,:)];    % horizontal finite diference transpose
L2T = @(x) [-x(1,:,:); x(1:end-2,:,:)-x(2:end-1,:,:); x(end-1,:,:)];    % vertical finite diference transpose


zeta = sqrt(8);%Lipschitz constant
tol = 1e-6;%tolerance
maxiter = 10000;%maximum number of iteration

% noisy image
load '20images';%load 20 blurred and noisy images


%stepsizes
kappa1 = 0.8;%Tunning parameter
kappa2 = 0.99;%Tunning parameter

for i=1:20
    b = BB{i};%select the i-th image
    [x_FRHFDI{i,1},t_FRHFDI(i,1),fo_FRHFDI(i,1),error_FRHFDI(i,1),iter_FRHFDI(i,1)] = FRHF_inertial_depending_on_n(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,kappa1,kappa2,1);%Experiment 1 for varying parameters
    [x_FRHFDI{i,2},t_FRHFDI(i,2),fo_FRHFDI(i,2),error_FRHFDI(i,2),iter_FRHFDI(i,2)] = FRHF_inertial_depending_on_n2(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,kappa1,kappa2,1);%Experiment 2 for varying parameters
    [x_FRHFDI{i,3},t_FRHFDI(i,3),fo_FRHFDI(i,3),error_FRHFDI(i,3),iter_FRHFDI(i,3)] = FRHF_inertial_depending_on_n3(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,kappa1,kappa2,1);%Experiment 3 for varying parameters
    
end

% for  iter = 1:100;
%  ka= kappa1 + kappa3/iter;
%  gam(iter) = ka*2/(4*zeta+1/mu);
%  A = 1-zeta*gam(iter)-gam(iter)/(2*mu);
%  a(iter) = 0.9999999*(2*A+1-sqrt((2*A+1)^2-4*(A-1)*(A-zeta*gam(iter))))/(2*(A-1));
% end

save('results_09_tol6_22_dependingn')%stored the results


%mean(iter_FRHFDI)

%mean(t_FRHFDI)

% system('shutdown -s')

