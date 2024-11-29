clear all;
close all;
clc;

rng('default');

%% original image
x_bar = double( imread('bird.tif') );


% blur operator
psf = fspecial('average', 9);
%psf = fspecial('motion');


%parameters
reg1 = 1;
reg2 = 5;

% data fidelity
K = @(x) imfilter(x, psf,'symmetric');
KT = @(x) imfilter(x, rot90(psf,2),'symmetric'); 
mu = reg1*sum(abs(psf(:)));
mu=1/mu;

%% forward finite differences (with Neumann boundary conditions)
L1= @(x) [x(:,2:end,:)-x(:,1:end-1,:), zeros(size(x,1),1,size(x,3))]; % horizontal
L2 = @(x) [x(2:end,:,:)-x(1:end-1,:,:); zeros(1,size(x,2),size(x,3))]; % vertical

% backward finite differences (with Neumann boundary conditions)
L1T = @(x) [-x(:,1,:), x(:,1:end-2,:)-x(:,2:end-1,:), x(:,end-1,:)];    % horizontal
L2T = @(x) [-x(1,:,:); x(1:end-2,:,:)-x(2:end-1,:,:); x(end-1,:,:)];    % vertical


zeta = sqrt(8);
tol = 1e-6;
maxiter = 10000;





%stepsizes
kappa = 0.99;
gam = kappa*2/(4*zeta+1/mu);
rng('default');

for i=1:20
    'kappa = 0.99' 
    i
    b =  imfilter(x_bar, psf) + 10 * randn( size(x_bar) );

    %%%%%%%%%%%%%%%%%%%%%%%% Inertial FHRB %%%%%%%%%%%%%%%%%%%%%%%
    % A = 1-zeta*gam-gam/(2*mu);
    % a1 = 0.9999999*(2*A+1-sqrt((2*A+1)^2-4*(A-1)*(A-zeta*gam)))/(2*(A-1));
    a1 = 0.2;
    re = [1000 2000 3000];
    %(1-a1).^2.*A-zeta*gam-a1.*(1+a1)
    for j=1:3
        j
        [x_FRHFRI1{i,j},t_FRHFRI1(i,j),fo_FRHFRI1(i,j),error_FRHFRI1(i,j),iter_FRHFRI1(i,j)] = FRHF_inertial_restart(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a1,re(j));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% FRHF %%%%%%%%%%%%%%%%%%%%%%%
    [x_FRHF{i,1},t_FRHF(i,1),fo_FRHF(i,1),error_FRHF(i,1),iter_FRHF(i,1)] = FRHF_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,0,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end







