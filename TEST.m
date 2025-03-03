clear all;
close all;
clc;

rng('default');%Random seed

%% original image
x_bar = double( imread('bird.tif') );%Load the image


% blur operator
psf = fspecial('average', 3);%Type and size of kernel for the blur
%psf = fspecial('motion');


%parameters
reg1 = 1;
reg2 = 5;

% data fidelity
K = @(x) imfilter(x, psf,'symmetric');
KT = @(x) imfilter(x, rot90(psf,2),'symmetric'); %K transpose
mu = reg1*sum(abs(psf(:)));
mu=1/mu;%cocoercivity constant

%% forward finite differences (with Neumann boundary conditions)
L1= @(x) [x(:,2:end,:)-x(:,1:end-1,:), zeros(size(x,1),1,size(x,3))]; % horizontal finite diference
L2 = @(x) [x(2:end,:,:)-x(1:end-1,:,:); zeros(1,size(x,2),size(x,3))]; % vertical finite diference

% backward finite differences (with Neumann boundary conditions)
L1T = @(x) [-x(:,1,:), x(:,1:end-2,:)-x(:,2:end-1,:), x(:,end-1,:)];    % horizontal finite diference transpose
L2T = @(x) [-x(1,:,:); x(1:end-2,:,:)-x(2:end-1,:,:); x(end-1,:,:)];    % vertical finite diference transpose


zeta = sqrt(8);%Lipschitz constant
tol = 1e-6;%tolerance 
maxiter = 10000;%maximum number of iteration


%stepsizes
kappa = 0.5;%Tunning parameter
gam = kappa*2/(4*zeta+1/mu);%step-size

for i=1:20
    'kappa = 0.5' 
    b =  imfilter(x_bar, psf) + 10 * randn( size(x_bar) );%Blurred and noisy immage.

    %%%%%%%%%%%%%%%%%%%%%%%% Inertial FHRB %%%%%%%%%%%%%%%%%%%%%%%
    A = 1-zeta*gam-gam/(2*mu);
    a1 = 0.9999999*(2*A+1-sqrt((2*A+1)^2-4*(A-1)*(A-zeta*gam)))/(2*(A-1));%maximum alpha
    a1 = a1*[1/3 2/3 1];%alpha range
    (1-a1).^2.*A-zeta*gam-a1.*(1+a1);
    for j=1:3
        [x_FRHFRI1{i,j},t_FRHFRI1(i,j),fo_FRHFRI1(i,j),error_FRHFRI1(i,j),iter_FRHFRI1(i,j)] = FRHF_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a1(j),1);
    end%run for three different alpha the algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% Double inertial FHRB %%%%%%%%%%%%%%%%%%%%%%%
    beta = 1;%beta
    alpha = 0.9999*(-(3-2*gam*zeta)+sqrt((3-2*gam*zeta)^2+4*gam*zeta*(1-2*gam*zeta-(1-beta).^2*gam/(2*mu))))/(2*gam*zeta);%maximum alpha
    alpha = alpha*[1/3 1/2 1];%alpha range
    theta = max(0.99999*(1-3*alpha -gam*(1-beta).^2/(2*mu)-gam*zeta*(1+(1-alpha).^2))/3,0.999*(gam*beta/(2*mu)-alpha+zeta*gam*alpha));%theta as large as possible
    %alpha+theta-gam*beta/(2*mu)-zeta*gam*alpha;
    %1-3*(alpha+theta)-gam*(1-beta).^2/(2*mu)-gam*zeta-zeta*gam*(1-alpha).^2;

    for j=1:3
        [x_FRHFDI{i,j},t_FRHFDI(i,j),fo_FRHFDI(i,j),error_FRHFDI(i,j),iter_FRHFDI(i,j)] = FRHF_double_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,alpha(j),beta,theta(j));
    end%run for three different alpha the algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%% SDI  FHRB %%%%%%%%%%%%%%%%%%%%%%%
    beta2 = 1;%beta=1 and alpha=0
    alpha2 = 0;
    theta2 = max(0.99999*(1-3*alpha2 -gam*(1-beta2).^2/(2*mu)-gam*zeta*(1+(1-alpha2).^2))/3,0.999*(gam*beta2/(2*mu)-alpha2+zeta*gam*alpha2));%theta as large as possible
    %alpha2+theta2-gam*beta2/(2*mu)-zeta*gam*alpha2;
    %1-3*(alpha2+theta2)-gam*(1-beta2).^2/(2*mu)-gam*zeta-zeta*gam*(1-alpha2).^2;

    [x_FRHFSDI{i,1},t_FRHFSDI(i,1),fo_FRHFSDI(i,1),error_FRHFSDI(i,1),iter_FRHFSDI(i,1)] = FRHF_double_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,alpha2,beta2,theta2);%run the algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%% FRHF %%%%%%%%%%%%%%%%%%%%%%%
    [x_FRHF{i,1},t_FRHF(i,1),fo_FRHF(i,1),error_FRHF(i,1),iter_FRHF(i,1)] = FRHF_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,0,1);%run the vanilla algortihm alpha=beta=theta=0 and lamba=1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%% Semi Inertial FHRB %%%%%%%%%%%%%%%%%%%%%%%
    a4 = (1-gam/(2*mu)-2*gam*zeta)/3;%maximum theta
    a4 = a4*[1/3 2/3 1];%range for theta
    for j=1:3
        [x_FRHFSI{i,j},t_FRHFSI(i,j),fo_FRHFSI(i,j),error_FRHFSI(i,j),iter_FRHFSI(i,j)] = FRHF_Semi_inertial(b,reg1,reg2,K,KT,L1,L2,L1T,L2T,mu,zeta,tol,maxiter,gam,a4(j));
    end%run for three different theta the algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




