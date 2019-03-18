% Approximating smooth funcitons (without noise)
% See Experiment 1

clear all, close all
N = 600; L = 200; 
ff = -1:0.01:1; 
lambda = 10^(-1); mu = ones(L+1,1);
Lmax = 150; LL = 1:Lmax;
% you can try different values of parameters N, L, lambda, mu and Lmax, but
% you should remember the Gauss quadrature rule which requires L <= 2N+1.

% basis index, include 1 (1st Chebyshev), 2 (Legendre). note that Legendre
% case requires much more run time than Chebyshev case due to different
% segments (chebpoly and legpoly) in Chebfun 

basis_idx = 1;

switch basis_idx
    case 1
        [f,w] = jacpts(N+1,-.5,-.5); 
    case 2
        [f,w] = legpts(N+1);
end

% example index, include 1, 2, 3, 4, 5, 6, 7, 8
% you can add another index here and add config in next part by yourself
example_idx = 1;

%% generate config 
switch example_idx
    case 1
        Y = tanh(20*sin(12*f)) + .02*exp(3*f).*sin(300*f); 
        YY = tanh(20*sin(12*ff)) + .02*exp(3*ff).*sin(300*ff); 
    case 2
        Y = tanh(4*f-1); YY = tanh(4*ff-1);
    case 3
        Y = abs(f); YY = abs(ff);
    case 4
        Y = 10*airy(30*f); YY =  10*airy(30*ff);
    case 5 
        Y = sin(6*f)+sign(sin(f+exp(2*f))); YY = sin(6*ff)+sign(sin(ff+exp(2*ff)));
    case 6 
        Y = 1./(1+25*f.^2); YY = 1./(1+25*ff.^2); % Runge's function
    case 7
        Y = f./(exp(f)-1); YY = ff./(exp(ff)-1);
    case 8
        Y = exp(-f.^2); YY = exp(-ff.^2);
end

%% Approximation
switch basis_idx
    case 1
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
    case 2
        for l = 0:L
            F = legpoly(l);
            A(:,l+1) = F(f)/sqrt(2/(2*l+1));
        end
end
beta2 = l2_beta(w,A,Y,lambda,L,mu);
beta1 = l1_beta(w,A,Y,lambda,L,mu);
beta0 = l2_beta(w,A,Y,0,L,mu);
p2 = zeros(201,1);
p1 = zeros(201,1);
p0 = zeros(201,1);
for l = 0:L
    switch basis_idx
        case 1
            if l == 0
                T = cos(l*acos(ff))/sqrt(pi);
            else
            T = cos(l*acos(ff))/sqrt(pi/2);
            end
            T = T';
        case 2
            F = legpoly(l)/sqrt(2/(2*l+1));
            T = F(ff');
    end
    p2 = p2+beta2(l+1)*T;
    p1 = p1+beta1(l+1)*T;
    p0 = p0+beta0(l+1)*T;
end

%% Plots 1,2
fontsize_baseline = 10;
subplot(2,2,1), plot(ff,p0,'-','linewidth',1), hold on, plot(ff,p2,'--','linewidth',1),plot(ff,p1,'-.','linewidth',1),plot(ff,YY,':','linewidth',1.3),...
      grid on, box on, xlabel('$x$', 'interpreter','latex','fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),....
      title(['Approximation, $N=$' num2str(N) ', $L=$' num2str(L) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...
      set(legend({'no regularization','$\ell_2-\ell_2$','$\ell_2-\ell_1$','original function'}),'interpreter','latex', 'fontsize', fontsize_baseline)
      axis([-1,1,-1.5,2.5]), set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
 subplot(2,2,2), plot(ff,abs(YY-p0'),'-','linewidth',1), hold on,plot(ff,abs(YY-p2'),'--','linewidth',1.2), plot(ff,abs(YY-p1'),'-.','linewidth',1),...
      grid on, box on, xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
      title(['Error, $N=$' num2str(N) ', $L=$' num2str(L) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...
      set(legend({'no regularization','$\ell_2-\ell_2$','$\ell_2-\ell_1$'}),'interpreter','latex', 'fontsize', fontsize_baseline)
      axis([-1,1,0,2]), set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),

%% Cases involving L = 1:Lmax      
for vL = 1:Lmax % variable L
    fprintf(['Current L (max L = ',num2str(Lmax), '): ',num2str(vL) '\n'])
    mu = ones(vL+1,1);
    A = [];
    switch basis_idx
        case 1
            for l = 0:vL
                for j = 0:N
                    A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);      
                end
            end
            A(:,1) = A(:,1)/sqrt(2);
        case 2
            for l = 0:vL
                F = legpoly(l);
                A(:,l+1) = F(f)/sqrt(2/(2*l+1));
            end
    end
    beta0 = l2_beta(w,A,Y,0,vL,mu);
    beta2 = l2_beta(w,A,Y,lambda,vL,mu);
    beta1 = l1_beta(w,A,Y,lambda,vL,mu); 
    p0 = zeros(201,1); p2 = zeros(201,1); p1 = zeros(201,1);
    P0 = zeros(N+1,1); P2 = zeros(N+1,1); P1 = zeros(N+1,1);
    for l = 0:vL
        switch basis_idx
            case 1
                if l == 0
                    T = cos(l*acos(ff))/sqrt(pi);
                    TT = cos(l*acos(f))/sqrt(pi);
                else
                    T = cos(l*acos(ff))/sqrt(pi/2);
                    TT = cos(l*acos(f))/sqrt(pi/2);
                end
                T = T';
            case 2
            F = legpoly(l);
            T = F(ff)/sqrt(2/(2*l+1)); T = T';
            TT = F(f)/sqrt(2/(2*l+1)); 
        end
    p0 = p0+beta0(l+1)*T; p2 = p2+beta2(l+1)*T; p1 = p1+beta1(l+1)*T;
    P0 = P0+beta0(l+1)*TT; P2 = P2+beta2(l+1)*TT; P1 = P1+beta1(l+1)*TT;
    end
    error0(vL) = norm(YY-p0',inf); error2(vL) = norm(YY-p2',inf); error1(vL) = norm(YY-p1',inf); % unifom error
    ERROR0(vL) = sqrt(w*(abs(Y-P0)).^2); ERROR2(vL) = sqrt(w*(abs(Y-P2)).^2); ERROR1(vL) = sqrt(w*(abs(Y-P1)).^2); % L_2 norm
end

%% Plots 3,4
subplot(2,2,3),plot(LL,error0,'-','linewidth',1),hold on,plot(LL,error2,'--','linewidth',1),hold on, plot(LL,error1,'-.','linewidth',1),...
    xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Uniform error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title(['Uniform errors with $N=$' num2str(N) ', $1\leq L\leq$' num2str(Lmax) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...    
    set(legend({'no regularization','uniform error for $\ell_{2}-\ell_{2}$','uniform error for $\ell_{2}-\ell_{1}$'}),'interpreter','latex', 'fontsize', fontsize_baseline),...
    grid on, box on, set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
subplot(2,2,4),plot(LL,ERROR0,'-','linewidth',1),hold on,plot(LL,ERROR2,'--','linewidth',1),hold on, plot(LL,ERROR1,'-.','linewidth',1),...
    xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline), ylabel('$L_{2}$ error','interpreter','latex', 'fontsize', fontsize_baseline),set(legend({'no regularization','$L_2$ error for $\ell_2-\ell_2$','$L_2$ error for $\ell_2-\ell_1$'}), 'interpreter','latex','fontsize', fontsize_baseline),
    title(['$L_{2}$ errors with $N=$' num2str(N) ', $1\leq L\leq$' num2str(Lmax) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...    
    grid on, box on, set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')