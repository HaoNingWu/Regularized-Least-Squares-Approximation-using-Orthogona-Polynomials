% Denosing
% See Experiment 3
clear all, close all
N = 100; L = 25; maxL = 100; 
ff = -1:0.01:1; 
Lambda = [0 10^(-1)];
for l = 0:1:L
    mu(l+1) = 1/Filter(L,l);
end
mu = mu';

% you can try different values of parameters N, L, lambda, mu and maxL, but
% you should always keep in mind that the Gauss quadrature rule requires 
% L <= 2N+1.

% output index, include 1, 2. if you require the result under the certain L
% and N, take 1. you will see growing error with L ranging from 1 to maxL 
% under case 2 

output_idx = 2;

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


% example index, include 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13
% case 2 denotes a highoscillated function, to obtain a detailed output
% figure, you may turn to highoscillateddenoising.m
% you can add another index here and add config in next part by yourself
example_idx = 1;

%% generate config 
switch example_idx
    case 1
        GG = sin(pi*ff*5)./(pi*ff); GG(101) = 5;
        G = sin(pi*f*5)./(pi*f);
        if mod(N+1,2) == 1
            G((N+2)/2) = 5;
        end
    case 2
        GG =  airy(40*ff); G= airy(40*f);    
    case 3    
        G = sin(40*f)+cos(20*f); GG = sin(40*ff)+cos(20*ff);
    case 4
        G = sin(6*f)+sign(sin(f+exp(2*f))); GG = sin(6*ff)+sign(sin(ff+exp(2*ff)));
    case 5
        G = sin(6*f)+sin(60*exp(f)); GG = sin(6*ff)+sin(60*exp(ff));
    case 6
        G = 1./(1+1000*(f+.5).^2)+1./sqrt(1+1000*(f-.5).^2); GG = 1./(1+1000*(ff+.5).^2)+1./sqrt(1+1000*(ff-.5).^2);
    case 7
        G = sign(f)-f./2; GG = sign(ff)-ff./2;
    case 8
        G = 1./(1+25*f.^2); GG = 1./(1+25*ff.^2);
    case 9
        G = airy(f); GG = airy(ff); 
    case 10
        G = exp(-f.^2); GG = exp(-ff.^2);
    case 11
        G = sin(pi*f)./(pi*f); GG = sin(pi*ff)./(pi*ff);
        if mod(N+1,2) == 1
            G((N+2)/2) = 5;
        end
        GG(101) = 5;
    case 12
        G = sin(pi*f*5)./(pi*f)+abs(pi*f); GG = sin(pi*ff*5)./(pi*ff)+abs(pi*ff);
        if mod(N+1,2) == 1
            G((N+2)/2) = 5;
        end
        GG(101) = 5;
    case 13
        G = tanh(20*sin(12*f)) + .02*exp(3*f).*sin(300*f); GG = tanh(20*sin(12*ff)) + .02*exp(3*ff).*sin(300*ff);
end
%%
% add noise
[YY,NOISE] = noisegen(GG,10); [Y,NOISE] = noisegen(G,10);

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
beta2 = l2_beta(w,A,Y,Lambda(2),L,mu);
beta1 = l1_beta(w,A,Y,Lambda(2),L,mu);
beta0 = l1_beta(w,A,Y,Lambda(1),L,mu);
p2 = zeros(201,1); p1 = zeros(201,1); p0 = zeros(201,1);
for l = 0:L        
    switch basis_idx
        case 1
            if l == 0
                F = chebpoly(l,1)/sqrt(pi);
            else
                F = chebpoly(l,1)/sqrt(pi/2);
            end
        case 2
            F = legpoly(l)/sqrt(2/(2*l+1));
    end
    p2 = p2+beta2(l+1)*F(ff');
    p1 = p1+beta1(l+1)*F(ff');
    p0 = p0+beta0(l+1)*F(ff');
end

%% Figure 1
Color = [215,25,28;
253,174,97;
254,204,92;
171,217,233;
44,123,182]/255;
fontsize_baseline = 10;
figure(1)
subplot(2,4,1), plot(ff,GG,'linewidth',1,'color',Color(1,:)), box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('Original function','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
subplot(2,4,5), plot(ff,YY,'linewidth',1,'color',Color(2,:)),box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('Noisy function','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
subplot(2,4,2),plot(ff,p0,'linewidth',1,'color',Color(3,:)), box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('$\lambda=0$','interpreter','latex', 'fontsize', fontsize_baseline), set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
subplot(2,4,6), hold on, plot(ff,abs(GG-p0'),'-','color',Color(3,:),'linewidth',1),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('Error, $\lambda=0$','interpreter','latex'),box on, set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
subplot(2,4,3),plot(ff,p2,'linewidth',1,'color',Color(4,:)),box on,...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('$\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
subplot(2,4,7), hold on, plot(ff,abs(GG-p2'),'-','color',Color(4,:),'linewidth',1),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error, $\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),box on,axis([-1,1,0,1.5]), set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
subplot(2,4,4),plot(ff,p1,'linewidth',1,'color',Color(5,:)), box on,...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baseline),title('$\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
subplot(2,4,8), hold on, plot(ff,abs(GG-p1'),'-','color',Color(5,:),'linewidth',1),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('Error: $\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),box on,axis([-1,1,0,1.5]), set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')

 

 %% Errors involving variable L
switch output_idx
    case 1
        return
    case 2
for vL = 1:maxL % variable L
    mu = []; A = [];
    for l = 0:1:vL
        mu(l+1) = 1/Filter(vL,l);
    end
    mu = mu';
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
    beta0 = l2_beta(w,A,Y,Lambda(1),vL,mu); 
    beta2 = l2_beta(w,A,Y,Lambda(2),vL,mu); 
    beta1 = l1_beta(w,A,Y,Lambda(2),vL,mu); 
    p0 = zeros(201,1); P0 = zeros(N+1,1); p2 = zeros(201,1); P2 = zeros(N+1,1); p1 = zeros(201,1); P1 = zeros(N+1,1);
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
    p0 = p0+beta0(l+1)*T; P0 = P0+beta0(l+1)*TT;
    p2 = p2+beta2(l+1)*T;P2 = P2+beta2(l+1)*TT;  
    p1 = p1+beta1(l+1)*T; P1 = P1+beta1(l+1)*TT;
    error0(vL) = norm(GG-p0',inf); ERROR0(vL) = sqrt(w*(abs(G-P0)).^2);
    error2(vL) = norm(GG-p2',inf); ERROR2(vL) = sqrt(w*(abs(G-P2)).^2);
    error1(vL) = norm(GG-p1',inf); ERROR1(vL) = sqrt(w*(abs(G-P1)).^2);
    end
    fprintf(['Current L (max L = ',num2str(maxL), '): ',num2str(vL) '\n'])
end
LL1 = 1:maxL;
figure(2)
subplot(1,2,1), plot(LL1,error0,'linestyle','none','marker','h'),xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),...
     hold on, set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
subplot(1,2,2), plot(LL1,ERROR0,'linestyle','none','marker','h'),xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),...
     hold on, set(gca, 'fontsize', fontsize_baseline), grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
subplot(1,2,1), hold on, plot(LL1,error2,'linestyle','none','marker','^'),hold on, plot(LL1,error1,'linestyle','none','marker','s'),xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Error','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('Uniform errors with $1\leq L\leq 100$','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(legend({'uniform error for model with $\lambda=0$','uniform error for $\ell_2-\ell_2$ model with $\lambda=10^{-1}$','uniform error for $\ell_2-\ell_1$ model with $\lambda=10^{-1}$'}),'interpreter','latex', 'fontsize', fontsize_baseline)
subplot(1,2,2), hold on, plot(LL1,ERROR2,'linestyle','none','marker','^'),plot(LL1,ERROR1,'linestyle','none','marker','s'),xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Error','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('$L_2$ errors with $1\leq L\leq 100$','interpreter','latex', 'fontsize', fontsize_baseline),...
     set(legend({'$L_2$ error for model with $\lambda=0$','$L_2$ error for $\ell_2-\ell_2$ model with $\lambda=10^{-1}$','$L_2$ error for $\ell_2-\ell_1$ model with $\lambda=10^{-1}$'}),'interpreter','latex', 'fontsize', fontsize_baseline)
end