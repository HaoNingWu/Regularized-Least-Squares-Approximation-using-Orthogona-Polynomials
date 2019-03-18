% This script is used to illustrate some theoritical results in
% approximating smooth and noisy functions.
% Congpei An and Hao-Ning Wu, Rugularized Least Squares Approximation by
% Orthogonal Poynomials, arXiv:1805.01140
% Modied by Hao-Ning Wu 2019-02
% Contact information:
% Hao-Ning Wu
% Department of Mathematics, Jinan University, Guangzhou, China 
% hwu@stu2015.jnu.edu.cn, haoning.wu@outlook.com 

%% 
% clear variables and windows
clear all, close all

example_idx = 3;

% The index of numerical example, example_idx, include 
% [1] smoothapproximation.m£ºapproximate a smooth function (This function 
% stems from L. N. Trefethen's Approximation Theory and Approximation 
% Practice)----example_idx = 1; 
% [2] reduction.m: recovery from a noisy function (spectral density of gate 
% signal)----example_idx = 2; 
% [3] highoscillateddenoising.m: recovery from a highly oscillatory 
% function with noise (function y=airy(40x))----example_idx = 3; 

switch example_idx
    case 1 %% Numerical Experiment 1       
        N = 600; L = 200; 
        ff = -1:0.01:1; 
        lambda = 10^(-1); mu = ones(L+1,1);
        Lmax = 150; LL = 1:Lmax;
        [f,w] = jacpts(N+1,-.5,-.5); 
        Y = tanh(20*sin(12*f)) + .02*exp(3*f).*sin(300*f); 
        YY = tanh(20*sin(12*ff)) + .02*exp(3*ff).*sin(300*ff); 
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
        beta2 = l2_beta(w,A,Y,lambda,L,mu);
        beta1 = l1_beta(w,A,Y,lambda,L,mu);
        beta0 = l2_beta(w,A,Y,0,L,mu);
        p2 = zeros(201,1);
        p1 = zeros(201,1);
        p0 = zeros(201,1);
        for l = 0:L
            if l == 0
                T = cos(l*acos(ff))/sqrt(pi);
            else
                T = cos(l*acos(ff))/sqrt(pi/2);
            end
            T = T';
            p2 = p2+beta2(l+1)*T;
            p1 = p1+beta1(l+1)*T;
            p0 = p0+beta0(l+1)*T;
        end
        fontsize_baseline = 10;
        subplot(2,2,1), plot(ff,p0,'-','linewidth',1), hold on, plot(ff,p2,'--','linewidth',1),plot(ff,p1,'-.','linewidth',1),plot(ff,YY,':','linewidth',1.3),...
            grid on, box on, xlabel('$x$', 'interpreter','latex','fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),....
            title(['Approximation, $N=$' num2str(N) ', $L=$' num2str(L) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...
            set(legend({'no regularization','$\ell_2-\ell_2$','$\ell_2-\ell_1$','original function'}),'interpreter','latex', 'fontsize', fontsize_baseline)
            axis([-1,1,-1.5,2.5]), set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        subplot(2,2,2), plot(ff,abs(YY-p0'),'-','linewidth',1), hold on,plot(ff,abs(YY-p2'),'--','linewidth',1.2), plot(ff,abs(YY-p1'),'-.','linewidth',1),...
            grid on, box on, xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title(['Error, $N=$' num2str(N) ', $L=$' num2str(L) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...
            set(legend({'no regularization','$\ell_2-\ell_2$','$\ell_2-\ell_1$'}),'interpreter','latex', 'fontsize', fontsize_baseline)
            axis([-1,1,0,2]), set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),   
        for vL = 1:Lmax % variable L
            fprintf(['Current L (max L = ',num2str(Lmax), '): ',num2str(vL) '\n'])
            mu = ones(vL+1,1);
            A = [];
            for l = 0:vL
                for j = 0:N
                    A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);      
                end
            end
            A(:,1) = A(:,1)/sqrt(2);
            beta0 = l2_beta(w,A,Y,0,vL,mu);
            beta2 = l2_beta(w,A,Y,lambda,vL,mu);
            beta1 = l1_beta(w,A,Y,lambda,vL,mu); 
            p0 = zeros(201,1); p2 = zeros(201,1); p1 = zeros(201,1);
            P0 = zeros(N+1,1); P2 = zeros(N+1,1); P1 = zeros(N+1,1);
            for l = 0:vL
                if l == 0
                    T = cos(l*acos(ff))/sqrt(pi);
                    TT = cos(l*acos(f))/sqrt(pi);
                else
                    T = cos(l*acos(ff))/sqrt(pi/2);
                    TT = cos(l*acos(f))/sqrt(pi/2);
                end
                T = T';
                p0 = p0+beta0(l+1)*T; p2 = p2+beta2(l+1)*T; p1 = p1+beta1(l+1)*T;
                P0 = P0+beta0(l+1)*TT; P2 = P2+beta2(l+1)*TT; P1 = P1+beta1(l+1)*TT;
            end
            error0(vL) = norm(YY-p0',inf); error2(vL) = norm(YY-p2',inf); error1(vL) = norm(YY-p1',inf); % unifom error
            ERROR0(vL) = sqrt(w*(abs(Y-P0)).^2); ERROR2(vL) = sqrt(w*(abs(Y-P2)).^2); ERROR1(vL) = sqrt(w*(abs(Y-P1)).^2); % L_2 norm
        end
        subplot(2,2,3),plot(LL,error0,'-','linewidth',1),hold on,plot(LL,error2,'--','linewidth',1),hold on, plot(LL,error1,'-.','linewidth',1),...
            xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Abosolute uniform error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title(['Uniform errors with $N=$' num2str(N) ', $1\leq L\leq$' num2str(Lmax) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...    
            set(legend({'no regularization','uniform error for $\ell_{2}-\ell_{2}$','uniform error for $\ell_{2}-\ell_{1}$'}),'interpreter','latex', 'fontsize', fontsize_baseline),...
            grid on, box on, set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        subplot(2,2,4),plot(LL,ERROR0,'-','linewidth',1),hold on,plot(LL,ERROR2,'--','linewidth',1),hold on, plot(LL,ERROR1,'-.','linewidth',1),...
            xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline), ylabel('Absolute $L_{2}$ error','interpreter','latex', 'fontsize', fontsize_baseline),set(legend({'no regularization','$L_2$ error for $\ell_2-\ell_2$','$L_2$ error for $\ell_2-\ell_1$'}), 'interpreter','latex','fontsize', fontsize_baseline),
            title(['$L_{2}$ errors with $N=$' num2str(N) ', $1\leq L\leq$' num2str(Lmax) ', $\lambda=$' num2str(lambda)],'interpreter','latex', 'fontsize', fontsize_baseline),...    
            grid on, box on, set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
    case 2
        N = 100; L = 35; maxL = 100;    
        ff = -1:0.01:1; 
        Lambda = [0 10^(-1)];
        for l = 0:1:L
            mu(l+1) = 1/Filter(L,l);
        end
        mu = mu';
        [f,w] = jacpts(N+1,-.5,-.5); 
        GG = sin(pi*ff*5)./(pi*ff); GG(101) = 5;
        G = sin(pi*f*5)./(pi*f);
        if mod(N+1,2) == 1
            G((N+2)/2) = 5;
        end
        [YY,NOISE] = noisegen(GG,10); [Y,NOISE] = noisegen(G,10);
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
        beta2 = l2_beta(w,A,Y,Lambda(2),L,mu);
        beta1 = l1_beta(w,A,Y,Lambda(2),L,mu);
        beta0 = l1_beta(w,A,Y,Lambda(1),L,mu);
        p2 = zeros(201,1); p1 = zeros(201,1); p0 = zeros(201,1);    
        for l = 0:L        
            if l == 0
                F = chebpoly(l,1)/sqrt(pi);
            else
                F = chebpoly(l,1)/sqrt(pi/2);
            end
        p2 = p2+beta2(l+1)*F(ff');
        p1 = p1+beta1(l+1)*F(ff');
        p0 = p0+beta0(l+1)*F(ff');
        end
        Color = [215,25,28; 253,174,97; 254,204,92;
            171,217,233; 44,123,182]/255;
        fontsize_baseline = 10;
        figure(1)
        subplot(2,4,1), plot(ff,GG,'linewidth',1.2,'color',Color(1,:)), box on,...
            xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('Original function','interpreter','latex', 'fontsize', fontsize_baseline),...
            set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
        subplot(2,4,5), plot(ff,YY,'linewidth',1.2,'color',Color(2,:)),box on,...
            xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('Noisy function','interpreter','latex', 'fontsize', fontsize_baseline),...
            set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
        subplot(2,4,2),plot(ff,p0,'linewidth',1.5,'color',Color(3,:)), box on,...
            xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('$\lambda=0$','interpreter','latex', 'fontsize', fontsize_baseline), set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
        subplot(2,4,6), hold on, plot(ff,abs(GG-p0'),'-','color',Color(3,:),'linewidth',1.5),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
            title('Error, $\lambda=0$','interpreter','latex'),box on, set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
        subplot(2,4,3),plot(ff,p2,'linewidth',1.5,'color',Color(4,:)),box on,...
            xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),title('$\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),...
            set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
        subplot(2,4,7), hold on, plot(ff,abs(GG-p2'),'-','color',Color(4,:),'linewidth',1.5),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error, $\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),box on,axis([-1,1,0,1.5]), set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        subplot(2,4,4),plot(ff,p1,'linewidth',1.5,'color',Color(5,:)), box on,...
            xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baseline),title('$\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),...
            set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,-2,6]),
        subplot(2,4,8), hold on, plot(ff,abs(GG-p1'),'-','color',Color(5,:),'linewidth',1),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error: $\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),box on,axis([-1,1,0,1.5]), set(gca, 'fontsize', fontsize_baseline), grid on,...
            set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        for vL = 1:maxL % variable L
            mu = []; A = [];
            for l = 0:1:vL
                mu(l+1) = 1/Filter(vL,l);
            end
            mu = mu';
            for l = 0:vL
                for j = 0:N
                    A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);      
                end
            end
            A(:,1) = A(:,1)/sqrt(2);
            beta0 = l2_beta(w,A,Y,Lambda(1),vL,mu); 
            beta2 = l2_beta(w,A,Y,Lambda(2),vL,mu); 
            beta1 = l1_beta(w,A,Y,Lambda(2),vL,mu); 
            p0 = zeros(201,1); P0 = zeros(N+1,1); p2 = zeros(201,1); 
            P2 = zeros(N+1,1); p1 = zeros(201,1); P1 = zeros(N+1,1);
            for l = 0:vL
                if l == 0
                    T = cos(l*acos(ff))/sqrt(pi);
                    TT = cos(l*acos(f))/sqrt(pi);
                else
                    T = cos(l*acos(ff))/sqrt(pi/2);
                    TT = cos(l*acos(f))/sqrt(pi/2);
                end
                T = T';
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
    case 3
        N = 500; L = 300; LL = 100;
        ff = -1:0.01:1; [f,w] = jacpts(N+1,-.5,-.5);
        GG =  airy(40*ff); G= airy(40*f);
        [YY,NOISE] = noisegen(GG,10); [Y,NOISE] = noisegen(G,10);
        Lambda = [0 10^(-1) 10^(-2) 10^(-5)];
        for l = 0:1:L
            mu(l+1) = 1/Filter(L,l);
        end
        mu = mu';
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
        Color = [215,25,28; 253,174,97; 254,204,92;
            171,217,233; 44,123,182]/255;
        fontsize_baseline = 10;
        subplot(4,4,1), plot(ff,GG), xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Original function','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,2), plot(ff,YY), xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Noisy function','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'),

        beta2 = l2_beta(w,A,Y,Lambda(1),L,mu);
        p2 = zeros(201,1);
        for l = 0:L        
            if l == 0
                F = chebpoly(l,1)/sqrt(pi);
            else
                F = chebpoly(l,1)/sqrt(pi/2);
            end
            p2 = p2+beta2(l+1)*F(ff');
        end
        subplot(4,4,5),plot(ff,p2),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('$\lambda=0$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,6), plot(ff,abs(GG-p2')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
            title('Error, $\lambda=0$','interpreter','latex'),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off'),
 
        beta21 = l2_beta(w,A,Y,Lambda(2),L,mu); beta11 = l1_beta(w,A,Y,Lambda(2),L,mu);
        beta22 = l2_beta(w,A,Y,Lambda(3),L,mu); beta12 = l1_beta(w,A,Y,Lambda(3),L,mu);
        beta23 = l2_beta(w,A,Y,Lambda(4),L,mu); beta13 = l1_beta(w,A,Y,Lambda(4),L,mu);
        p21 = zeros(201,1); p11 = zeros(201,1);
        p22 = zeros(201,1); p12 = zeros(201,1);
        p23 = zeros(201,1); p13 = zeros(201,1);
        for l = 0:L
            if l == 0
                F = chebpoly(l,1)/sqrt(pi);
            else
                F = chebpoly(l,1)/sqrt(pi/2);
            end
            p21 = p21+beta21(l+1)*F(ff');
            p11 = p11+beta11(l+1)*F(ff');
            p22 = p22+beta22(l+1)*F(ff');
            p12 = p12+beta12(l+1)*F(ff');
            p23 = p23+beta23(l+1)*F(ff');
            p13 = p13+beta13(l+1)*F(ff');
        end

        subplot(4,4,9),plot(ff,p21),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,10),plot(ff,abs(GG-p21')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error, $\ell_2-\ell_2$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        subplot(4,4,13),plot(ff,p11),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,14),plot(ff,abs(GG-p11')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error: $\ell_2-\ell_1$ model, $\lambda=10^{-1}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
        subplot(4,4,3),plot(ff,p22),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_2$ model, $\lambda=10^{-2}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,4),plot(ff,abs(GG-p22')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error, $\ell_2-\ell_2$ model, $\lambda=10^{-2}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off')
        subplot(4,4,7),plot(ff,p12),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_1$ model, $\lambda=10^{-2}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,8), plot(ff,abs(GG-p12')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error: $\ell_2-\ell_1$ model, $\lambda=10^{-2}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]), set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),    
        subplot(4,4,11),plot(ff,p23),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_2$ model, $\lambda=10^{-5}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline), set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,12), plot(ff,abs(GG-p23')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error, $\ell_2-\ell_2$ model, $\lambda=10^{-5}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off')
        subplot(4,4,15),plot(ff,p13),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$f(x)$', 'interpreter','latex','fontsize', fontsize_baseline),...
            title('$\ell_2-\ell_1$ model, $\lambda=10^{-5}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,-.6,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
        subplot(4,4,16), plot(ff,abs(GG-p13')),xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
            title('Error: $\ell_2-\ell_1$ model, $\lambda=10^{-5}$','interpreter','latex', 'fontsize', fontsize_baseline),grid on,box on,...
            axis([-1,1,0,.6]),set(gca, 'fontsize', fontsize_baseline),set(gca, 'XMinorGrid', 'off'),set(gca, 'YMinorGrid', 'off')
end
