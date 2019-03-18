clear all
close all 
result = [];
N = 500;
maxL = 60;
ff = -1:0.01:1;
[f,w] = jacpts(N+1,-.5,-.5);
lambda = 1e-1;
fontsize_baseline = 10;
% colormap autumn
color = [0,128,255;
127,0,255;
255,0,0]/255;
setno = 1e-12;


for k = 1:4
    result = [];
    example_idx = k;
    switch example_idx
        case 1
            Y = 10*airy(10*f); 
        case 2 
            Y = sin(pi*f)+cos(pi*f);
        case 3 
            Y = randn(N+1,1);
        case 4
            Y = sin(pi*f)./(pi*f);
            Y(N/2+1) = 1;
    end
    for L = 1:maxL
        A = [];
        mu = ones(1,L+1);
        for l = 0:L
            for j = 0:N
                A(j+1,l+1) = cos(l*acos(f(j+1)))/sqrt(pi/2);
            end
        end
        A(:,1) = A(:,1)/sqrt(2);
        beta = l1_beta(w,A,Y,lambda,L,mu);
        AWY = A'*diag(w)*Y;
        estimateupper = 0;
        for i = 1:L+1
            if abs(AWY(i))>setno
                estimateupper = estimateupper+1;
            end
        end
        count = 0;
        for i = 1:L+1
            if abs(beta(i))==0 & abs(AWY(i))>setno
            count = count + 1;   
            end
        end
        number = 0;
        for i = 1:L+1
            if abs(beta(i))>setno
            number = number+1;
            end
        end
        b = [L estimateupper-count number estimateupper]; 
        result = [result;b];
        
    end
    LL = 1:maxL; U = result(maxL,4);
    subplot(2,2,k), plot(LL,result(:,4),'bp'), hold on, plot(LL,result(:,3),'k<'), plot(LL,result(:,2),'r*'), ...
        xlabel('Degree $L$','interpreter','latex', 'fontsize', fontsize_baseline), ylabel('Number','interpreter','latex', 'fontsize', fontsize_baseline),set(legend({'$\|{\rm\bf{A}}_L^T{\bf{\Lambda}}{\rm{\bf{f}}}\|_0$','$\|\beta\|_0$','$\|{\rm\bf{A}}_L^T{\bf{\Lambda}}{\rm{\bf{f}}}\|_0-\#$'}),'interpreter','latex','location','NorthWest', 'fontsize', fontsize_baseline),...
        set(gca, 'fontsize', fontsize_baseline), grid on, box on, axis([1,maxL,0,U*1.5]),...
        set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off')
    switch example_idx
        case 1
            title('Airy function $f(x) = 10{\rm Airy}(10x)$','interpreter','latex', 'fontsize', fontsize_baseline),
        case 2
            title('$f(x) = \sin(\pi x)+\cos(\pi x)$','interpreter','latex', 'fontsize', fontsize_baseline),
        case 3
            title('Random sampling $\verb"Y = randn(N+1,1)"$','interpreter','latex', 'fontsize', fontsize_baseline),
        case 4
            title('Sinc function $f(x) = \sin(\pi x))/(\pi x)$','interpreter','latex', 'fontsize', fontsize_baseline),
    end            
end
            
        