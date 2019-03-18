clear all
close all
N =500; L= N;
[x,w,v] = jacpts(N+1,-.5,-.5);
xx = linspace(-1, 1, 1000);
lambda = 10^(-.5); mu = ones(L+1,1); 
for l = 0:L
    for j = 0:N
        A(j+1,l+1) = cos(l*acos(x(j+1)))/sqrt(pi/2);
    end
end
A(:,1) = A(:,1)/sqrt(2);
for k = 1:4
    example_idx = k;
    switch example_idx
        case 1
            G = 1./( 1 + 25*x.^2 );
            GG = 1./( 1 + 25*xx.^2 );
        case 2
            G = airy(20*x);  
            GG = airy(20*xx); 
        case 3
            G = exp(-x.^2);
            GG = exp(-xx.^2);
        case 4
            G = tanh(20*sin(12*x)) + .02*exp(3*x).*sin(300*x); 
            GG = tanh(20*sin(12*xx)) + .02*exp(3*xx).*sin(300*xx); 
    end
   [Y,NOISE] = noisegen(G,10);
% Y = G;
subplot(4,4,k), plot(x,Y), title('original function')
beta1 = l1_beta(w,A,Y,lambda,L,mu);
p1 = zeros(1000,1);
for l = 0:L        
            if l == 0
                T = cos(l*acos(xx))/sqrt(pi);
            else
            T = cos(l*acos(xx))/sqrt(pi/2);
            end
            T = T';
    p1 = p1+beta1(l+1)*T;
end
ff1 = p1;
subplot(4,4,k+4), plot(xx,ff1),title('$p_{1}$: closed-form solution','interpreter','latex')
a = zeros(L+1,1);
Ax = A'*(w'.*Y);
for i = 0:L
    if 2*Ax(i+1)>lambda*mu(i+1)
        a(i+1) = -lambda*mu(i+1)/2;
    elseif 2*Ax(i+1)<-lambda*mu(i+1)
        a(i+1) = lambda*mu(i+1)/2;
    else
        a(i+1) = -Ax(i+1);
    end
end
AA = a'*A';
ff2 = bary(xx, Y+AA',x, v);
subplot(4,4,k+8),plot(xx,ff2), title('$p_{2}$: barycentric interpolant','interpreter','latex')
subplot(4,4,k+12), plot(xx,abs(ff1-ff2')), title('$|p_{1}-p_{2}|$','interpreter','latex')
end