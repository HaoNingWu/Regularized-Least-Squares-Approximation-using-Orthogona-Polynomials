function beta = l2_beta(w,A,f,lambda,L,mu)
 Af = A'*diag(w)*f;
% beta = conjgrad(ALA,Af);
%  for i = 1:L+1
%       beta(i) = Af(i)/(1+lambda*mu(i)^2);
%  end
beta = Af./(1+lambda*mu.^2);
%beta = (diag(ones([L+1,1]))+lambda*diag(ones([L+1,1])))\(A'*diag(w)*f);
%beta = (A'*diag(w)*A+lambda*diag(ones([L+1,1])))\(A'*diag(w)*f);
%beta = (A'*diag(w)*A+lambda*A'*diag(w)*A)\(A'*diag(w)*f);
% beta = (A'*diag(w)*A+lambda*ones(L+1))\(A'*diag(w)*f);
end
