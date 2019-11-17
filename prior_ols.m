function [Beta,VBOLS,sigma_u] = prior_ols(data,p)

[T,n]=size(data);

%OLS estimation
x=[];
for i= 1:p
for j = 1:n
    x = [x data(p+1-i:end-i,j)];
end
end
x=[ones(T-p,1) x];

y= data(1+p:end,:);
Beta= inv(x'*x)*x'*y;
u= y-x*Beta;
sigma_u =(u'*u)/T;
VBOLS=kron(sigma_u,inv(x'*x));

end