function F=Aiyagari1994_ReturnFn(aprime, a, z,alpha,delta,gamma,r)
% Action space: aprime,a,z
% Rest are parameters

F=-Inf; % placeholder

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*z+(1+r)*a-aprime; 

if c>0
    if gamma==1
        F=log(c);
    else
        F=(c^(1-gamma) -1)/(1-gamma);
    end
end

end