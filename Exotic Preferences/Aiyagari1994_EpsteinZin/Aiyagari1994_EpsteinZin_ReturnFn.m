function F=Aiyagari1994_EpsteinZin_ReturnFn(aprime, a, z,alpha,delta,r)
% Action space: aprime,a,z
% Rest are parameters

F=-Inf; % placeholder

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

c=w*z+(1+r)*a-aprime; 

if c>0
    F=c; % As Epstein-Zin
end

end