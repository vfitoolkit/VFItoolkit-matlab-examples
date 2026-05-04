function F=Aiyagari1994_EndoLabor_EpsteinZin_ReturnFn(l, aprime, a, z,chi,r,w)
% Action space: l,aprime,a,z
% Rest are parameters

F=-Inf; % placeholder

c=w*l*z+(1+r)*a-aprime; % budget constraint

if c>0 && l<1
    F=((c^(1-chi)) *((1-l)^chi));
end

end