function F=EpsteinZinPreferences_ReturnFn(l,kprime,k,z, zeta, delta, upsilon)

c=exp(z)*(k^zeta)*(l^(1-zeta))+(1-delta)*k-kprime;

F=-Inf;
if c>0
    F=(c^upsilon)*(1-l).^(1-upsilon);
end

end