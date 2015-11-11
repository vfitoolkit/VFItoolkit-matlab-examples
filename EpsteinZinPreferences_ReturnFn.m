function F=CaldaraFVillaverdeRRamirezYao2012_FmatrixFn(l_val,kprime_val,k_val,z_val, zeta, delta, upsilon)

c=exp(z_val)*(k_val^zeta)*(l_val^(1-zeta))+(1-delta)*k_val-kprime_val;

F=-Inf;
if c>0
    F=(c^upsilon)*(1-l_val).^(1-upsilon);
end

end