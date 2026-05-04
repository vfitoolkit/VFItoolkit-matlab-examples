function F=Aiyagari1994_EndoLabor_ReturnFn(l, aprime, a, z,gamma_c,gamma_l,chi,r,w,nonSeperableUtility)
% Action space: l,aprime,a,z
% Rest are parameters

F=-Inf; % placeholder

c=w*l*z+(1+r)*a-aprime; % budget constraint

if c>0 && l<1
    if nonSeperableUtility==0
        if gamma_c==1 && gamma_l==1
            F=log(c)+chi*log(1-l);
        elseif gamma_c==1
            F=log(c)+chi*((1-l)^(1-gamma_l) -1)/(1-gamma_l);
        elseif gamma_l==1
            F=(c^(1-gamma_c) -1)/(1-gamma_c)+chi*log(1-l);
        else
            F=(c^(1-gamma_c) -1)/(1-gamma_c) +chi*((1-l)^(1-gamma_l) -1)/(1-gamma_l);
        end
    else % nonSeperableUtility==1
        % gamma_c plays the role of gamma
        F=(((c^(1-chi))*((1-l)^chi))^(1-gamma_c))/(1-gamma_c);
    end
end

end