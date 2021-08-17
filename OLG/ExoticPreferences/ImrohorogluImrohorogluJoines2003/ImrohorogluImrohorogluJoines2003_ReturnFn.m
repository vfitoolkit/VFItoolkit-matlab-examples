function F=ImrohorogluImrohorogluJoines2003_ReturnFn(l,kprime,bprime,k,ebar,b,z,r,tau_c, tau_k, tau_l,tau_u, tau_s,gamma, epsilon_j,alpha,delta, A,Tr_beq,g,agej,Jr,ebarspacing,ebendpt1,ebendpt2,brate1,brate2,brate3) %i

w=(1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1));
% Yhat=(rhat+delta)/(1-alpha);

earnings=z*w*l*epsilon_j;
u=phi*w*epsilon_j*l; % Unemployment benefits are at replacement rate phi;
unemploymentbenefit=(1-z)*u;

b_ebar=0; % b_ebar is defined in HV2000, pg 376
if ebar>0
    if ebar<ebendpt1
        b_ebar=b_ebar+brate1*ebar;
    else
        b_ebar=b_ebar+brate1*ebendpt1;
        if ebar<ebendpt2
            b_ebar=b_ebar+brate2*(ebar-ebendpt1);
        else
            b_ebar=b_ebar+brate2*(ebendpt2-ebendpt1)+brate3*(ebar-ebendpt2);
        end
    end
end
socialsecuritybenefit=b+b_ebar/((1+g)^(agej-Jr));

Tax_except_tau_c=tau_k*r*k+(tau_l+tau_s+tau_u)*w*l*epsilon_j*z; % IIJ working paper (1999) does not include that the labor income taxes are payed only when labor income is actually earned (their eqn 6), but I assume this is a typo

c=(1+r)*k+earnings-Taxes_except_tau_c+unemploymentbenefit+b*socialsecuritybenefit+Tr_beq-kprime; % Before consumption tax
c=c/(1+tau_c); % Consumption tax

F=-Inf;

if c>0
    F=((c^(varphi) *(1-l)^(1-varphi)).^(1-gamma))/(1-gamma);
end

% Borrowing constaint is imposed via grid on assets

% Endogenous retirement decision
if agej<Jr-1
    if bprime==1
        F=-Inf; % Cannot retire before reaching retirement age Jr (since this is next period you can choose retirement at Jr-1)
    end
else
    if b==1 && bprime~=1
        F=-Inf; % Once retired you cannot unretire
    end
end

% % Determine ebar, the average income 
% if b==1 % retired
%     if ebarprime==ebar
%         % ebar must remain constant during retirement
%     else 
%         F=-Inf;
%     end
% else % b=0 not yet retired
%     ebardist=abs(ebarprime-(ebar+min(earnings,emax)/agej));
%     % Note that earnings is zero if not working. IIJ2003 say they follow
%     % Huggett & Venture (1999) to define the evolution of ebar. But there
%     % is no unemployment in ebar, so unclear how exactly to treat this.
%     % I could follow the social security more closely (average of best X
%     % years) but this would require another variable that counts number of
%     % years worked.
%     if ebardist>ebarspacing/2 % If eprime is not the closed point on ebar_grid to what it should be
%         F=-Inf;
%     end
% end

end