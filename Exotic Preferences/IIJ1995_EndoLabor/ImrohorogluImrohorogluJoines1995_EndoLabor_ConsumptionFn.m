function c=ImrohorogluImrohorogluJoines1995_EndoLabor_ConsumptionFn(l,kprime,k,z,r,w,tau_u, tau_s,h,zeta,theta,epsilon_j,I_j,SSdivw, Tr_beq,workinglifeincome,g,agej,MedicalShock,LumpSum)

earnings=z*w*l*h*epsilon_j*I_j;
u=zeta*w*h; % Unemployment benefits are at replacement rate phi;
unemploymentbenefit=(l==0)*u*I_j;
socialsecuritybenefit=w*SSdivw*(1-I_j);
if g>0 % So not using baseline model
    socialsecuritybenefit=w*SSdivw*(1-I_j)*workinglifeincome*(1/((1+g)^agej));
    % Note: the actual social security benefit is constant, but the model
    % has been detrended by (1+g)^t and it is for these reason it here
    % appears to decrease with age.
end
% q=earnings+unemploymentbenefit+socialsecuritybenefit; % This is the notation used by IIJ1995 and is just here to ease comprehension. 

% IIJ1995 includes an extension to 'medical shocks' (see pg 110 of IIJ1995)
medicalexpense=0;
if MedicalShock>0
    % In retirement z now becomes the 'cost of illness' (when non-zero valued, the household is ill, when zero valued the household is healthy)
    % The cost of illness is measured as a fraction of employed wage (which is w*h)
    medicalexpense=(1-I_j)*z*w*h; % Note that this will be zero when working age (as I_j will be 1).
end

% Fixed cost of work
fc=0;
if l>0
    fc=theta;
end

c=(1+r)*k+(1-tau_u-tau_s)*earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq-medicalexpense-kprime+LumpSum-fc;
    
end