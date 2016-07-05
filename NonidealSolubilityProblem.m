function [ err ] = NonidealSolubilityProblem( A,Hfus,Tm,T,xreal )
% In: Guess for the Porter parameter A
% Enthalpy of fusion Hfus [J/mol]
% Melting temperature Tm [K]
% Temperature T [K]
% xreal 

% Out: Error/deviation of function from data

err = 0; %initializing the total error for the chosen A
R = 8.3145; %[J/mol/K] Gas constant
for k=1:numel(T)
Funi = @(x) exp(Hfus./R.*(1./Tm-1./T(k))-A*(1-x)^2)-x; % This expression is zero for the correct value of x.
x_sol = fzero(Funi,0.02);
%x_sol2 = NewtonScalar(Funi,0.5)
err = err+ abs(x_sol - xreal(k))*1000;
end
end

