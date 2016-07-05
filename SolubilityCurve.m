
% Demonstration of solubility estimation and comparison to experimental
% data

clc
%clear all
close all

DataExp = [0	35.65;
10	35.72;
20	35.89;
30	36.09;
40	36.37;
50	36.69;
60	37.04;
70	37.46;
80	37.93;
90	38.47;
100	38.99]; % experimental data constisting of temoerature (°C, first column)
            % and solubility (grams solute per kilogram of water)
            
TempExp = DataExp(:,1);  % temperature [°C]
SolubExp = DataExp(:,2); % solubility [g/kg]

plot(TempExp,SolubExp,'o','LineWidth',2)
xlabel('Solubility / [g solute / kg water]')
ylabel('Temperature / °C')

MW_NaCl = 58.44; %[g/mol] Molar mass sodium chloride
MW_H2O = 18.02;  %[g/mol] Molar mass water

AmountH2O = 1000/MW_H2O; % [mol] amount of water
AmountNaCl = SolubExp/MW_NaCl; % [mol] amount of sodium chloride
MolfracExp = AmountNaCl./(AmountNaCl+AmountH2O); %mole fraction of sodium chloride
TempExp_K = TempExp+273.15; %Experimental temperature in Kelvin
figure
plot(TempExp_K,MolfracExp,'LineWidth',2)
xlabel('Temperature / K')
ylabel('x_{NaCl}')

%% Plot of the ideal solubility as a function of temperature
Hfus_NaCl = 27950; %[J/mol] fusion enthalpy of sodium chloride
R = 8.3145;        %[J/mol/K] gas constant
Tm_NaCl = 801 + 273; %[K] Melting temperature of NaCl
x_ideal = exp(Hfus_NaCl./R.*(1./Tm_NaCl-1./TempExp_K));

hold on
plot(TempExp_K,x_ideal,'r','LineWidth',2)

%% Testing the function
val = NonidealSolubilityProblem( 0.7,Hfus_NaCl,Tm_NaCl,TempExp_K,MolfracExp);
Asol = fminsearch(@(A)NonidealSolubilityProblem(A,Hfus_NaCl,Tm_NaCl,TempExp_K,MolfracExp),0.9); %Finding A parameter for the Porter expression activity coefficient
x_sol = zeros(1,length(TempExp_K));
for k=1:numel(TempExp_K)
Funi = @(x) exp(Hfus_NaCl./R.*(1./Tm_NaCl-1./TempExp_K(k))-Asol*(1-x)^2)-x; % This expression is zero for the correct value of x.
x_sol(k) = fzero(Funi,0.02);
%x_sol2 = NewtonScalar(Funi,0.5)
%err = err+ (x_sol - xreal(k))^2;
end
plot(TempExp_K,x_sol)
%% Fitting the Porter parameter to the data
