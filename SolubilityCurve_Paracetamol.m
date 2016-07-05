
% Demonstration of solubility estimation and comparison to experimental
% data

clc
%clear all
close all


DataExp = [-5	174.48;
0	191.48;
5	215.09;
10	239.6;
15	265.43;
20	297.81;
25	332.11;
30	371.61]; % experimental data consisting of temperature (°C, first column)
             % and solubility (grams solute per kilogram of water)

            
            
TempExp = DataExp(:,1);  % temperature [°C]
SolubExp = DataExp(:,2); % solubility [g/kg]

plot(TempExp,SolubExp,'o','LineWidth',2)
ylabel('Solubility / [g Paracetamol / kg Methanol]')
xlabel('Temperature / °C')
axis([-7.5 32.5 150 400])

MW_PM = 151.163; %[g/mol] Molar mass paracetamol
MW_MeOH = 32.04;  %[g/mol] Molar mass methanol

AmountMeOH = 1000/MW_MeOH; % [mol] amount of water
AmountPM = SolubExp/MW_PM; % [mol] amount of sodium chloride
MolfracExp = AmountPM./(AmountMeOH+AmountPM); %mole fraction of sodium chloride
TempExp_K = TempExp+273.15; %Experimental temperature in Kelvin
figure
plot(TempExp_K,MolfracExp,'o','LineWidth',2)
xlabel('Temperature / K')
ylabel('x_{Paracetamol}')

%% Plot of the ideal solubility as a function of temperature
Hfus_PM = 27100; %[J/mol] fusion enthalpy of paracetamol
R = 8.3145;        %[J/mol/K] gas constant
Tm_PM = 443.6; %[K] Melting temperature of NaCl
x_ideal = exp(Hfus_PM./R.*(1./Tm_PM-1./TempExp_K));

hold on
plot(TempExp_K,x_ideal,'r','LineWidth',2)

%% Testing the function
val = NonidealSolubilityProblem( 0.7,Hfus_PM,Tm_PM,TempExp_K,MolfracExp);
Asol = fminsearch(@(A)NonidealSolubilityProblem(A,Hfus_PM,Tm_PM,TempExp_K,MolfracExp),0); %Finding A parameter for the Porter expression activity coefficient
x_sol = zeros(1,length(TempExp_K));
deviation = x_sol;
for k=1:numel(TempExp_K)
Funi = @(x) exp(Hfus_PM./R.*(1./Tm_PM-1./TempExp_K(k))-Asol*(1-x)^2)-x; % This expression is zero for the correct value of x.
x_sol(k) = fzero(Funi,0.02);
exp(Asol*(1-x_sol(k))^2)
deviation(k) = x_sol(k) - MolfracExp(k);
end
plot(TempExp_K,x_sol,'LineWidth',2)
legend('Experimental','Ideal','Nonideal')

figure
bar(TempExp_K,deviation)
xlabel('Temperature / K')
ylabel('Error')

