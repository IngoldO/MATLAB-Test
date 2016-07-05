function [ sol ] = NewtonScalar( funName, x3 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x2 = x3;
x1 = x2;


Tol = 0.001;

while abs(funName(x3))>Tol
x3 = x2 - funName(x2)*(x2-x1)/(funName(x2)-funName(x1));
x1 = x2;
x2 = x3;
end
 sol = x3;   
end

