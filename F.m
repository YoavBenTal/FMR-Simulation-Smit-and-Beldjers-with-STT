function [F] = F(x, z)
%CALCULATEFANDG calculates F and G according to formula in notebook

%calculates the F and G functions defined in notebook #2
%08-06-09

F = x .* log(x.^2 + z.^2) + 2 .* z .* atan(x./z);


