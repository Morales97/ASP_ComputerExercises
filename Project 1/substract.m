function [Num,Den] = substract(Num1,Den1,Num2,Den2)
%
% SUBSTRACT 	Substract Z-transforms (H1(z) - H2(z))
%
% [Num,Den]=substract(Num1,Den1,Num2,Den2)
%

% Make numerators and denominator transforms have equal length

[Num1,Den1]=eqsize(Num1,Den1);
[Num2,Den2]=eqsize(Num2,Den2);

% Substact

Num=conv(Num1,Den2) - conv(Den1,Num2);
Den=conv(Den1,Den2);

% Remove zeros

[Num,Den]=rmczeros(Num,Den);
