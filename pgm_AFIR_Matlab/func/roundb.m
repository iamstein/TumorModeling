function[output] = roundb(input,n)
%ROUNDB - rounds input to the nearest nth digit
%if n=2, round(123) = 100
%if n=-2, round(123.456) = 123.45

output = round(input*10.^-n)*10.^n;