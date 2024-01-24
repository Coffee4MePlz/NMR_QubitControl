function [summ,prodd,subb] = sumprod(x1,x2)
% usage:
% [summ, prodd] = sumprod(x1,x2)
% x1 = a complex number
% x2 = another complex number
% summ = sum of x1 and x2
% prodd = product pf x1 and x2
prodd = (x1)*(x2);
summ = x1 + x2;
subb = x1 - x2;