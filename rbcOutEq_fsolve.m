function [r,J] = rbcOutEq_fsolve(p,y)
r = rbcOutEq_f(p, y);
J = rbcOutEq_ff(p, y);
