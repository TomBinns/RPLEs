function F = fScore(FP,TP,FN,beta)

nom = (1+beta^2)*TP;
den = (1+beta^2)*TP + beta^2*FN + FP;

F = nom./den;


end