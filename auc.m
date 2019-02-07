% vectors of 1 with the same length as x and y, respectively
%AUC function (using Simpson’s method);
function [auc] =auc(roctrue,grid) %input ROC curve vector as roctrue.
auc=1/3*grid*(roctrue(1)+roctrue(length(roctrue))+2*sum(roctrue(2:(length(roctrue)-1)))+2*sum(roctrue(2:2:(length(roctrue)-1))));
end