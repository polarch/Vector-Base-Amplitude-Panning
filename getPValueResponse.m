function pValue = getPValueResponse(f, DTT)
%GETPVALUERESPONSE Summary of this function goes here
%   Detailed explanation goes here

a1 = 4.5e-4;
a2 = 8.5e-5;
p0 = 1.5-0.5*cos( 4.7*tanh(a1*f) ) .* max(0, 1-a2*f);
pValue = (p0-2) * sqrt(DTT) + 2;

end
