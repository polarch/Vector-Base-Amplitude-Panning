function pValue = getPValueResponse(f, DTT)
%GETPVALUERESPONSE Returns the normalization power (p-value) of the VBAP
% gains for a certain frequency, based on the experimental formula of
%
%   Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014). 
%   Gain normalization in amplitude panning as a function of frequency and room reverberance. 
%   55th International Conference of the AES. Helsinki, Finland.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 1/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = 4.5e-4;
a2 = 8.5e-5;
p0 = 1.5-0.5*cos( 4.7*tanh(a1*f) ) .* max(0, 1-a2*f);
pValue = (p0-2) * sqrt(DTT) + 2;

end
