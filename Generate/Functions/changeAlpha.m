function noise = changeAlpha(data,alpha)
% changeAlpha - Converts Gaussian pseudo-random numbers into noise with a 
%     specific spectral density as determined by the alpha value of the 
%     noise. Adapted from Schurger (2018)'s code
% 
% Argument(s):
%  data -     column vector of Gaussian pseudo-random numbers
%  alpha -    desired alpha value of the PSD of the noise
%             
% Returns:
%  noise -    row vector of noise with a specific spectral density
%  
% Author(s): Aaron Schurger; Thomas Binns, 2020
% 
% Reference(s): Schurger A (2018). Specific relationship between the shape 
%     of the readiness potential, subjective decision time, and waiting 
%     time predicted by an accumulator model with temporally autocorrelated 
%     input noise. Eneuro 5.


N2 = floor(size(data,1)/2)-1;
f = (2:(N2+1))';
A2 = 1./(f.^(alpha/2));
d = [1;A2;1/((N2+2)^alpha);flipud(A2)];
z = fft(data);
z = z.*sqrt(d.^2);
noise = ifft(z);

end