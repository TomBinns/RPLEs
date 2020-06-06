function noise = genNoise(alpha,base)
% genNoise - Converts Gaussian pseudo-random numbers into noise with a 
%     specific spectral density as determined by the alpha value of the 
%     noise. Adapted from Schurger (2018)'s code
% 
% Argument(s):
%  alpha -    desired alpha value of the PSD of the noise
%  base -     colums vector of Gaussian pseudo-random numbers representing 
%             a whole experimental-condition's worth of trials 
%             (single-channel data)
%             
% Returns:
%  noise -    column vector of noise with a specific spectral density
%  
% Author(s): Aaron Schurger; Thomas Binns, 2020
% 
% Reference(s): Schurger A (2018). Specific relationship between the shape 
%     of the readiness potential, subjective decision time, and waiting 
%     time predicted by an accumulator model with temporally autocorrelated 
%     input noise. Eneuro 5.


noise = nan(size(base));

if alpha > 0
    for aa = 1:size(base,2)
        noise(:,aa) = changeAlpha(base(:,aa),alpha);
    end
else
    noise = base;
end

noise = double(noise);

end