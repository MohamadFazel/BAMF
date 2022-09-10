function produceIm(InIm,Name)
%produceIm() takes an iput image and convert it to a nice coloful image.
%
% INPUTS:
%   InIm:    Input image.
%   Name:    A string containing the name that the user want the image to
%            be saved with that.
%
% OUTPUTS:
%  
% REQUIRES:
%   MATLAB 2014a or higher.
%
% CITATION:
%   Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf,
%   Marjolein M.B.M. Meddens and Keith A. Lidek, "Bayesian Multiple Emitter
%   Fitting using Reversible Jump Markov Chain Monte Carlo".
%
% Created by:
%   Mohamadreza Fazel (Lidke Lab 2018)
%

     P = prctile(InIm(InIm>0),99.8);
     InIm(InIm>P) = P;
     InIm = InIm/P;
     Str = sprintf('%s.png',Name);
     imwrite(InIm*255,hot(256),Str);
     %The final image will be saved in your current directory.
end
 
 
 