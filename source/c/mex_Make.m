%% This acts as a make file for mex functions. 

% This should be run to generate mex files for each architecture

%% cRJMCMC

%using visual studio to compile the mexFunction. You need to have visual
%studio installed
%mex(fullfile(sourcePath,'cRJMCMC.cpp'),fullfile(sourcePath,'LibRJMCMC.cpp'));  
mex('cRJMCMC.cpp','LibRJMCMC.cpp');  
%using g++ to compile mexFunction. You need to have MinGW installed
%cd(fullfile(sourcePath,'g++'))
%mex -v GPP='C:/MinGW/bin/g++' cRJMCMC.cpp LibRJMCMC.cpp;  

%movefile(fullfile(basePath,['cRJMCMC.' mexext]), mexFilePath);
movefile(['cRJMCMC.' mexext], fullfile('..','..','mex'));