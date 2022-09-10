addpath('.')
if ispc
   % Adding system path for nvcc to compile with nvcc
   setenv('PATH', [getenv('PATH') ';C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\bin']);
   % Adding system path for VS2013 to compile with cl
   setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin']);
else % Linux/MacOS
   % Adding system path for nvcc to compile with nvcc
   setenv('PATH', [getenv('PATH') ':/usr/local/cuda-8.0/bin']);
end

cd source/c
fprintf('Running mex_Make\n');
mex_Make
cd ../..
cd source/cuda
fprintf('Running cuda_Make\n');
cuda_Make
cd ../..
