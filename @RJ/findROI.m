      function [SMD, ROIStack, Res]=findROI(Image,Sigma1,Sigma2,Boxsize,Minval)
            %   findROI() gets a sequence of super-resolution data as 
            %   input, which is processed in data2Photons() by taking into
            %   account the gain and offset. It finds particles and then
            %   crops and return them. The cropped region is a box with the
            %   size of boxsize. The code uses gpu function to speed up 
            %   calculations.
            %   Three nested functions are included:
            %   gauss_inplace() implements the gaussian filters.
            %   subtract_inplace() subtracts the outcomes of the
            %   gaussian_inplace functions with two different Gaussian
            %   filters.
            %   local_max() finds the maximums of the outcome of the
            %   previous nested function (gaussian_inplace), which are a
            %   rough estimation of the particles' locations.
            %   
            % INPUTS:
            %   SMD:        A structure that some fileds will be added to 
            %               it. (No default)
            %   Image:      The sequence of super-resolution data.
            %               (No default).
            %   Sigma1:     The sigma of the first Gaussian filter.(Pixels)
            %               (Default = 1.3).
            %   Sigma2:     The sigma of the second Gaussina filter. It is 
            %               usually two times of sigma1. (Pixels)
            %               (Default = 2*1.3)
            %   Boxsize:    The size of the returned boxes. (Pixels) 
            %               (Default = 7)
            %   Minval:     A limit related to the intesity.
            %               (Default = (1/4)*2000*(1/(2*pi*2*1.3^2))).   
            %
            % OUTPUTS:
            %   SMD:        The input structure, which the following fields
            %               are added to it:
            %   XBoxCorner: The X-coordinate of the upper left corner of 
            %               the box.(Pixels),(Number of found particles x1)
            %   YBoxCorner: The Y-ccordinate of the upper left corner of 
            %               the box.(Pixels),(Number of found particles x1)
            %   FrameNum:   The corresponding frame numbers where the 
            %               particles were found.
            %               (Number of found particles x1)
            %
            %   ROIStack:   The found boxes. 
            %               (Boxsize x Boxsize x Number of found particles)
            %   Res:        A binary image the same size as the input image
            %               where the bright pixels represent the positions 
            %               of the found particles.
            %               (The same size of the input sequence)
            %
            % REQUIRES:
            %   MATLAB 2014a or later versions
            %   Parallel Procesing Toolbox
            %   NVidia GPU
            %   cuda_FindROI.ptx
            %   cuda_FindROI.cu
            %
            % CITATION:
            %   Mohamadreza Fazel (Lidke Lab, 2017)
            
            %Default values
            if nargin < 1
               error('The two first inputs are required.'); 
            end
            if nargin < 2
                Sigma1 = 1.3;
            end
            if nargin < 3
                Sigma2 = 2*1.3;
            end
            if nargin < 4
                Boxsize = 7;
            end
            if nargin < 5
               Minval = (1/4)*2000*(1/(2*pi*2*1.3^2));
            end
           
            %Creating GPU CUDA kernel objects from PTX and CU code
            k1 = parallel.gpu.CUDAKernel('cuda_FindROI.ptx','cuda_FindROI.cu','gaussX');
            k2 = parallel.gpu.CUDAKernel('cuda_FindROI.ptx','cuda_FindROI.cu','gaussY');
            k3 = parallel.gpu.CUDAKernel('cuda_FindROI.ptx','cuda_FindROI.cu','subtract');
            k4 = parallel.gpu.CUDAKernel('cuda_FindROI.ptx','cuda_FindROI.cu','maxX');
            k5 = parallel.gpu.CUDAKernel('cuda_FindROI.ptx','cuda_FindROI.cu','maxY');

            %getting size of the input image
            Im_size = size(Image);
            Ysize = size(Image,1);
            Xsize = size(Image,2);
            Zsize = size(Image,3);
            %number of pixels in the input sequence
            Nelem = Xsize*Ysize*Zsize;

            %gpuDevice() gives the info of your hardware (graphic card)
            g = gpuDevice;
            %The whole set of data cannot be sent into the gpu at the same
            %time because of the limitation of the gpu memory. Here, we are
            %finding the number of chuncks of data based on the available
            %memory (g.Totalmemory) to send them to the gpu.
            Nloops = ceil(4*4*Nelem/(g.TotalMemory));
            %number of frames in each chunk.
            Sz2_inst = floor(Zsize/Nloops);
            Res = zeros(Im_size);

            %The loop sends one chunk of data at each iteration to the gpu.
            for ii = 1:Nloops

                Aa = ii*Sz2_inst;
                if ii == Nloops
                    Aa = Zsize;
                end
                %The chunk of data that are sent to the gpu
                ImageChunk = Image(:,:,(ii - 1)*Sz2_inst + 1:Aa);
                %Here are some nested functions, where the number of
                %threads are set up and the gpu codes are called.
                [D_A] = gauss_inplace(ImageChunk,Sigma1,Xsize,Ysize,Sz2_inst,k1,k2);
                [D_B] = gauss_inplace(ImageChunk,Sigma2,Xsize,Ysize,Sz2_inst,k1,k2);
                [SubIm] = subtract_inplace(D_A,D_B,Xsize,Ysize,Sz2_inst,k3);
                [ResChunk] = local_max(SubIm,D_A,Xsize,Ysize,Sz2_inst,Boxsize,Minval,k4,k5);
                %Using the gather command to retrieve output from gpu memory.
                %The output is binary image where bright pixels show the
                %particle positions.
                Res(:,:,(ii-1)*Sz2_inst+1:Aa) = gather(ResChunk);

            end
            %Finding the bright pixels.
            [Ab] = find(Res);
            [Ii,Jj,Kk] = ind2sub(size(Res),Ab);
            %centers = [Ii,Jj,Kk];
            
            %finding the boxes
            HalfBox = floor(Boxsize/2+0.5);
            StartRow = Ii-HalfBox+1;
            minusRow = find(StartRow<1);
            StartRow(minusRow) = 1;
            StartCol = Jj-HalfBox+1;
            MinusCol = find(StartCol<1);
            StartCol(MinusCol) = 1;

            EndRow = StartRow + Boxsize-1;
            EndRow(minusRow) = Boxsize;
            LargeRow = find(EndRow > Ysize);
            EndRow(LargeRow) = Ysize;
            StartRow(LargeRow) = Ysize - Boxsize +1;

            EndCol = StartCol + Boxsize-1;
            EndCol(MinusCol) = Boxsize;
            LargeCol = find(EndCol > Xsize);
            EndCol(LargeCol) = Xsize;
            StartCol(LargeCol) = Xsize - Boxsize + 1;
            
            %Outputs:
            %the coordinates of the left top corner of the boxes with the
            %number of frames.
            startcoords = [StartRow, StartCol, Kk];
            %number found boxes
            N = size(startcoords,1);
            %found boxes
            ROIStack = zeros(Boxsize,Boxsize,N);
            for bb = 1:N
                ROIStack(:,:,bb) = Image(StartRow(bb):EndRow(bb), StartCol(bb):EndCol(bb), Kk(bb));
            end
            SMD.XBoxCorner = StartCol;
            SMD.YBoxCorner = StartRow;
            SMD.FrameNum = Kk;
            
            %Nested functions:
            function [Im] = gauss_inplace(Image, Sigma, Xsize, Ysize, Zsize_sub, k1, k2)
                %gauss_inplace() apply a Gaussian to the image (along X and Y), where the
                %sigma of the Gaussina is an input.
                if Sigma > 2.5
                   Q = 0.98711*Sigma-0.96330; 
                else
                   Q = 3.97156-4.14554*sqrt(1-0.26891*Sigma);
                end

                Qq = Q*Q;
                Qqq = Qq*Q;

                B0 = 1.57825+2.44413*Q+1.4281*Qq+0.422205*Qqq;
                B1 = 2.44413*Q+2.85619*Qq+1.26661*Qqq;
                B2 = -1.4281*Qq-1.26661*Qqq;
                B3 = 0.422205*Qqq;
                B = 1 - (B1+B2+B3)/B0;
                
                %Setting up number of threads and blocks
                k1.GridSize(1) = Zsize_sub;
                k1.ThreadBlockSize(1) = Ysize;
                %Calling the gpu code to apply Gaussian filter along the
                %X-axis.
                Im0 = feval(k1,Image,Xsize,B0,B1,B2,B3,B);
                Im0 = gather(Im0);
                %Setting up the number of threads and blocks for the
                %gpu-code that applies the Gaussian function along the
                %Y-axis.
                k2.GridSize = Zsize_sub;
                k2.ThreadBlockSize = Xsize;
                %Calling the gpu-code that applies the Gaussian function
                %along the Y-axis.
                Im = feval(k2,Im0,Ysize,B0,B1,B2,B3,B);
                Im = gather(Im);
           end
           function [DIm]=subtract_inplace(Image1,Image2,Xsize,Ysize,Zsize_sub,k3)
                %subtract_inplace() gets the two filtered images with two
                %different Gaussians and subtract them to find the edges.
                k3.GridSize = [Ysize,Zsize_sub,1];
                k3.ThreadBlockSize(1) = Xsize;
                DIm = feval(k3,Image1,Image2);
                DIm = gather(DIm);

           end
           function [Res]=local_max(SubIm,Image1,Xsize,Ysize,Zsize_sub,Boxsize,Minval,k4,k5)
                %local_max() takes the output of the subtract-inplace
                %function and finds the particle positions using some
                %citeria like the intensity of the pixels.
                k4.GridSize = [Ysize,Zsize_sub,1];
                k4.ThreadBlockSize(1) = Xsize;
                [~,Im1] = feval(k4,SubIm,Image1,Boxsize,Minval);
                Im1 = gather(Im1);

                k5.GridSize = [Xsize,Zsize_sub,1];
                k5.ThreadBlockSize(1) = Ysize;
                Res = feval(k5,SubIm,Im1,Boxsize,Minval);
                Res = gather(Res);

           end
        end