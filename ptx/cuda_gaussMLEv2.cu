/*!
 * \file GPUgaussMLEv2.cu
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief This file contains all of the Cuda kernels.  The helper functions
 * are defined in GPUgaussLib.cuh
 */

#include "definitions.h"
#include "MatInvLib.h"
#include "GPUgaussLib.cuh"
#include "GPUgaussMLEv2.h"

//*******************************************************************************************
//theta is: {x,y,N,bg}
__global__ void kernel_MLEFit_XYNB_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[NV_P*NV_P], Diag[NV_P], Minv[NV_P*NV_P];
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    const int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_P;
    float dudt[NV_P];
    float d2udt2[NV_P];
    float NR_Numerator[NV_P], NR_Denominator[NV_P];
    float theta[NV_P];
    float maxjump[NV_P]={1e0f, 1e0f, 1e2f, 2e0f};
    float Nmax;

    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_P*NV_P*sizeof(float));
	memset(Minv,0,NV_P*NV_P*sizeof(float));
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma); //Added 2* on 8.9.16 to account for smoothing filter.
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
		memset(NR_Numerator,0,NV_P*sizeof(float));
		memset(NR_Denominator,0,NV_P*sizeof(float));

        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
            PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], &d2udt2[1]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
         
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
        PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBS_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    
    //__shared__ float s_data[MEM];
    float M[NV_PS*NV_PS], Diag[NV_PS], Minv[NV_PS*NV_PS];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_PS;
    float dudt[NV_PS];
    float d2udt2[NV_PS];
    float NR_Numerator[NV_PS], NR_Denominator[NV_PS];
    float theta[NV_PS];
    float maxjump[NV_PS]={1e0f, 1e0f, 1e2f, 2e0f, 5e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_PS*NV_PS*sizeof(float));
	memset(Minv,0,NV_PS*NV_PS*sizeof(float));      
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV_PS*sizeof(float));
		memset(NR_Denominator,0,NV_PS*sizeof(float));
      
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], &d2udt2[4]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain Sigma
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]); //bug fix 8.9.16 
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]); //bug fix 8.9.16 
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBZ_(const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[5*5], Diag[5], Minv[5*5];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=5;
    float dudt[5];
    float d2udt2[5];
    float NR_Numerator[5], NR_Denominator[5];
    float theta[5];
    float maxjump[5]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
    theta[4]=0;
   
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, dudt, d2udt2);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating remaining derivatives
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -maxjump[4]), maxjump[4]);
        
        // Any other constraints
        theta[2]=max(theta[2], 1.0f);
        theta[3]=max(theta[3], 0.01f);
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay, Bx,By, gamma, d, &PSFx, &PSFy, dudt, NULL);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating remaining derivatives
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
       
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) 
    kernel_MatInvN(M, Minv, Diag, NV);
  
   //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBSXSY_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
 
    //__shared__ float s_data[MEM];
    float M[6*6], Diag[6], Minv[6*6];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=6;
    float dudt[6];
    float d2udt2[6];
    float NR_Numerator[6], NR_Denominator[6];
    float theta[6];
    float maxjump[6]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f,1e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    
	//initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    theta[5]=PSFSigma;
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
   
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], &d2udt2[4]);
            kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], &d2udt2[5]);
            
            
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);
        theta[5]-=min(max(NR_Numerator[5]/NR_Denominator[5], -theta[5]), theta[5]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain SigmaX
        theta[5]=max(theta[5], 0.5f); //Constrain Sigma
        theta[5]=min(theta[5], sz/2.0f); //Constrain SigmaX
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], NULL);
        kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
   
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

// SCMOS Versions---------------------------------------

__global__ void kernel_MLEFit_SCMOSXYNB_(const float *d_data, const float *d_Coords, const float *d_GainRatio, 
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[NV_P*NV_P], Diag[NV_P], Minv[NV_P*NV_P];
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    const int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_P;
    float dudt[NV_P];
    float d2udt2[NV_P];
    float NR_Numerator[NV_P], NR_Denominator[NV_P];
    float theta[NV_P];
    float maxjump[NV_P]={1e0f, 1e0f, 1e2f, 2e0f};
    float Nmax;
	float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_P*NV_P*sizeof(float));
	memset(Minv,0,NV_P*NV_P*sizeof(float));
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
	const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
		memset(NR_Numerator,0,NV_P*sizeof(float));
		memset(NR_Denominator,0,NV_P*sizeof(float));

        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
            PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			//GRind=(int)s_Coords[0];
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], &d2udt2[0]);//x
            kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], &d2udt2[1]);//y
            dudt[2] = PSFx*PSFy;// I
            d2udt2[2] = 0.0f;// I
            dudt[3] = 1.0f;// bg
            d2udt2[3] = 0.0f;// bg
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR); // add variance-gain ratio: v/g^2
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2); // add variance-gain ratio: v/g^2
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
        PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/(model+gainR);// add gain ratio
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;// add gain ratio
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}


//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBS_(const float *d_data, const float *d_Coords, const float *d_GainRatio,
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view.
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    
    //__shared__ float s_data[MEM];
    float M[NV_PS*NV_PS], Diag[NV_PS], Minv[NV_PS*NV_PS];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_PS;
    float dudt[NV_PS];
    float d2udt2[NV_PS];
    float NR_Numerator[NV_PS], NR_Denominator[NV_PS];
    float theta[NV_PS];
    float maxjump[NV_PS]={1e0f, 1e0f, 1e2f, 2e0f, 5e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_PS*NV_PS*sizeof(float));
	memset(Minv,0,NV_PS*NV_PS*sizeof(float));      
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV_PS*sizeof(float));
		memset(NR_Denominator,0,NV_PS*sizeof(float));
      
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], &d2udt2[4]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain Sigma
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/(model+gainR);
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBZ_(const float *d_data, const float *d_Coords, const float *d_GainRatio, const float *d_x0,
		const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[5*5], Diag[5], Minv[5*5];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=5;
    float dudt[5];
    float d2udt2[5];
    float NR_Numerator[5], NR_Denominator[5];
    float theta[5];
    float maxjump[5]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
	const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	const float *z_initial = d_x0+(bx*BlockSize+tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
    theta[4]=z_initial[0];
   
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, dudt, d2udt2);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating remaining derivatives
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -maxjump[4]), maxjump[4]);
        
        // Any other constraints
        theta[2]=max(theta[2], 1.0f);
        theta[3]=max(theta[3], 0.01f);
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay, Bx,By, gamma, d, &PSFx, &PSFy, dudt, NULL);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating remaining derivatives
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
       
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) 
    kernel_MatInvN(M, Minv, Diag, NV);
  
   //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBSXSY_(const float *d_data, const float *d_Coords, const float *d_GainRatio, 
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
 
    //__shared__ float s_data[MEM];
    float M[6*6], Diag[6], Minv[6*6];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=6;
    float dudt[6];
    float d2udt2[6];
    float NR_Numerator[6], NR_Denominator[6];
    float theta[6];
    float maxjump[6]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f,1e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	//initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    theta[5]=PSFSigma;
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], &d2udt2[4]);
            kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], &d2udt2[5]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);
        theta[5]-=min(max(NR_Numerator[5]/NR_Denominator[5], -theta[5]), theta[5]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain SigmaX
        theta[5]=max(theta[5], 0.5f); //Constrain Sigma
        theta[5]=min(theta[5], sz/2.0f); //Constrain SigmaX
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], NULL);
        kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
   
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

