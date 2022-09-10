function SMR=threhsSMA(SMD,SMAStruct,BoxSize)
if ~isfield(SMAStruct,'MaxXY_SE')
    SMAStruct.MaxXY_SE = 0.2;
end
Pvalue = pValue('Basic',BoxSize,SMD.LogLikelihood);
MinMax.XY_SE=[0 SMAStruct.MaxXY_SE];
MinMax.Pvalue=[SMAStruct.MinPValue 1000000000];
MinMax.Photons = [SMAStruct.MinPhotons 100000000];
Ind = SMD.X_SE < MinMax.XY_SE(2) & isreal(SMD.X_SE) & SMD.Y_SE < MinMax.XY_SE(2) & ...
    isreal(SMD.Y_SE) & Pvalue > MinMax.Pvalue(1) & isreal(Pvalue) & SMD.Photons > MinMax.Photons(1) & ...
    isreal(SMD.Photons) & isreal(SMD.X) & isreal(SMD.Y) & SMD.X > 0 & SMD.Y > 0 & ...
    SMD.X < BoxSize & SMD.Y < BoxSize & SMD.BG > 0;
SMR.Photons = SMD.Photons(Ind);
SMR.Bg = SMD.BG(Ind); 
end

function [Pvalue]=pValue(FitType,FitBoxSize,LL)
X2_CDF=inline('gammainc(x/2,k/2)','k','x');
  switch FitType
      case 'Basic'
           k=FitBoxSize^2-4;
      case {'Sigma','Astig','SampledZ'}
           k=FitBoxSize^2-5;
      otherwise
           error('SR_demo:Threshold','Unknown FitType');
  end
X2=-2*LL;
Pvalue=1-X2_CDF(k,X2);
end