function MCMC=threshRJ(ClustInfo,SZ,Zoom)
Ind = ClustInfo.X_SE < 0.2 & ClustInfo.X_SE > 0.005 & ClustInfo.Y_SE < 0.2 ...
      & ClustInfo.Y_SE > 0.005;
      SMD.X = Zoom*ClustInfo.X(Ind);
      SMD.Y = Zoom*ClustInfo.Y(Ind);
      SMD.I = ClustInfo.I(Ind);
      SMD.X_SE = Zoom*ClustInfo.X_SE(Ind);
      SMD.Y_SE = Zoom*ClustInfo.Y_SE(Ind);
      MCMC = zeros(Zoom*SZ,'single');
      [Xg,Yg]=meshgrid((0.5:SZ*Zoom),(0.5:SZ*Zoom));
      for nn = 1:length(SMD.X)
          MCMC = MCMC+SMD.I(nn)*normpdf(Xg,SMD.X(nn),SMD.X_SE(nn)).*normpdf(Yg,SMD.Y(nn),SMD.Y_SE(nn));
      end
end