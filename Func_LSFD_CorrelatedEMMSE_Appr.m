function SINRLSFDCorrelatedEMMSE_Appr = Func_LSFD_CorrelatedEMMSE_Appr(lossovernoise, EstError,EstPhi, CorrelatedFading,DataPowerMatix, PilotPowerMatrix, nbrBSs,K,tau, NBScases)
%This function computes SINR values of all users in case of element-wise
%MMSE estimator and approximated LSFD vectors
%
%This Matlab function was developed to generate simulation results in
%
%Trinh Van Chien, Christopher Mollen and Emil Bjornson,
%"Large-Scale-Fading Decoding in Cellular Massive MIMO Systems with
%Spatially Correlated Channels", IEEE Transactions on Communications,
%Accepted for publication.
%
%This is version 1.0 (Last edited: 2018-12-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


SINRLSFDCorrelatedEMMSE_Appr = zeros(nbrBSs,K);% Space for SINRvalues;
Bmatrix = zeros(nbrBSs, nbrBSs, K);
Csmallmatrix =  zeros(nbrBSs, K, nbrBSs, K); %c_{jk}^it in the paper
Dmatrix = zeros(nbrBSs, K);

% Compute Bmatrix
for j = 1 : nbrBSs
    for i = 1: nbrBSs
        for k = 1:K
            Bmatrix(j,i,k) = sqrt(tau)*EstError(j,j,k)*EstError(j,i,k)*abs(trace(squeeze(EstPhi(j,k,:,:))));
        end
    end
end

% Compute Cmatrix
for i = 1 : nbrBSs
    for t = 1 : K
        for j = 1 : nbrBSs
            for k = 1 : K
                Csmallmatrix(i,t,j,k) = ((EstError(j,j,k))^2)*abs(trace(squeeze(CorrelatedFading(j,i,t,:,:))*squeeze(EstPhi(j,k,:,:))));
            end
        end
    end
end

% Compute Dmatrix
for j = 1 : nbrBSs
    for k = 1 : K
        Dmatrix(j,k) = ((EstError(j,j,k))^2)*abs(trace(squeeze(EstPhi(j,k,:,:))));
    end
end

% Compute Approximation
Bmatrix_Appr = zeros(nbrBSs, nbrBSs, K);
Csmallmatrix_Appr =  zeros(nbrBSs, K, nbrBSs, K); %c_{jk}^it in the paper
Dmatrix_Appr = zeros(nbrBSs, K);

% Compute Approximatation of Bmatrix
for j = 1 : nbrBSs
    for i = 1: nbrBSs
        for k = 1:K
            Bmatrix_Appr(j,i,k) = sqrt(NBScases*tau*PilotPowerMatrix(i,k))*EstError(j,j,k)*lossovernoise(j,i,k);
        end
    end
end

% Compute Approximation of Cmatrix
for i = 1 : nbrBSs
    for t = 1 : K
        for j = 1 : nbrBSs
            for k = 1 : K
                Csmallmatrix_Appr(i,t,j,k) = (EstError(j,j,k))*sqrt(PilotPowerMatrix(j,k))*lossovernoise(j,j,k)*lossovernoise(j,i,t);
            end
        end
    end
end

% Compute Approximation of Dmatrix
for j = 1 : nbrBSs
    for k = 1 : K
        Dmatrix_Appr(j,k) = (EstError(j,j,k))*sqrt(PilotPowerMatrix(j,k))*lossovernoise(j,j,k);
    end
end

% Compute SINR values
for l = 1 : nbrBSs
    for k = 1:K
        % compute the first part of Clk
        Clk = zeros(nbrBSs,nbrBSs);
        for i = [1:l-1,l+1:nbrBSs]
            Clk = Clk + DataPowerMatix(i,k)*squeeze(Bmatrix(:,i,k))*squeeze(Bmatrix(:,i,k))';
        end
        
        % compute the second part of Clk
        ctemplk = zeros(nbrBSs,1);
        for s = 1: nbrBSs
            ctemplk(s) = Dmatrix(s,k);
            for i =1 : nbrBSs
                for t = 1: K
                    ctemplk(s) = ctemplk(s) + DataPowerMatix(i,k)*Csmallmatrix(i,t,s,k);
                end
            end
        end
        Clk = Clk + diag(ctemplk);
        % ======Compute the suboptimal a_lk ==========================
        % compute the first part of Approximated Clk
        Clk_Appr = zeros(nbrBSs,nbrBSs);
        for i = [1:l-1,l+1:nbrBSs]
            Clk_Appr = Clk_Appr + DataPowerMatix(i,k)*squeeze(Bmatrix_Appr(:,i,k))*squeeze(Bmatrix_Appr(:,i,k))';
        end
        
        % compute the second part of Approximated Clk
        ctemplk_Appr = zeros(nbrBSs,1);
        for s = 1: nbrBSs
            ctemplk_Appr(s) = Dmatrix_Appr(s,k);
            for i =1 : nbrBSs
                for t = 1: K
                    ctemplk_Appr(s) = ctemplk_Appr(s) + DataPowerMatix(i,k)*Csmallmatrix_Appr(i,t,s,k);
                end
            end
        end
        Clk_Appr = Clk_Appr + diag(ctemplk_Appr);
        alk_Appr = Clk_Appr\squeeze((Bmatrix_Appr(:,l,k)));
        Numerator_Appr = DataPowerMatix(l,k)*abs(alk_Appr'*squeeze((Bmatrix(:,l,k))))^2;
        Denumerator_Appr = squeeze(alk_Appr'*Clk*alk_Appr);
        SINRLSFDCorrelatedEMMSE_Appr(l,k) = Numerator_Appr/Denumerator_Appr;
        
    end
end
