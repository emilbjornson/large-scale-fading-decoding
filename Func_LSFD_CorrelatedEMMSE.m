function SINRLSFDCorrelatedElementMMSE = Func_LSFD_CorrelatedEMMSE(EstError,EstPhi, CorrelatedFading,DataPowerMatix, nbrBSs,K,tau)
%This function computes SINR values of all users in case of element-wise
%MMSE estimator and optimal LFSD vectors
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


SINRLSFDCorrelatedElementMMSE = zeros(nbrBSs,K);% Space for SINRvalues;
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
        SINRLSFDCorrelatedElementMMSE(l,k) = DataPowerMatix(l,k)*squeeze(Bmatrix(:,l,k))'*((Clk)\squeeze(Bmatrix(:,l,k)));
        
    end
end
