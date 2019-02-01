function SINRCorrelatedRealMMSE = Func_CorrelatedSMMSE(InveseEstPhi, CorrelatedFading, DataPowerMatix, PilotPowerMatrix,nbrBSs,K,tau)
%This function computes SINR values of all users in case of standard (real)
%MMSE estimator
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


SINRCorrelatedRealMMSE = zeros(nbrBSs,K);% Space for SINRvalues;
Bmatrix = zeros(nbrBSs, nbrBSs, K);
Csmallmatrix =  zeros(nbrBSs, K, nbrBSs, K); %c_{jk}^it in the paper
Dmatrix = zeros(nbrBSs, K);

% Compute Bmatrix
for j = 1 : nbrBSs
    for i = 1: nbrBSs
        for k = 1:K
            Bmatrix(j,i,k) = sqrt(tau*PilotPowerMatrix(i,k)*PilotPowerMatrix(j,k))*abs(trace(squeeze(InveseEstPhi(j,k,:,:))*squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(CorrelatedFading(j,i,k,:,:))));
        end
    end
end

% Compute Cmatrix
for i = 1 : nbrBSs
    for t = 1 : K
        for j = 1 : nbrBSs
            for k = 1 : K
                Csmallmatrix(i,t,j,k) = PilotPowerMatrix(i,k)*abs(trace(squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(InveseEstPhi(j,k,:,:))*squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(CorrelatedFading(j,i,t,:,:))));
            end
        end
    end
end

% Compute Dmatrix
for j = 1 : nbrBSs
    for k = 1 : K
        Dmatrix(j,k) = PilotPowerMatrix(j,k)*abs(trace(squeeze(InveseEstPhi(j,k,:,:))*squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(CorrelatedFading(j,j,k,:,:))));
    end
end


for l = 1 : nbrBSs
    for k = 1:K
        Numerator = DataPowerMatix(l,k)*Bmatrix(l,l,k)^2;
        Denominator = Dmatrix(l,k);
        for i = [1:l-1,l+1:nbrBSs]
            Denominator = Denominator + DataPowerMatix(i,k)*Bmatrix(l,i,k)^2;
        end
        for i = 1: nbrBSs
            for  t = 1 :K
                Denominator  = Denominator + DataPowerMatix(i,t)*Csmallmatrix(i,t,l,k);
            end
        end
        SINRCorrelatedRealMMSE(l,k) = Numerator/Denominator;
    end
end
