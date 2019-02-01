function [SINROptCorrelatedEMMSE, DataPowermatrixOptMMSE, Asave] =  Func_Opt_CorrelatedEMMSE(IntDataPowerMatrix,EstError,EstPhi, CorrelatedFading, DataPowerMax,nbrBSs,K,tau,NumIter)
%This function optimizes SINR values of all users in case of element-wise
%MMSE estimation
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


SINROptCorrelatedEMMSE = zeros(nbrBSs,K);% Space for SINRvalues;
Bmatrix = Compute_Bmatrix(tau,EstError,EstPhi,nbrBSs,K);
Csmallmatrix = Compute_Cmatrix(EstError,EstPhi,CorrelatedFading,nbrBSs,K);
Dmatrix = Compute_Dmatrix(EstError,EstPhi,nbrBSs,K);
SqrtDataPowermatrix = IntDataPowerMatrix;
Asave = zeros(NumIter,1);
DataPowermatrixSave = zeros(nbrBSs,K,NumIter);
SINROptCorrelatedEMMSE_Save = zeros(nbrBSs,K,NumIter);

for Iter = 1:NumIter
    
    % Update U matrix
    [Umatrix, ShareTerm] = Compute_Umatrix(SqrtDataPowermatrix, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K);
    
    % Update W matrix
    Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTerm,Bmatrix,nbrBSs, K);
    SqrtDataPowermatrix = Compute_Rho(Umatrix,Wmatrix,Bmatrix,Csmallmatrix,DataPowerMax,nbrBSs, K);
    DataPowermatrixOptMMSE = abs(SqrtDataPowermatrix).^2;
    DataPowermatrixSave(:,:,Iter) = DataPowermatrixOptMMSE;
    for l = 1 : nbrBSs
        for k = 1:K
            Numerator = DataPowermatrixOptMMSE(l,k)*Bmatrix(l,l,k)^2;
            Denominator = Dmatrix(l,k);
            for i = [1:l-1,l+1:nbrBSs]
                Denominator = Denominator + DataPowermatrixOptMMSE(i,k)*Bmatrix(l,i,k)^2;
            end
            for i = 1: nbrBSs
                for  t = 1 :K
                    Denominator  = Denominator + DataPowermatrixOptMMSE(i,t)*Csmallmatrix(i,t,l,k);
                end
            end
            SINROptCorrelatedEMMSE(l,k) = Numerator/Denominator;
            clear Numerator Denominator;
        end
    end
    SINROptCorrelatedEMMSE_Save(:,:,Iter)= SINROptCorrelatedEMMSE;
    Asave(Iter)= sum(log2(1+SINROptCorrelatedEMMSE(:)));
end

SINROptCorrelatedEMMSE = squeeze(SINROptCorrelatedEMMSE_Save(:,:,end));
DataPowermatrixOptMMSE = squeeze(DataPowermatrixSave(:,:,end));

end

function Bmatrix = Compute_Bmatrix(tau,EstError,EstPhi,L,K)
Bmatrix = zeros(L, L, K);
for j = 1 : L
    for i = 1: L
        for k = 1:K
            Bmatrix(j,i,k) = sqrt(tau)*EstError(j,j,k)*EstError(j,i,k)*abs(trace(squeeze(EstPhi(j,k,:,:))));
        end
    end
end
end

function Csmallmatrix = Compute_Cmatrix(EstError,EstPhi,CorrelatedFadingFun,L,K)
Csmallmatrix =  zeros(L, K, L, K); %c_{jk}^it in the paper

% Compute Cmatrix
for i = 1 : L
    for t = 1 : K
        for j = 1 : L
            for k = 1 : K
                Csmallmatrix(i,t,j,k) = ((EstError(j,j,k))^2)*abs(trace(squeeze(CorrelatedFadingFun(j,i,t,:,:))*squeeze(EstPhi(j,k,:,:))));
            end
        end
    end
end
end

function Dmatrix = Compute_Dmatrix(EstError,EstPhi,L,K)
Dmatrix = zeros(L, K);
% Compute Dmatrix
for j = 1 : L
    for k = 1 : K
        Dmatrix(j,k) = ((EstError(j,j,k))^2)*abs(trace(squeeze(EstPhi(j,k,:,:))));
    end
end
end

function [Umatrix, ShareTerm] = Compute_Umatrix(SqrtDataPowermatrix, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K)
Umatrix = zeros(nbrBSs,K);
ShareTerm = zeros(nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1: K
        Numerator_ulk = SqrtDataPowermatrix(l,k)*squeeze(Bmatrix(l,l,k));
        Denominator_ulk = Dmatrix(l,k);
        for i = 1 : nbrBSs
            for t = 1 :K
                Denominator_ulk = Denominator_ulk + (SqrtDataPowermatrix(i,t)^2)*squeeze(Csmallmatrix(i,t,l,k));
            end
        end
        for i = 1 :nbrBSs
            Denominator_ulk = Denominator_ulk + (SqrtDataPowermatrix(i,k)^2)*squeeze(Bmatrix(l,i,k))^2;
        end
        Umatrix(l,k) = Numerator_ulk/Denominator_ulk;
        ShareTerm(l,k) = Denominator_ulk;
    end
end

end

function Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTerm,Bmatrix,nbrBSs, K)
Wmatrix = zeros(nbrBSs, K);
for l = 1 : nbrBSs
    for k = 1 : K
        elk = (Umatrix(l,k)^2)*ShareTerm(l,k) - 2*Umatrix(l,k)*SqrtDataPowermatrix(l,k)*Bmatrix(l,l,k) + 1;
        Wmatrix(l,k) = 1/elk;
    end
end
end

function SqrtDataPowermatrix = Compute_Rho(Umatrix,Wmatrix,Bmatrix,Csmallmatrix,DataPowerMax,nbrBSs, K)
SqrtDataPowermatrix = zeros(nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1 : K
        Numerator_rholk = Wmatrix(l,k)*Umatrix(l,k)*squeeze(Bmatrix(l,l,k));
        Denominator_rholk = 0;
        for i = 1 : nbrBSs
            for t = 1: K
                Denominator_rholk = Denominator_rholk + Wmatrix(i,t)*(Umatrix(i,t)^2)*squeeze(Csmallmatrix(l,k,i,t));
            end
        end
        for i = 1 : nbrBSs
            Denominator_rholk = Denominator_rholk + Wmatrix(i,k)*(Umatrix(i,k)^2)*squeeze(Bmatrix(i,l,k))^2;
        end
        SqrtDataPowermatrix(l,k) = min(Numerator_rholk/Denominator_rholk,sqrt(DataPowerMax));
        
    end
end

end
