function [SINROptLSFDCorrelatedEMMSE_Appr, DataPowermatrix, Asave] =  Func_OptLSFD_CorrelatedEMMSE_Appr(IntDataPowerMatrix,lossovernoise, EstError,EstPhi, CorrelatedFading, PilotPowerMatrix, DataPowerMax, nbrBSs,K,tau,NBScases, NumIter)
%This function optimizes SINR values of all users in case of element-wise
%MMSE estimation and approximated LSFD vectors
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


SINROptLSFDCorrelatedEMMSE_Appr = zeros(nbrBSs,K);% Space for SINRvalues;
Bmatrix = Compute_Bmatrix(tau,EstError,EstPhi,nbrBSs,K);
Csmallmatrix = Compute_Cmatrix(EstError,EstPhi,CorrelatedFading,nbrBSs,K);
Dmatrix = Compute_Dmatrix(EstError,EstPhi,nbrBSs,K);
SqrtDataPowermatrix =IntDataPowerMatrix;

% Compute Approximation
Bmatrix_Appr = Compute_Bmatrix_Appr(NBScases, tau,PilotPowerMatrix,EstError,lossovernoise,nbrBSs,K);
Csmallmatrix_Appr = Compute_Cmatrix_Appr(PilotPowerMatrix,lossovernoise,EstError,nbrBSs,K);
Dmatrix_Appr =Compute_Dmatrix_Appr(PilotPowerMatrix,lossovernoise,EstError,nbrBSs,K);
Amatrix_Appr = zeros(nbrBSs,nbrBSs,K);

FlagNote = 0;
for l = 1 : nbrBSs
    for k = 1:K
        Clk_Appr  = Compute_Clk(Bmatrix_Appr,Dmatrix_Appr,Csmallmatrix_Appr,SqrtDataPowermatrix.^2,nbrBSs,l,K,k, FlagNote);
        Amatrix_Appr(:,l,k) = ((Clk_Appr)\squeeze(Bmatrix_Appr(:,l,k)));
    end
end
Asave = zeros(NumIter,1);
DataPowermatrixSave = zeros(nbrBSs,K,NumIter);
SINROptLSFDCorrelatedEMMSE_ApprSave = zeros(nbrBSs,K,NumIter);
for Iter = 1:NumIter
    
    % Update U matrix
    [Umatrix, ShareTerm] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix_Appr, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K);
    
    % Update W matrix
    Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTerm,Amatrix_Appr,Bmatrix,nbrBSs, K);
    
    % Update A matrix
    Amatrix_Appr  = Compute_Amatrix_Appr(SqrtDataPowermatrix,Bmatrix_Appr, Csmallmatrix_Appr, Dmatrix_Appr, nbrBSs, K);
    SqrtDataPowermatrix = Compute_Rho(Umatrix,Wmatrix,Amatrix_Appr,Bmatrix,Csmallmatrix,DataPowerMax,nbrBSs, K);
    DataPowermatrix = SqrtDataPowermatrix.^2;
    DataPowermatrixSave(:,:,Iter) = DataPowermatrix;
    for l = 1 : nbrBSs
        for k = 1:K
            Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,DataPowermatrix,nbrBSs,l,K,k,0);
            Numerator = DataPowermatrix(l,k)*abs(squeeze(Amatrix_Appr(:,l,k))'*squeeze(Bmatrix(:,l,k)))^2;
            Denominator = squeeze(Amatrix_Appr(:,l,k))'*Clk*squeeze(Amatrix_Appr(:,l,k));
            SINROptLSFDCorrelatedEMMSE_Appr(l,k) = Numerator/Denominator;
        end
    end
    Asave(Iter)= sum(log2(1+SINROptLSFDCorrelatedEMMSE_Appr(:)));
    SINROptLSFDCorrelatedEMMSE_ApprSave(:,:,Iter)= SINROptLSFDCorrelatedEMMSE_Appr;
end

[~,IterMax]= max(Asave);
SINROptLSFDCorrelatedEMMSE_Appr = squeeze(SINROptLSFDCorrelatedEMMSE_ApprSave(:,:,IterMax));
DataPowermatrix = squeeze(DataPowermatrixSave(:,:,IterMax));

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

function Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,DataPower,L,l,K,k, FlagNote)

% compute the first part of Clk
Clk = zeros(L,L);
if FlagNote == 0
    for i = [1:l-1,l+1:L]
        Clk = Clk + DataPower(i,k)*squeeze(Bmatrix(:,i,k))*squeeze(Bmatrix(:,i,k))';
    end
else
    for i = 1 : L
        Clk = Clk + DataPower(i,k)*squeeze(Bmatrix(:,i,k))*squeeze(Bmatrix(:,i,k))';
    end
end

% compute the second part of Clk
ctemplk = zeros(L,1);
for s = 1: L
    ctemplk(s) = Dmatrix(s,k);
    for i =1 : L
        for t = 1: K
            ctemplk(s) = ctemplk(s) + DataPower(i,t)*Csmallmatrix(i,t,s,k);
        end
    end
end
Clk = Clk + diag(ctemplk);
end

function [Umatrix, ShareTerm] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K)

Umatrix = zeros(nbrBSs,K);
ShareTerm = zeros(nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1: K
        Numerator_ulk = SqrtDataPowermatrix(l,k)*squeeze(Amatrix(:,l,k))'*squeeze(Bmatrix(:,l,k));
        Denominator_ulk = (squeeze(Amatrix(:,l,k)).^2)'*Dmatrix(:,k);
        for i = 1 : nbrBSs
            for t = 1 :K
                Denominator_ulk = Denominator_ulk + (SqrtDataPowermatrix(i,t)^2)*(squeeze(Amatrix(:,l,k)).^2)'*squeeze(Csmallmatrix(i,t,:,k));
            end
        end
        for i = 1 :nbrBSs
            Denominator_ulk = Denominator_ulk + (SqrtDataPowermatrix(i,k)^2)*(squeeze(Amatrix(:,l,k))'*squeeze(Bmatrix(:,i,k)))^2;
        end
        Umatrix(l,k) = Numerator_ulk/Denominator_ulk;
        ShareTerm(l,k) = Denominator_ulk;
    end
end

end

function Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTerm,Amatrix,Bmatrix,nbrBSs, K)

Wmatrix = zeros(nbrBSs, K);
for l = 1 : nbrBSs
    for k = 1 : K
        elk = (Umatrix(l,k)^2)*ShareTerm(l,k) - 2*Umatrix(l,k)*SqrtDataPowermatrix(l,k)*abs(squeeze(Amatrix(:,l,k))'*Bmatrix(:,l,k)) + 1;
        Wmatrix(l,k) = 1/elk;
    end
end
end

function  Amatrix_Appr  = Compute_Amatrix_Appr(SqrtDataPowermatrix,Bmatrix_Appr, Csmallmatrix_Appr, Dmatrix_Appr, nbrBSs, K)

Amatrix_Appr = zeros(nbrBSs,nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1 : K
        Clk_Appr  = Compute_Clk(Bmatrix_Appr,Dmatrix_Appr,Csmallmatrix_Appr,SqrtDataPowermatrix.^2,nbrBSs,l,K,k,0);
        Amatrix_Appr(:,l,k) = ((Clk_Appr)\squeeze(Bmatrix_Appr(:,l,k)));
    end
end
end

function SqrtDataPowermatrix = Compute_Rho(Umatrix,Wmatrix,Amatrix,Bmatrix,Csmallmatrix,DataPowerMax,nbrBSs, K)

SqrtDataPowermatrix = zeros(nbrBSs,K);

for l = 1 : nbrBSs
    for k = 1 : K
        Numerator_rholk = Wmatrix(l,k)*Umatrix(l,k)*abs(squeeze(Amatrix(:,l,k))'*squeeze(Bmatrix(:,l,k)));
        Denominator_rholk = 0;
        for i = 1 : nbrBSs
            for t = 1: K
                Denominator_rholk = Denominator_rholk + Wmatrix(i,t)*(Umatrix(i,t)^2)*(squeeze(Amatrix(:,i,t)).^2)'*squeeze(Csmallmatrix(l,k,:,t));
            end
        end
        for i = 1 : nbrBSs
            Denominator_rholk = Denominator_rholk + Wmatrix(i,k)*(Umatrix(i,k)^2)*(squeeze(Amatrix(:,i,k))'*squeeze(Bmatrix(:,l,k)))^2;
        end
        SqrtDataPowermatrix(l,k) = min(Numerator_rholk/Denominator_rholk,sqrt(DataPowerMax));
        
    end
end
end

function Bmatrix_Appr = Compute_Bmatrix_Appr(NBScases, tau,PilotPower,EstError,lossovernoise,L,K)
Bmatrix_Appr = zeros(L, L, K);

% Compute Approximatation of Bmatrix
for j = 1 : L
    for i = 1: L
        for k = 1:K
            Bmatrix_Appr(j,i,k) = sqrt(NBScases*tau*PilotPower(i,k))*EstError(j,j,k)*lossovernoise(j,i,k);
        end
    end
end
end

function Csmallmatrix_Appr = Compute_Cmatrix_Appr(PilotPower,lossovernoise,EstError,L,K)
Csmallmatrix_Appr =  zeros(L, K, L, K); %c_{jk}^it in the paper

% Compute Approximation of Cmatrix
for i = 1 : L
    for t = 1 : K
        for j = 1 : L
            for k = 1 : K
                Csmallmatrix_Appr(i,t,j,k) = (EstError(j,j,k))*sqrt(PilotPower(j,k))*lossovernoise(j,j,k)*lossovernoise(j,i,t);
            end
        end
    end
end
end

function   Dmatrix_Appr = Compute_Dmatrix_Appr(PilotPower,lossovernoise,EstError,L,K)
Dmatrix_Appr = zeros(L, K);

% Compute Approximation of Dmatrix
for j = 1 : L
    for k = 1 : K
        Dmatrix_Appr(j,k) = (EstError(j,j,k))*sqrt(PilotPower(j,k))*lossovernoise(j,j,k);
    end
end

end
