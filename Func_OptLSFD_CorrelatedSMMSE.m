function [SINROptLSFDCorrelatedSMMSE, DataPowermatrix, Asave] =  Func_OptLSFD_CorrelatedSMMSE(IntDataPowerMatrix,InveseEstPhi, CorrelatedFading, DataPowerMax, PilotPowerMatrix,nbrBSs,K,tau, NumIter)
%This function computes SINR values of all users in case of standard (real)
%MMSE estimation and optimal LSFD vectors 
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


SINROptLSFDCorrelatedSMMSE = zeros(nbrBSs,K);% Space for SINRvalues;
Bmatrix = Compute_Bmatrix(tau,PilotPowerMatrix,InveseEstPhi,CorrelatedFading, nbrBSs, K);
Csmallmatrix = Compute_Cmatrix(PilotPowerMatrix,InveseEstPhi,CorrelatedFading,nbrBSs,K);
Dmatrix = Compute_Dmatrix(PilotPowerMatrix,InveseEstPhi,CorrelatedFading,nbrBSs,K);
SqrtDataPowermatrix = IntDataPowerMatrix;
Amatrix = zeros(nbrBSs,nbrBSs,K);


FlagNote = 0;
for l = 1 : nbrBSs
    for k = 1:K
        Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,SqrtDataPowermatrix.^2,nbrBSs,l,K,k, FlagNote);
        Amatrix(:,l,k) = ((Clk)\squeeze(Bmatrix(:,l,k)));
    end
end
FlagNote = 1;

Asave = zeros(NumIter,1);
DataPowermatrixSave = zeros(nbrBSs,K,NumIter);
SINROptLSFDCorrelatedSMMSE_Save = zeros(nbrBSs,K,NumIter);
for Iter = 1:NumIter
    
    % Update U matrix
    [Umatrix, ShareTermForE, ShareTermForA] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K);
    
    % Update W matrix
    Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTermForE,Amatrix,Bmatrix,nbrBSs, K);
    
    % Update A matrix
    Amatrix  = Compute_Amatrix(SqrtDataPowermatrix,Bmatrix, Csmallmatrix, Dmatrix, ShareTermForA, nbrBSs, K,FlagNote);
    SqrtDataPowermatrix = Compute_Rho(Umatrix,Wmatrix,Amatrix,Bmatrix,Csmallmatrix,DataPowerMax,nbrBSs, K);
    DataPowermatrix = SqrtDataPowermatrix.^2;
    DataPowermatrixSave(:,:,Iter) = DataPowermatrix;
    for l = 1 : nbrBSs
        for k = 1:K
            Clk  = Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,DataPowermatrix,nbrBSs,l,K,k, 0);
            SINROptLSFDCorrelatedSMMSE(l,k) = DataPowermatrix(l,k)*squeeze(Bmatrix(:,l,k))'*((Clk)\squeeze(Bmatrix(:,l,k)));
        end
    end
    Asave(Iter)= sum(log2(1+SINROptLSFDCorrelatedSMMSE(:)));
    SINROptLSFDCorrelatedSMMSE_Save(:,:,Iter)= SINROptLSFDCorrelatedSMMSE;
end

%[~,IterMax]= max(Asave);
SINROptLSFDCorrelatedSMMSE = squeeze(SINROptLSFDCorrelatedSMMSE_Save(:,:,end));
DataPowermatrix = squeeze(DataPowermatrixSave(:,:,end));

end

function Bmatrix = Compute_Bmatrix(tau,PilotPower,InveseEstPhiFun,CorrelatedFadingFun,L,K)

Bmatrix = zeros(L, L, K);
for j = 1 : L
    for i = 1: L
        for k = 1:K
            Bmatrix(j,i,k) = sqrt(tau*PilotPower(i,k)*PilotPower(j,k))*abs(trace(squeeze(InveseEstPhiFun(j,k,:,:))*squeeze(CorrelatedFadingFun(j,j,k,:,:))*squeeze(CorrelatedFadingFun(j,i,k,:,:))));
        end
    end
end
end

function Csmallmatrix = Compute_Cmatrix(PilotPower,InveseEstPhiFun,CorrelatedFadingFun,L,K)

Csmallmatrix =  zeros(L, K, L, K); %c_{jk}^it in the paper
% Compute Cmatrix
for i = 1 : L
    for t = 1 : K
        for j = 1 : L
            for k = 1 : K
                %m = PilotPowerMatrix(i,k)*trace(squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(InveseEstPhi(j,k,:,:))*squeeze(CorrelatedFading(j,j,k,:,:))*squeeze(CorrelatedFading(j,i,t,:,:)))
                Csmallmatrix(i,t,j,k) = PilotPower(i,k)*abs(trace(squeeze(CorrelatedFadingFun(j,j,k,:,:))*squeeze(InveseEstPhiFun(j,k,:,:))*squeeze(CorrelatedFadingFun(j,j,k,:,:))*squeeze(CorrelatedFadingFun(j,i,t,:,:))));
            end
        end
    end
end
end

function Dmatrix = Compute_Dmatrix(PilotPower,InveseEstPhiFun,CorrelatedFadingFun,L,K)

Dmatrix = zeros(L, K);
% Compute Dmatrix
for j = 1 : L
    for k = 1 : K
        Dmatrix(j,k) = PilotPower(j,k)*abs(trace(squeeze(InveseEstPhiFun(j,k,:,:))*squeeze(CorrelatedFadingFun(j,j,k,:,:))*squeeze(CorrelatedFadingFun(j,j,k,:,:))));
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

function [Umatrix, ShareTermForE, ShareTermForA] = Compute_Umatrix(SqrtDataPowermatrix,Amatrix, Bmatrix,Csmallmatrix,Dmatrix, nbrBSs, K)

Umatrix = zeros(nbrBSs,K);
ShareTermForE = zeros(nbrBSs,K);
ShareTermForA = zeros(nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1: K
        Numerator_ulk = SqrtDataPowermatrix(l,k)*squeeze(Amatrix(:,l,k))'*squeeze(Bmatrix(:,l,k));
        Numerator_ulkForA = squeeze(Amatrix(:,l,k))'*squeeze(Bmatrix(:,l,k));
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
        ShareTermForE(l,k) = Denominator_ulk;
        ShareTermForA(l,k) = Numerator_ulkForA/Denominator_ulk;
    end
end

end

function Wmatrix  = Compute_Wmatrix(SqrtDataPowermatrix,Umatrix,ShareTermForE,Amatrix,Bmatrix,nbrBSs, K)

Wmatrix = zeros(nbrBSs, K);
for l = 1 : nbrBSs
    for k = 1 : K
        elk = (Umatrix(l,k)^2)*ShareTermForE(l,k) - 2*Umatrix(l,k)*SqrtDataPowermatrix(l,k)*abs(squeeze(Amatrix(:,l,k))'*Bmatrix(:,l,k)) + 1;
        Wmatrix(l,k) = 1/elk;
    end
end
end

function  Amatrix  = Compute_Amatrix(SqrtDataPowermatrix,Bmatrix, Csmallmatrix, Dmatrix, ShareTermForA, nbrBSs, K,FlagNote)

Amatrix = zeros(nbrBSs,nbrBSs,K);
for l = 1 : nbrBSs
    for k = 1 : K
        Clk =  Compute_Clk(Bmatrix,Dmatrix,Csmallmatrix,SqrtDataPowermatrix.^2,nbrBSs,l,K,k, FlagNote);
        Amatrix(:,l,k) =(1/ShareTermForA(l,k))*(Clk\squeeze(Bmatrix(:,l,k)));
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
