function  [SINRMRCLSFDmatrix, SINRMRCmatrix, SINRRZFLSFDmatrix,SINRRZFmatrix] = Func_RateDiffLinear(nbrSmallScaleRealizations, tau,K,nbrBSs,NBScases, PowerSymbol, InveseEstPhi,Cmatrix, ChannelRealizations,NoiseMatrix,CorrelatedFading)
%This function computes SINR values of all users in case of standard (real)
%MMSE estimation and MRC/RZF precoding vectors in the first layer
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


Gmatrix = zeros(nbrBSs,K,NBScases,nbrSmallScaleRealizations);
for l = 1:nbrBSs
    for k = 1: K
        hlksum = squeeze(sum(squeeze(ChannelRealizations(l,:,k,:,:)),1));
        Noiselk = squeeze(NoiseMatrix(l,l,k,:,:));
        for m = 1: nbrSmallScaleRealizations
            Gmatrix(l,k,:,m) = sqrt(PowerSymbol)*squeeze(CorrelatedFading(l,l,k,:,:))*squeeze(InveseEstPhi(l,k,:,:))*(tau*sqrt(PowerSymbol)*hlksum(:,m) + Noiselk(:,m));
        end
        
    end
end

% Define MR combining vectors
VMRCmatrix = Gmatrix;

% Define RZF combining vectors
Idmatrix = eye(K);
PjInversematrix = Idmatrix/PowerSymbol;
VRZFmatrix = zeros(size(Gmatrix));
for l=1:nbrBSs
    for m=1:nbrSmallScaleRealizations
        Tempmatrix  = (squeeze(Gmatrix(l,:,:,m))).';
        VRZFmatrix(l,:,:,m) = (Tempmatrix/((Tempmatrix')*Tempmatrix + PjInversematrix)).';
    end
end


[SINRMRCLSFDmatrix,SINRMRCmatrix] = Func_ComputeSINRmatrix(VMRCmatrix, ChannelRealizations, nbrSmallScaleRealizations, PowerSymbol,nbrBSs,K);
[SINRRZFLSFDmatrix,SINRRZFmatrix] = Func_ComputeSINRmatrix(VRZFmatrix, ChannelRealizations, nbrSmallScaleRealizations, PowerSymbol,nbrBSs,K);

end

function Bmatrix = Func_ComputeBmatrix(Vmatrix,ChannelRealizations,K,nbrBSs)

Bmatrix = zeros(nbrBSs,K,nbrBSs);
for l=1:nbrBSs
    for k = 1:K
        for i=1: nbrBSs
            hilk = squeeze(ChannelRealizations(i,l,k,:,:));
            vik = squeeze(Vmatrix(i,k,:,:));
            Bmatrix(l,k,i) = mean(sum(conj(vik).*hilk,1));
        end
    end
end
end

function Conematrix = Func_ComputeConematrix(Bmatrix,PowerSymbol,K,nbrBSs)

Conematrix = zeros(nbrBSs,K,nbrBSs,nbrBSs);
for l = 1: nbrBSs
    for k = 1: K
        TempVal = zeros(nbrBSs,nbrBSs);
        for i = [1:l-1, l+1:nbrBSs]
            Btemp = squeeze(Bmatrix(i,k,:));
            TempVal = TempVal + PowerSymbol*Btemp*(Btemp');
        end
        Conematrix(l,k,:,:) = TempVal;
    end
end

end

function Ctwomatrix = Func_ComputeCtwomatrix(Vmatrix,ChannelRealizations,PowerSymbol, K,nbrBSs,nbrSmallScaleRealizations)

Ctwomatrix =  zeros(nbrBSs,K,nbrBSs,nbrBSs);
for l = 1:nbrBSs
    for k = 1:K
        VectiTemp = zeros(nbrBSs,nbrBSs);
        for i = 1: nbrBSs
            VectjTemp = zeros(nbrBSs,nbrSmallScaleRealizations);
            for j=1:nbrBSs
                VvecTemp = squeeze(Vmatrix(j,k,:,:));
                Tempval = conj(VvecTemp).*squeeze(ChannelRealizations(j,i,k,:,:));
                TempValv1 = sum(Tempval,1);
                TempValv2 = abs(mean(sum(Tempval,1)));
                VectjTemp(j,:) = TempValv1 - TempValv2;
            end
            VectiAvg = zeros(nbrBSs,nbrBSs);
            for m= 1: nbrSmallScaleRealizations
                VectiAvg = VectiAvg + VectjTemp(:,m)*(VectjTemp(:,m)');
            end
            VectiAvg = VectiAvg/nbrSmallScaleRealizations;
            VectiTemp = VectiTemp + PowerSymbol*VectiAvg;
        end
        Ctwomatrix(l,k,:,:) =VectiTemp;
    end
end
end

function Cthreematrix = Func_ComputeCthreematrix(Vmatrix,ChannelRealizations,PowerSymbol, K,nbrBSs)

Cthreematrix =  zeros(nbrBSs,K,nbrBSs,nbrBSs);
for l = 1:nbrBSs
    for k = 1:K
        VectiTemp = zeros(nbrBSs,nbrBSs);
        for i = 1: nbrBSs
            for t = [1:k-1, k+1:K]
                VectjTemp = zeros(nbrBSs,1);
                for j=1:nbrBSs
                    VvecTemp = squeeze(Vmatrix(j,k,:,:));
                    ChannelTemp = squeeze(ChannelRealizations(j,i,t,:,:));
                    Tempval = conj(VvecTemp).*ChannelTemp;
                    VectjTemp(j) = PowerSymbol*mean(abs(sum(Tempval,1)).^2);
                end
                VectiTemp = VectiTemp +  diag(VectjTemp);
            end
        end
        Cthreematrix(l,k,:,:) =VectiTemp;
    end
end
end

function Cfourmatrix = Func_ComputeCfourmatrix(Vmatrix, K,nbrBSs)

Cfourmatrix = zeros(K,nbrBSs,nbrBSs);
for k = 1:K
    Vvalue = zeros(nbrBSs,1);
    for j = 1 : nbrBSs
        vjk = squeeze(Vmatrix(j,k,:,:));
        Vvalue(j) = mean(sum(conj(vjk).*vjk,1));
    end
    Cfourmatrix(k,:,:) = diag(Vvalue);
    
end
end

function [SINRLFSDmatrix, SINRmatrix] = Func_ComputeSINRmatrix(Vmatrix, ChannelRealizations, nbrSmallScaleRealizations, PowerSymbol,nbrBSs,K)

Bmatrix = Func_ComputeBmatrix(Vmatrix,ChannelRealizations,K,nbrBSs);
Conematrix = Func_ComputeConematrix(Bmatrix,PowerSymbol,K,nbrBSs);
Ctwomatrix = Func_ComputeCtwomatrix(Vmatrix,ChannelRealizations,PowerSymbol, K,nbrBSs,nbrSmallScaleRealizations);
Cthreematrix = Func_ComputeCthreematrix(Vmatrix,ChannelRealizations,PowerSymbol, K,nbrBSs);
Cfourmatrix = Func_ComputeCfourmatrix(Vmatrix, K,nbrBSs);
SINRLFSDmatrix = zeros(nbrBSs,K);
SINRmatrix = zeros(nbrBSs,K);

for l = 1 : nbrBSs
    for k = 1:K
        Bmatrixlk = squeeze(Bmatrix(l,k,:));
        Conematrixlk = squeeze(Conematrix(l,k,:,:));
        Ctwomatrixlk = squeeze(Ctwomatrix(l,k,:,:));
        Cthreematrixlk = squeeze(Cthreematrix(l,k,:,:));
        Cfourmatrixlk = squeeze(Cfourmatrix(k,:,:));
        Clk = Conematrixlk + Ctwomatrixlk + Cthreematrixlk + Cfourmatrixlk;
        SINRLFSDmatrix(l,k) = PowerSymbol*((Bmatrixlk')/Clk)*Bmatrixlk;
        
        alk = zeros(nbrBSs,1);
        alk(l)=1;
        SINRmatrix(l,k) = PowerSymbol*((abs((alk')*Bmatrixlk))^2)/((alk')*Clk*alk);
        
    end
    
end

SINRLFSDmatrix = abs(SINRLFSDmatrix);
SINRmatrix = abs(SINRmatrix);

end