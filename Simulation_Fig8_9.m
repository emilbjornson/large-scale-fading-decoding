%This Matlab script can be used to reproduce Figures 8,9 in the paper:
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


clc; clear; close all;
rng('shuffle'); %Initiate the random number generators with a random seed

nbrOfRealizations =300; % Number of channels realizations
NBScases = 200; % Number of antannas per BS
squareLength = 1.00; % Define the coverage area (as a square)
nbrBSs = 4; % Number of BS in the area
nbrUsers = 2:2:10; % Possible Number of user per BS% Constraint power
rho_BS = 1;
Bandwidth = 20e6; % Bandwidth in Hz
noiseFigure = 5;%dB
noiseFloordBm = -174+10*log10(Bandwidth) + noiseFigure;
PowerSymbol= 0.2e3;

nbrBSsPerDim = sqrt(nbrBSs); % BSs per dimension
interSiteDistance = squareLength/nbrBSsPerDim; %Distance between BSs in vertical/horizontal direction
% Put out BSs on the grid
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:)+1i*locationsGridVertical(:);
% Compute the exact dimension of the square where the users are located
maxDistance = interSiteDistance;

%Compute alterantive BS locations by using wrap arround
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions, [1 length(wrapLocations)]) + repmat(wrapLocations, [nbrBSs 1]);

%Channel Profile
pentrationlossdB = 20; %dB
ShadowFadingBS = 7;
%Vector for correlated Rayleigh channels
ArrayVector = 0:1:(NBScases-1);
CorreFactor = 0.5;
Im = eye(NBScases);
NumIter = 500; %Number of Iterations for Optimization

SumRateLSFDCorrelatedSMMSE= zeros(nbrOfRealizations,length(nbrUsers));
SumRateCorrelatedSMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateLSFDCorrelatedSMMSE_Appr = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptLSFDCorrelatedSMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptCorrelatedSMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptLSFDCorrelatedSMMSE_Appr = zeros(nbrOfRealizations,length(nbrUsers));
SumRateLSFDCorrelatedEMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateCorrelatedEMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateLSFDCorrelatedEMMSE_Appr = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptLSFDCorrelatedEMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptCorrelatedEMMSE = zeros(nbrOfRealizations,length(nbrUsers));
SumRateOptLSFDCorrelatedEMMSE_Appr = zeros(nbrOfRealizations,length(nbrUsers));
ObjFunLFSD_SMMSESave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
ObjFun_SMMSESave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
ObjFunLFSD_EMMSESave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
ObjFun_EMMSESave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
ObjFunLFSD_SMMSE_ApprSave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
ObjFunLFSD_EMMSE_ApprSave = zeros(length(nbrUsers),NumIter,nbrOfRealizations);
for Iterk = 1: length(nbrUsers)
    fprintf('%d users per cell \n', nbrUsers(Iterk));
    K = nbrUsers(Iterk);
    tau= K;
    DataPowerMatix = PowerSymbol*ones(nbrBSs,K);
    PilotPowerMatrix = PowerSymbol*ones(nbrBSs,K);
    for iter = 1: nbrOfRealizations
        %iter
        %Prepare to put out UEs in the cells
        UEpositions = zeros(K,nbrBSs);
        perBS = zeros(nbrBSs,1);

        AngleNet = zeros(nbrBSs, nbrBSs, K);
        lossovernoise = zeros(nbrBSs, nbrBSs, K);
        EstError = zeros(nbrBSs, nbrBSs, K);
        CorrelatedFading  = zeros(nbrBSs, nbrBSs,K,NBScases,NBScases);
        EstPhi = zeros(nbrBSs, K, NBScases,NBScases);
        InveseEstPhi = zeros(nbrBSs, K, NBScases,NBScases);
        Amatrix = zeros(nbrBSs, nbrBSs, K);

        for l=1:nbrBSs
            while min(perBS(l)) < K
                for k=1: K
                    posX = rand(1)*maxDistance + real(BSpositions(l)) - maxDistance/2;
                    posY = rand(1)*maxDistance + imag(BSpositions(l)) - maxDistance/2;
                    PositionTemp = posX + 1i*posY;
                    while abs(BSpositions(l) - PositionTemp) < 0.035
                        posX = rand(1)*maxDistance + real(BSpositions(l)) - maxDistance/2;
                        posY = rand(1)*maxDistance + imag(BSpositions(l)) - maxDistance/2;
                        PositionTemp = posX + 1i*posY;
                    end
                    UEpositions(k,l) = posX + 1i*posY;
                end
                perBS(l) = K;
            end
        end
        % Generate Large Scale Fading and Angles of users to BSs
        for l=1:nbrBSs
            BSSelectedPosK = zeros(nbrBSs,1);
            for k = 1: K
                [distancesSquaredBSj,Position] = min(abs(repmat(UEpositions(k,l),size(BSpositionsWrapped))- BSpositionsWrapped),[],2);
                % Compute the angles of user and BSs
                for BSSelectedIter = 1: nbrBSs
                    BSSelectedPosK(BSSelectedIter) = BSpositionsWrapped(BSSelectedIter, Position(BSSelectedIter));
                end
                VectorK = UEpositions(k,l) - BSSelectedPosK;
                AngleK = angle(VectorK);
                AngleNet(:,l,k) = AngleK;
                loss = 128.1 + 37.6*log10(distancesSquaredBSj) + ShadowFadingBS*randn(nbrBSs,1) + noiseFloordBm + pentrationlossdB;
                 % Check to ensure that home BS having largest BS
                 [Value, Location] = min(loss);
                 while Location ~= l
                     Temploss = 128.1 + 37.6*log10(distancesSquaredBSj(l)) + ShadowFadingBS*randn(1) + noiseFloordBm + pentrationlossdB;
                     loss(l) = Temploss;
                     [Value, Location] = min(loss);
                 end
                 lossovernoise(:,l,k) = rho_BS./(10.^(loss/10));
            end
        end
        % Compute autocorrelation matrices 
        for l = 1 : nbrBSs
            for i = 1 : nbrBSs
                for k = 1 : K
                    % Compute autocorrelation matrix
                    CorreVec = lossovernoise(l,i,k)*(CorreFactor*exp(1i*AngleNet(l,i,k))).^ArrayVector;
                    CorrelatedFading(l,i,k,:,:) = toeplitz(CorreVec);


                    % Compute MMSE estimation factors
                    NumFactor = (sqrt(PilotPowerMatrix(i,k))*lossovernoise(l,i,k));
                    DenFactor =  1 + sum(tau*PilotPowerMatrix(:,k).*squeeze(lossovernoise(l,:,k))');
                    EstError(l,i,k) = NumFactor/DenFactor;
                end
            end
        end

        for l = 1 : nbrBSs
            for k = 1 : K
                 PhiTemp = Im;
                 for s = 1 : nbrBSs
                      PhiTemp = PhiTemp + tau*PilotPowerMatrix(s,k)*squeeze(CorrelatedFading(l,s,k,:,:));
                 end
                 EstPhi(l, k,:,:) = PhiTemp;
                 InveseEstPhi(l,k,:,:) = inv(PhiTemp);
            end
        end
       IntDataPowerMatrix = sqrt(PowerSymbol)*rand(nbrBSs,K);
       % Estimate SINRs with real MMSE
       SINRLSFDCorrelatedSMMSE = Func_LSFD_CorrelatedSMMSE(InveseEstPhi, CorrelatedFading,DataPowerMatix, PilotPowerMatrix,nbrBSs,K,tau);
       SINRCorrelatedSMMSE = Func_CorrelatedSMMSE(InveseEstPhi, CorrelatedFading, DataPowerMatix, PilotPowerMatrix,nbrBSs,K,tau);
       SINRLSFDCorrelatedSMMSE_Appr = Func_LSFD_CorrelatedSMMSE_Appr(lossovernoise,EstError, InveseEstPhi, CorrelatedFading, DataPowerMatix, PilotPowerMatrix,nbrBSs,K,tau, NBScases);
       [SINROptLSFDCorrelatedSMMSE, OptDataPowermatrixLFSD_SMMSE, ObjFunLFSD_SMMSE] = Func_OptLSFD_CorrelatedSMMSE(IntDataPowerMatrix,InveseEstPhi, CorrelatedFading, PowerSymbol, PilotPowerMatrix,nbrBSs,K,tau, NumIter);
       [SINROptCorrelatedSMMSE, OptDataPowermatrix_SMMSE, ObjFun_SMMSE] = Func_Opt_CorrelatedSMMSE(IntDataPowerMatrix,InveseEstPhi, CorrelatedFading, PowerSymbol, PilotPowerMatrix,nbrBSs,K,tau, NumIter);
       [SINROptLSFDCorrelatedSMMSE_Appr, OptDataPowermatrixLFSD_SMMSE_Appr, ObjFunLFSD_SMMSE_Appr] =  Func_OptLSFD_CorrelatedSMMSE_Appr(IntDataPowerMatrix,lossovernoise, EstError, InveseEstPhi, CorrelatedFading, PowerSymbol, PilotPowerMatrix,nbrBSs,K,tau, NBScases,NumIter);
       % Estimate SINRs with element-wise MMSE
       SINRLSFDCorrelatedEMMSE = Func_LSFD_CorrelatedEMMSE(EstError,EstPhi, CorrelatedFading,DataPowerMatix, nbrBSs,K,tau);
       SINRCorrelatedEMMSE = Func_CorrelatedEMMSE(EstError, EstPhi,CorrelatedFading, DataPowerMatix, nbrBSs,K,tau);
       SINRLSFDCorrelatedEMMSE_Appr = Func_LSFD_CorrelatedEMMSE_Appr(lossovernoise, EstError,EstPhi, CorrelatedFading,DataPowerMatix, PilotPowerMatrix, nbrBSs,K,tau, NBScases);
       [SINROptLSFDCorrelatedEMMSE, OptDataPowermatrixLFSD_EMMSE, ObjFunLFSD_EMMSE] = Func_OptLSFD_CorrelatedEMMSE(IntDataPowerMatrix,EstError,EstPhi, CorrelatedFading, PowerSymbol, nbrBSs,K,tau, NumIter);
       [SINROptCorrelatedEMMSE, OptDataPowermatrix_EMMSE, ObjFun_EMMSE] = Func_Opt_CorrelatedEMMSE(IntDataPowerMatrix,EstError,EstPhi, CorrelatedFading, PowerSymbol,nbrBSs,K,tau, NumIter);
       [SINROptLSFDCorrelatedEMMSE_Appr, OptDataPowermatrixLFSD_EMMSE_Appr, ObjFunLFSD_EMMSE_Appr] = Func_OptLSFD_CorrelatedEMMSE_Appr(IntDataPowerMatrix,lossovernoise, EstError,EstPhi, CorrelatedFading, PilotPowerMatrix, PowerSymbol, nbrBSs,K,tau,NBScases, NumIter);

       % Save Sum rate of all methods
       SumRateLSFDCorrelatedSMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINRLSFDCorrelatedSMMSE(:)));
       SumRateCorrelatedSMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINRCorrelatedSMMSE(:)));
       SumRateLSFDCorrelatedSMMSE_Appr(iter,Iterk) =(1-K/200)*sum(log2(1+SINRLSFDCorrelatedSMMSE_Appr(:)));
       SumRateOptLSFDCorrelatedSMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptLSFDCorrelatedSMMSE(:)));
       SumRateOptCorrelatedSMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptCorrelatedSMMSE(:)));
       SumRateOptLSFDCorrelatedSMMSE_Appr(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptLSFDCorrelatedSMMSE_Appr(:)));
       SumRateLSFDCorrelatedEMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINRLSFDCorrelatedEMMSE(:)));
       SumRateCorrelatedEMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINRCorrelatedEMMSE(:)));
       SumRateLSFDCorrelatedEMMSE_Appr(iter,Iterk) = (1-K/200)*sum(log2(1+SINRLSFDCorrelatedEMMSE_Appr(:)));
       SumRateOptLSFDCorrelatedEMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptLSFDCorrelatedEMMSE(:)));
       SumRateOptCorrelatedEMMSE(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptCorrelatedEMMSE(:)));
       SumRateOptLSFDCorrelatedEMMSE_Appr(iter,Iterk) = (1-K/200)*sum(log2(1+SINROptLSFDCorrelatedEMMSE_Appr(:)));
       ObjFunLFSD_SMMSESave(Iterk,:,iter) = (1-K/200)*ObjFunLFSD_SMMSE;
       ObjFun_SMMSESave(Iterk,:,iter) = (1-K/200)*ObjFun_SMMSE;
       ObjFunLFSD_EMMSESave(Iterk,:,iter) = (1-K/200)*ObjFunLFSD_EMMSE;
       ObjFun_EMMSESave(Iterk,:,iter) = (1-K/200)*ObjFun_EMMSE;
       ObjFunLFSD_SMMSE_ApprSave(Iterk,:,iter) = (1-K/200)*ObjFunLFSD_SMMSE_Appr;
       ObjFunLFSD_EMMSE_ApprSave(Iterk,:,iter) = (1-K/200)*ObjFunLFSD_EMMSE_Appr;
    end
    clear CorrelatedFading  EstPhi  InveseEstPhi;
end
ResultsFigs89.MeanSumRateLSFDCorrelatedSMMSE = mean(SumRateLSFDCorrelatedSMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateCorrelatedSMMSE = mean(SumRateCorrelatedSMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateLSFDCorrelatedSMMSE_Appr = mean(SumRateLSFDCorrelatedSMMSE_Appr,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptLSFDCorrelatedSMMSE =mean(SumRateOptLSFDCorrelatedSMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptCorrelatedSMMSE = mean(SumRateOptCorrelatedSMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptLSFDCorrelatedSMMSE_Appr =mean(SumRateOptLSFDCorrelatedSMMSE_Appr,1)/nbrBSs;
ResultsFigs89.MeanSumRateLSFDCorrelatedEMMSE = mean(SumRateLSFDCorrelatedEMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateCorrelatedEMMSE = mean(SumRateCorrelatedEMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateLSFDCorrelatedEMMSE_Appr = mean(SumRateLSFDCorrelatedEMMSE_Appr,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptLSFDCorrelatedEMMSE = mean(SumRateOptLSFDCorrelatedEMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptCorrelatedEMMSE = mean(SumRateOptCorrelatedEMMSE,1)/nbrBSs;
ResultsFigs89.MeanSumRateOptLSFDCorrelatedEMMSE_Appr = mean(SumRateOptLSFDCorrelatedEMMSE_Appr,1)/nbrBSs;
save Results_Figs89.mat ResultsFigs89 nbrUsers;

%=============Plot Fig 8 =============================
figure; box on; % With MMSE estimator
plot(nbrUsers,ResultsFigs89.MeanSumRateOptLSFDCorrelatedSMMSE,'b-');
hold on
plot(nbrUsers,ResultsFigs89.MeanSumRateOptLSFDCorrelatedSMMSE_Appr,'b-*');
plot(nbrUsers,ResultsFigs89.MeanSumRateOptCorrelatedSMMSE,'b-s');
plot(nbrUsers,ResultsFigs89.MeanSumRateLSFDCorrelatedSMMSE,'b-');
plot(nbrUsers,ResultsFigs89.MeanSumRateLSFDCorrelatedSMMSE_Appr,'b-*');
plot(nbrUsers,ResultsFigs89.MeanSumRateCorrelatedSMMSE,'b-s');
xlabel('Number of antennas per base station');
ylabel('Sum SE per cell [b/s/Hz]');
legend('Two-layer decoding','Two-layer decoding, Approx.', 'Single-layer decoding');
hold off
grid on;
%=============Plot Fig 9===============================
figure; box on; % With EMMSE estimator
plot(nbrUsers,ResultsFigs89.MeanSumRateOptLSFDCorrelatedEMMSE,'b-');
hold on
plot(nbrUsers,ResultsFigs89.MeanSumRateOptLSFDCorrelatedEMMSE_Appr,'b-*');
plot(nbrUsers,ResultsFigs89.MeanSumRateOptCorrelatedEMMSE,'b-s');
plot(nbrUsers,ResultsFigs89.MeanSumRateLSFDCorrelatedEMMSE,'b-');
plot(nbrUsers,ResultsFigs89.MeanSumRateLSFDCorrelatedEMMSE_Appr,'b-*');
plot(nbrUsers,ResultsFigs89.MeanSumRateCorrelatedEMMSE,'b-s');
xlabel('Number of users per cell');
ylabel('Sum SE per cell [b/s/Hz]');
legend('Two-layer decoding','Two-layer decoding, Approx.', 'Single-layer decoding');
hold off
grid on