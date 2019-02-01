%This Matlab script can be used to reproduce Figure 10 in the paper:
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

nbrOfRealizations = 3000; % Number of channels realizations
nbrSmallScaleRealizations = 1000;
NBScases = 200; % Number of antannas per BS
squareLength = 1.00; % Define the coverage area (as a square)
nbrBSs = 4; % Number of BS in the area
K =5; % Number of user per BS% Constraint power
tau= K;
rho_BS = 1;
Bandwidth = 20e6; % Bandwidth in Hz
noiseFigure = 5;%dB
noiseFloordBm = -174+10*log10(Bandwidth) + noiseFigure;
PowerSymbol= 0.2e3;
DataPowerMatix = PowerSymbol*ones(nbrBSs,K);
PilotPowerMatrix = PowerSymbol*ones(nbrBSs,K);
nbrBSsPerDim = sqrt(nbrBSs); % BSs per dimension
interSiteDistance = squareLength/nbrBSsPerDim; %Distance between BSs in vertical/horizontal direction
% Put out BSs on the grid
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:)+1i*locationsGridVertical(:);

% Compute the exact dimension of the square where the users are located
maxDistance = interSiteDistance;

%Compute alterantive BS locations by using wrap around
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

lossovernoiseSave = zeros(nbrBSs, nbrBSs,K, nbrOfRealizations); % Store: Large Scale Fading
RateLSFDMRCMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);
RateMRCMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);
RateLSFDRZFMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);
RateRZFMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);
RateLSFDSMMSEMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);
RateSMMSEMonteCarlosave = zeros(nbrBSs,K,nbrOfRealizations);


%% Run Monte Carlo realizations
for iter = 1: nbrOfRealizations
    fprintf('Iteration %d of %d iterations \n', iter, nbrOfRealizations);
    
    % Small Scale fading
    ChannelRealizations = zeros(nbrBSs, nbrBSs,K, NBScases,nbrSmallScaleRealizations);
    NoiseMatrix = sqrt(tau/2)*(randn(size(ChannelRealizations)) + 1i*randn(size(ChannelRealizations)));
    
    %Prepare to put out UEs in the cells
    UEpositions = zeros(K,nbrBSs);
    perBS = zeros(nbrBSs,1);
    
    AngleNet = zeros(nbrBSs, nbrBSs, K);
    lossovernoise = zeros(nbrBSs, nbrBSs, K);
    EstError = zeros(nbrBSs, nbrBSs, K);
    CorrelatedFading  = zeros(nbrBSs, nbrBSs,K,NBScases,NBScases);
    Cmatrix  = zeros(nbrBSs, K,NBScases,NBScases);
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
            AngleK = angle(VectorK); % Compute angles
            AngleNet(:,l,k) = AngleK;
            loss = 128.1 + 37.6*log10(distancesSquaredBSj) + ShadowFadingBS*randn(nbrBSs,1) + noiseFloordBm + pentrationlossdB;
            % Check to ensure that home BS having largest BS
            [Value, Location] = min(loss);
            while Location ~= l
                Temploss = 128.1 + 37.6*log10(distancesSquaredBSj(l)) + ShadowFadingBS*randn(1) + noiseFloordBm + pentrationlossdB;
                loss(l) = Temploss;
                [Value, Location] = min(loss);
            end
            lossovernoiseSave(:,l,k,iter) = rho_BS./(10.^(loss/10));
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
    
    % Compute matrices for MMSE estimation
    for l = 1 : nbrBSs
        for k = 1 : K
            PhiTemp = Im;
            for s = 1 : nbrBSs
                PhiTemp = PhiTemp + tau*PilotPowerMatrix(s,k)*squeeze(CorrelatedFading(l,s,k,:,:));
            end
            EstPhi(l, k,:,:) = PhiTemp;
            InPhiTemp = inv(PhiTemp);
            InveseEstPhi(l,k,:,:) = InPhiTemp;
            Cmatrix(l,k,:,:) = squeeze(CorrelatedFading(l,l,k,:,:)) - PilotPowerMatrix(l,k)*tau*squeeze(CorrelatedFading(l,l,k,:,:))*InPhiTemp*squeeze(CorrelatedFading(l,l,k,:,:));
        end
    end
    
    for l = 1 : nbrBSs
        for i = 1 : nbrBSs
            for k = 1 : K
                CorrelatedFadingTemp = sqrtm(squeeze(CorrelatedFading(l,i,k,:,:)));
                for IterSmallScale = 1:nbrSmallScaleRealizations
                    SmallFadinglik = (randn(NBScases,1) + 1i*randn(NBScases,1))/sqrt(2);
                    ChannelRealizations(l,i,k,:, IterSmallScale) = CorrelatedFadingTemp*SmallFadinglik;
                    
                end
            end
        end
    end
    [SINRLSFDMRCMonteCarlo, SINRMRCMonteCarlo, SINRLSFDRZFMonteCarlo, SINRRZFMonteCarlo] = Func_RateDiffLinear(nbrSmallScaleRealizations, tau,K,nbrBSs,NBScases, PowerSymbol, InveseEstPhi,Cmatrix, ChannelRealizations,NoiseMatrix,CorrelatedFading);
    
    RateLSFDMRCMonteCarlosave(:,:,iter) = (1 -tau/200)*log2(1+SINRLSFDMRCMonteCarlo);
    RateMRCMonteCarlosave(:,:,iter) = (1 -tau/200)*log2(1+SINRMRCMonteCarlo);
    RateLSFDRZFMonteCarlosave(:,:,iter) = (1 -tau/200)*log2(1+SINRLSFDRZFMonteCarlo);
    RateRZFMonteCarlosave(:,:,iter) = (1 -tau/200)*log2(1+SINRRZFMonteCarlo);
end
ResultsFig10.SumRateLSFDMRCMonteCarlosave=sum(squeeze(sum(RateLSFDMRCMonteCarlosave,1)),1)/nbrBSs;
ResultsFig10.SumRateMRCMonteCarlosave=sum(squeeze(sum(RateMRCMonteCarlosave,1)),1)/nbrBSs;

ResultsFig10.SumRateLSFDRZFMonteCarlosave=sum(squeeze(sum(RateLSFDRZFMonteCarlosave,1)),1)/nbrBSs;
ResultsFig10.SumRateRZFMonteCarlosave=sum(squeeze(sum(RateRZFMonteCarlosave,1)),1)/nbrBSs;

save Results_Fig10.mat ResultsFig10;

%% ==========Plot Fig. 10=============================
F1=cdfplot(ResultsFig10.SumRateLSFDRZFMonteCarlosave);
set(F1,'Color','b', 'LineStyle','-');
hold on;
F11 =cdfplot(ResultsFig10.SumRateRZFMonteCarlosave);
set(F11,'Color','b', 'LineStyle','--');
F2=cdfplot(ResultsFig10.SumRateLSFDMRCMonteCarlosave);
set(F2,'Color','b', 'LineStyle','-');
F22=cdfplot(ResultsFig10.SumRateMRCMonteCarlosave);
set(F22,'Color','b', 'LineStyle','--');
legend('Two-layer decoding', 'Single-layer decoding');
xlabel('Sum SE per cell [b/s/Hz]');
ylabel('CDF');
xlim([10,36])