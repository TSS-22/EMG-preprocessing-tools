function EMG_proper = EMGccaFilt(sigToFilter, fqAcq)
%% Based on Al Harrach, M., Boudaoud, S., Hassan, M., Ayachi, F. S., Gamet, D., Grosset, J. F., & Marin, F. (2017). Denoising of HD-sEMG signals using canonical correlation analysis. Medical & biological engineering & computing, 55(3), 375-388.
% This function remove movement artifact, white noise (WGN) and power line inference (PLI)
% Contact: lucien.robinault@protonmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%/!\IMPORTANT/!\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The function need the first 0.5 seconds of signals to contain no muscle
% activity.
% ECG contamination prevent successful use most of the time, it is vividly
% recommended to filter them out before. A function serving this purpose is
% available in this repository
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% line interference (PLI) from HD EMG signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigToFilter --> Signals to filter, with COLUMNS = Samples & ROWS = Sources
% fqAcq     --> Acquisition frequency of the EMG signals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMG_proper --> EMG signals filtered, with COLUMNS = Samples & ROWS = Sources


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART I: CANONICAL CORRELATION ANALSYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rowSig, colSig] = size(sigToFilter);

% Create a lagged signal to improve quality of the CCA result
laggedSigToFilter = [sigToFilter(:, end) sigToFilter(:, 1:(end-1))];

% CCA operation
[A, B, r, U, V] = canoncorr(sigToFilter', laggedSigToFilter');
[lenU, compU] = size(U);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART II: CCA COMPONENT THRESHOLDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intensityRatio = zeros(1,compU);

%% Calculus of the Inensity Ratio (IR)
for i = 1:1:compU
    % Equation (6) from Al Harrach et al., 2016
    intensityRatio(i) =  (sum(abs(U(:,i)))/length(U))*...
        ((fqAcq*0.5)/sum(abs(U(1:round(fqAcq*0.5),i))));
    
end

% Check which ratio threshold is the best suited, from 1 to 4
for threshold = 10:1:40
    % Check which CCA components are below or equal to the IR
    compBelowRatio{threshold} = intensityRatio<=threshold/10;
    
    % Put the CCA components who are <= threshold to zero
    tempU = U;
    tempU(:,compBelowRatio{threshold}) = 0;
    % reconstruction of the EMG signals (temporary) to test the correlation
    % of these new signals to the non filtered ones
    cleanSignalTest = tempU / A;
    
    % Test the correlation between these filtered signals and the
    % unfiltered ones.
    corrSignals{threshold} = corrcoef(cleanSignalTest,sigToFilter');
    % If correlation is < 0.8 the CCA filtering stop and we take the last
    % CCA filtering setup that had a correlation of >= 0.8.
    % If we go up to a threshold of 4 we keep the threshold 4 (therefore
    % the last.
    if ~(corrSignals{threshold}(1,2)>=0.8)
        finalThreshold = threshold-1;
        
        break
    else
        finalThreshold = threshold;
    end
end

% Filtering the HD EMG signals by removing the selected components
tempU = U;
tempU(:,compBelowRatio{finalThreshold}) = 0;
cleanSignal = (tempU / A)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART III: SELECTIVE CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if filtered signals got better PNR than the non filtered ones, in
% order to prevent noise contamination during reconstruction of the
% signals.
% This is done channel/electrode by channel/electrode
for i = 1:1:rowSig
    % Calculus of the non filtered channel PNR
    dirtySNR = 20*log10(((sum(abs(sigToFilter(i,:))))/colSig)*...
        ((0.5*fqAcq)/(sum(abs(sigToFilter(i,1:round(fqAcq*0.5)))))));
    
    % Calculus of the filtered channel PNR
    cleanSNR = 20*log10(((sum(abs(cleanSignal(i,:))))/colSig)*...
        ((0.5*fqAcq)/(sum(abs(cleanSignal(i,1:round(fqAcq*0.5)))))));
   
    % Compare PNR between those two and add the highest one to the
    % definitive filtered HD EMG matrix (EMG_proper)
    if dirtySNR>=cleanSNR
        EMG_proper(i,:) = sigToFilter(i,:);
        
    else
        EMG_proper(i,:) = cleanSignal(i,:);
    end

end
