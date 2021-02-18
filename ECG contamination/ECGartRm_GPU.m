function [EMG_proper, ECG_proper] = ECGartRm_GPU(sigToFilt, fqAcq)
%% Based on Mak, J. N., Hu, Y., & Luk, K. D. (2010). An automated ECG-artifact removal method for trunk muscle surface EMG recordings. Medical engineering & physics, 32(8), 840-848.
% This function filter out ECG componant from High Density EMG signal
% Contact: lucien.robinault@protonmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigToFilt --> Signals to filter, with COLUMNS = Samples & ROWS = Sources
% fqAcq     --> Acquisition frequency of the EMG signals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMG_proper --> EMG signals filtered, with COLUMNS = Samples & ROWS = Sources
% ECG_Proper --> ECG signals taken out, with COLUMNS = Samples & ROWS = Sources

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: Signal pre-processing and fastICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fastICA
[ICAsig, sepMat, ~] = fastica(gpuArray(sigToFilt));
ICAsig = gather(ICAsig);
sepMat = gather(sepMat);

[rowIcaSig, colIcaSig] = size(ICAsig);
EMGvsECG = zeros(rowIcaSig, 1);

% Pre processing based on
% Myers, L. ., Lowery, M., O’Malley, M., Vaughan, C. ., Heneghan, C., St Clair Gibson, A., … Sreenivasan, R. (2003). Rectification and non-linear pre-processing of EMG signals for cortico-muscular analysis. Journal of Neuroscience Methods, 124(2), 157–165. doi:10.1016/s0165-0270(03)00004-9 
% as seen in Mak et al., 2010
for i=1:1:rowIcaSig
    % Hilbert transform 
    sigToFiltH(i,:) = hilbert(ICAsig(i,:)./max(ICAsig(i,:)));
    
    % get the a value of the complex number
    sigToFiltH(i,:) = abs(sigToFiltH(i,:));
    
    % Median filter order 50
    sigToFiltH(i,:) = medfilt1(sigToFiltH(i,:),50);
    
    % THIS PART HAS BEEN ADDED BY EMPIRICAL EVIDENCE
    % IN ORDER FOR THE FILTER TO WORK IN THE FIELD
    % Moving average filter to smooth potential artifact on the ECG signal
    % (multiple crests) that can appear in field conditions
    windowWidth = round(fqAcq*0.05); 
    kernel = ones(windowWidth,1) / windowWidth;
    sigToFiltH(i,:) = filter(kernel, 1, sigToFiltH(i,:));
end

S1 = sigToFiltH;

for i = 1:1:rowIcaSig
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 2: Peak detection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Searching  for  peaks  in  a  source  component  signals(n)was
    % accomplished by a simple peak detection algorithm as follows:
    
    % 1. Scan the signals(n) which may be expected to contain a series
    %   of peaks and determine the maximum value Smax.
    Smax = max(sigToFiltH(i,:));
    
    % 2.Define a threshold as a fraction of the maximum,Th = 0.6 Smax
    threshold = 0.6*Smax;
    
    % 3. Convert the signal into binary format:
    %       if S1(n)>= Th, S(n)=1
    %       if S1(n)< Th, S(n)=0
    S1(i,sigToFiltH(i,:)>=threshold) = 1;
    S1(i,sigToFiltH(i,:)<threshold) = 0;
    
    % 4. Calculate the rate of change of signal S1(n),
    %   i.e. the first derivativeof S(n), which is approximated as:
    %   s2(n)=s1(n)-s1(n-1),n=2,3, ..., N
    %   where N is the number of samples
    S2(i,1) = 0;
    for j = 2:1:colIcaSig
        S2(i,j) = S1(i,j) - S1(i,j-1);
    end
    
    % 5. Select those samples for which the corresponding S2(n) value is equal
    %   to one, that is, having a positive rate of change:
    %   P={n|s2(n)=1}
    %   The setPdefine as above contains the indices of the peaks ins(n).
    P{i} = find(S2(i,:) == 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 3: Identification of ECG source components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Number of peaks:
    %   (200 BPM/60 s)·d>=|P|>=(40 BPM/60 s)*d
    %   where |P| indicates the number of elements in the set P, that is
    %   the number of peaks detected,
    %   and d represents the length of source component signal (in second).
    limPeaksLow = (40/60)*(colIcaSig/fqAcq);
    limPeaksHigh = (200/60)*(colIcaSig/fqAcq);
    
    % 2. RR interval:
    %   1.5 s>=P(n+1)-P(n)>=0.3s, n=1,2,...,N
    %   where P(n) represents the time information of thenth detected peak.
    %   Nis the number of peaks detected. 1.5 s is the averaged RR interval
    %   value with a heart rate of 40 BPM and 0.3 s is
    %   the averaged RR interval value with a heart rate of 200 BPM.
    limRRlow = 0.3*fqAcq;
    limRRhigh = 1.5*fqAcq;
    
    % 3. Variance of RR intervals:
    %   [P(n+2)-P(n+1)]-[P(n+1)-P(n)]<=R*(1.5s),n=1,2,..., N
    %   where 1.5 s is the upper limit of the RR interval value.
    %   A scaling factor R of 0.5 was adopted in this study.
    R = 0.5;
    limVarRR = 0.5*1.5*fqAcq;
    
    % Is the ICA signal a ECG component?
    check1 = false;
    check2 = false;
    check3 = false;
    
    % Check Number of peaks
    if (length(P{i})<=limPeaksHigh)&&(length(P{i})>=limPeaksLow)
        check1 = true;
    end
    
    % Check RR interval
    intervalCheck = zeros(length(P{i})-1,1);
    for j = 1:1:length(P{i})-1
        test2 = P{i}(j+1)- P{i}(j);
        if (test2<=limRRhigh)&&(test2>=limRRlow)
            intervalCheck(j) = 1;
        end
    end
    
    % THE 10% ERROR ALLOWED HAS BEEN ADDED BY EMPIRICAL EVIDENCE
    % IN ORDER FOR THE FILTER TO WORK IN THE FIELD
    % It might be optional though
    if sum(intervalCheck)>=(floor((length(P{i})-1)*0.90))
        check2 = true;
    end

    % Check Variance of RR intervals
    varianceCheck = zeros(length(P{i})-2,1);
    for j = 1:1:length(P{i})-2
        test3 = (P{i}(j+2)-P{i}(j+1)) - (P{i}(j+1)-P{i}(1));
        if test3 <= (R*1.5*fqAcq)
            varianceCheck(j) = 1;
        end
    end
    
    if sum(varianceCheck)==(length(P{i})-2)
        check3 = true;
    end
    
    % Check all the conditions
    % EMGvsECG = 1 --> ECG component
    % EMGvsECG = 0 --> EMG component
    if check1 && check2 && check3
        EMGvsECG(i) = 1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4: Reconstruct EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EMG_proper = sepMat(:,~EMGvsECG)*ICAsig(~EMGvsECG,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3: Reconstruct ECG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ECG_proper = sepMat(:,logical(EMGvsECG))*ICAsig(logical(EMGvsECG),:);

end
