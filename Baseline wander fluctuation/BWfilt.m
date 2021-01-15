function EMG_proper = BWfilt(sigToFilt, varargin)
%% Based on Fasano, A., & Villani, V. (2014). Baseline wander removal for bioelectrical signals by quadratic variation reduction. Signal Processing, 99, 48-57.
% This function filter out baseline wander also called baseline fluctation
% from classic EMG signal
% Contact: lucien.robinault@protonmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   sigToFilt --> Signals to filter (must be of a shape N*1
%       setpL --> The step used in the gradient descent function to find the
%                 optimal lambda. The early test showed good result with
%                 lambda around the [600:700] mark.
%                 default value = 25;
% formulaUsed --> Formula used for the gradient descent (more details
%                 below)
%                 1 = s^x
%                 2 = tanh(a*s^b)
%                 3 = -log10(1+s)
%                 default value = 3;
%           x --> parameter from the sparsity formula 1
%                 default value = 2;
%           a --> parameter from the sparsity formula 2
%                 default value = 1;
%           b --> parameter from the sparsity formula 2
%                 default value = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMG_proper --> EMG signals filtered

%% DISCLAIMERS
% CALCUL OPTIMIZATION:
%     This filter uses a matrix inversion of a matrix of the size N*N, thus
%     being extremely costly memory and computationally wise. Because the 
%     matrix is a tridiagonal, symmetric, positive-definite system it is 
%     possible to “linearize” the calculus. To do so the equation (16) have
%     been developed,and the linear system solver “\” has been used to “get rid of”
%     the standardinverse calculation. This is possible due to the 
%     aforementioned specific properties of the system.
%     More details and explanation can be found in the books from Golub on 
%     matrix computation.
%     Nonetheless, clear optimization process not being detailed in
%     Fasano et al, 2014, the optimization process implemented here might not be
%     the one used by the author of the article

% GRADIENT DESCENT:
%     The process to find the optimal lambda for the filter is not detailed 
%     in Fasano et al., 2014. Therefore, the process implemented here is of 
%     my own design.
%     The logic behind it is, we are facing a signal away from the baseline
%     it should be on, which is 0. The filter is trying to get the signal 
%     back on this baseline. As the signal get filtered closer and closer to
%     the baseline, the number of 0 value increase.
%     We therefore can formulate the problem as a sparsity problem. 
%     The gradient descent work on the following logic: 
%     if sparsity of sFilt(n-1) - sparsity of sFilt(n) < 0 then stop
%     we found the optimal lambda.
%     The formula used have been chosen from Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity. IEEE Transactions on Information Theory, 55(10), 4723-4741, and selected for their behavior that seemed to get close to the optimal lambda found experimentally.
%     The default formula (3) seems to be the best. 
%     Customization of the parameters of the formula (1) and (2) have been
%     made possible but is not recommended.
%     Multiple formula had been chosen in the first place due to them being
%     used on an unoptimized algorithm and therefore were showing different
%     computation time. This is not relevant anymore and formula (1) and (2)
%     might be deleted if extensive testing show that formula (3) does the job.

%% INITIALIZATION
% Default value for the parameters
defaultStepL = 25;
defaultFormulaUsed = 3;
expectedFormula = [1, 2 ,3];
defaulta = 1;
defaultb = 2;
defaultx = 2;

% Function handle for the testing if parameter value is numeric
validNumericParameter = @(x) isnumeric(x);

% Create the input parser (IP) object
inParser = inputParser;

% Requiered input
addRequired(inParser,'sigToFilt');

% Main parameters
addParameter(inParser,'stepL',defaultStepL,validNumericParameter);
addParameter(inParser,'formulaUsed',defaultFormulaUsed,@(formulaChoice) (ismember(formulaChoice,expectedFormula)));

% Formula 1 parameter
addParameter(inParser,'x',defaultx,validNumericParameter);

% Formula 2 parameters
addParameter(inParser,'a',defaulta,validNumericParameter);
addParameter(inParser,'b',defaultb,validNumericParameter);

% Parse the argument into the IP object
parse(inParser,sigToFilt,varargin{:});

lengthSignal = length(inParser.Results.sigToFilt);


switch inParser.Results.formulaUsed
    %% FORMULA 1
    case 1
       % First step of the filter
        lambda = 0.1;
        for stepInit = 1:1:3
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Compute the sparsity value of the filtered signal
            sparsity(stepInit) = sqrt(sum(EMG_proper.^inParser.Results.x));
            
            % Increase lambda by stepL
            lambda = lambda + inParser.Results.stepL;
        end
        
        % Filtering
        while (diff([(diff([sparsity(1) sparsity(2)]).*sparsity(1)) (diff([sparsity(2) sparsity(3)]).*sparsity(2))])>0)
            % Increase lambda by stepL
            lambda = lambda+inParser.Results.stepL;
            
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Refresh the sparsity values to n-2, n-1, n
            sparsity(1) = sparsity(2);
            sparsity(2) = sparsity(3);
            sparsity(3) = sqrt(sum(EMG_proper.^inParser.Results.x));
        end
    
    %% FORMULA 2
    case 2
        % First step of the filter
        lambda = 0.1;
        for stepInit = 1:1:3
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Compute the sparsity value of the filtered signal
            sparsity(stepInit) = -sum(tanh((inParser.Results.a*EMG_proper).^inParser.Results.b));
            
            % Increase lambda by stepL
            lambda = lambda + inParser.Results.stepL;
        end
        % Filtering
        while (diff([(diff([sparsity(1) sparsity(2)]).*sparsity(1)) (diff([sparsity(2) sparsity(3)]).*sparsity(2))])>0)
            % Increase lambda by stepL
            lambda = lambda+inParser.Results.stepL;
            
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Refresh the sparsity values to n-2, n-1, n
            sparsity(1) = sparsity(2);
            sparsity(2) = sparsity(3);
            sparsity(3) = -sum(tanh((inParser.Results.a*EMG_proper).^inParser.Results.b));
        end

    %% FORMULA 3
    case 3
        % First step of the filter
        lambda = 0.1;
        for stepInit = 1:1:3
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Compute the sparsity value of the filtered signal
            sparsity(stepInit) = sum(-log(1+EMG_proper));
            
            % Increase lambda by stepL
            lambda = lambda + inParser.Results.stepL;
        end
        
        % Filtering
        while (diff([(diff([sparsity(1) sparsity(2)]).*sparsity(1)) (diff([sparsity(2) sparsity(3)]).*sparsity(2))])>0)
            % Increase lambda by stepL
            lambda = lambda+inParser.Results.stepL;
            
            % Creating the compressed matrix I+L*D'D from (16)
            systemD = spdiags([-lambda*ones(1,lengthSignal); [lambda+1, 2*lambda*ones(1,lengthSignal-2)+1, lambda+1]; -lambda*ones(1,lengthSignal)]',-1:1,lengthSignal,lengthSignal);
            % Eq (16) developped and optimized from Fasano et al., 2014
            % This is the actual filtering process
            EMG_proper = inParser.Results.sigToFilt-systemD\inParser.Results.sigToFilt;
            
            % Refresh the sparsity values to n-2, n-1, n
            sparsity(1) = sparsity(2);
            sparsity(2) = sparsity(3);
            sparsity(3) = sum(-log(1+EMG_proper));
        end
end
end
