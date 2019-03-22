function [StoppingCriteria stoppingCondition] = terminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    %If StoppingCriteria is returned 1 then the code terminates otherwise not
    StoppingCriteria = 0;
    stoppingCondition = [];
    [improvementStoppingCriteria improvementStoppingCondition] = improvementTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters);
    if improvementStoppingCriteria==1
        StoppingCriteria = 1;
        stoppingCondition = improvementStoppingCondition;
    end
    [varianceStoppingCriteria varianceStoppingCondition] = varianceTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters);
    if varianceStoppingCriteria==1
        StoppingCriteria = 1;
        stoppingCondition = varianceStoppingCondition;
    end
    
    
function [StoppingCriteria stoppingCondition] = varianceTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    ulEpsilonStopping = stoppingParameters.ulEpsilonStopping;
    llEpsilonStopping = stoppingParameters.llEpsilonStopping;
    alphaStoppingInitial = stoppingParameters.alphaStoppingInitial;
    eliteConstrViolation = stoppingParameters.eliteConstrViolation;
    
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    if sum(tag.ulPop==1) <= 1
        alphaStopping = Inf;
    else
        alphaStopping = sum(var([ulPop(tag.ulPop==1,:) llPop(tag.ulPop==1,:)]))/(ulDim+llDim);
        alphaStopping = alphaStopping/alphaStoppingInitial;
    end
    if alphaStopping<ulEpsilonStopping & eliteConstrViolation<=0
        StoppingCriteria = 1;
        disp('Variance based termination criterion at upper level met.')
        stoppingCondition = 'Variance based'
    end
    
    
function [StoppingCriteria stoppingCondition] = improvementTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    ulEpsilonStopping = stoppingParameters.ulEpsilonStopping;
    llEpsilonStopping = stoppingParameters.llEpsilonStopping;
    eliteFunctionValueAtGen = stoppingParameters.eliteFunctionValueAtGen;
    eliteConstrViolation = stoppingParameters.eliteConstrViolation;
    
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    %Stops if normalized improvement is less than ulEpsilonStopping in ulImprovementGenDiff generations
    ulEpsilonStoppingImprovement = ulEpsilonStopping;
    ulImprovementGenDiff = 200;

    if (gen>ulImprovementGenDiff)
        if abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff)) == 0
            betaStopping = 0;
        elseif abs(eliteFunctionValueAtGen(gen)+eliteFunctionValueAtGen(1)) == 0
            betaStopping = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff));
        else
            betaStopping = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff))/abs(eliteFunctionValueAtGen(gen)+eliteFunctionValueAtGen(1));
        end
    else
        betaStopping = Inf;
    end
    if betaStopping<ulEpsilonStoppingImprovement & eliteConstrViolation<=0
        StoppingCriteria = 1;
        display('Lower level search done. Improvement based termination condition met.');
        stoppingCondition = 'Improvement based'
    end