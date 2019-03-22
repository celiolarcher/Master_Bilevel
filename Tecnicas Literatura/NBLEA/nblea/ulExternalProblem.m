function [functionValue equalityConstrVals inequalityConstrVals] = ulExternalProblem(xu, xl)

    %Upper level TP1 implemented
    functionValue = (xu(1)-30)^2 + (xu(2)-20)^2 - 20*xl(1) + 20*xl(2);
        
    functionValue = -functionValue;
    
    equalityConstrVals = [];
    inequalityConstrVals(1) = 30 - xu(1) - 2*xu(2);
    inequalityConstrVals(2) = xu(1) + xu(2) - 25;
