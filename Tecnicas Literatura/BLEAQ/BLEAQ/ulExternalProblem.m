function [functionValue equalityConstrVals inequalityConstrVals] = ulExternalProblem(xu, xl)

    %Note: All equality and inequality constraints have to be initialized
    %in the form of an array. If there are no equality or inequality
    %constraints then the variables have to be marked as empty. The upper
    %level vector is represented by xu and the lower level vector is
    %represented by xl,

    %Upper level TP1 implemented
    functionValue = (xu(1)-30)^2 + (xu(2)-20)^2 - 20*xl(1) + 20*xl(2);
        
    %Putting a negative sign because the function has to be minimized but
    %the BLEAQ code is for maximization.  Following line will not be
    %required if you are maximizing the upper level
    functionValue = -functionValue;
    
    %If there are no equality constraints make the following variable empty
    equalityConstrVals = [];
    
    %Assign the inequality constraints through an array
    inequalityConstrVals(1) = 30 - xu(1) - 2*xu(2);
    inequalityConstrVals(2) = xu(1) + xu(2) - 25;