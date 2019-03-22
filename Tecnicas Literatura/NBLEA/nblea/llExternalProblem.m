function [functionValue equalityConstrVals inequalityConstrVals] = llExternalProblem(xu, xl)

    %Lower level TP1 implemented
    %xu is an upper level decision vector
    %xl is a lower level decision vector
    functionValue = (xu(1) - xl(1))^2 + (xu(2) - xl(2))^2;

    functionValue = -functionValue;

    equalityConstrVals = [];
    inequalityConstrVals = [];
