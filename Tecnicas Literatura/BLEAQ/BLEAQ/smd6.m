function [] = smd6(ulPopSize, llPopSize, ulMaxGens, llMaxGens, ulDim, llDim)

problemName = 'smd6';             % Test problem name

if nargin==0
    ulPopSize=50;                 % Size of UL population
    ulMaxGens=2000;               % Maximum number of generations allowed at UL
    ulDim=5;                      % Number of UL dimensions
    llPopSize=50;                 % Size of LL population
    llMaxGens=2000;               % Maximum number of generations allowed at LL
    llDim=5;                      % Number of LL dimensions
end

r = floor(ulDim/2);
p = ulDim - r;
q = floor((llDim - r)/2 - eps);
s = ceil((llDim - r)/2 + eps);

size_xu1 = p;
size_xu2 = r;
size_xl1 = q+s;
size_xl2 = r;

ulDimMin = -5*ones(1,ulDim);      % Minimum bound accross dimensions
ulDimMax = 10*ones(1,ulDim);      % Maximum bound accross dimensions

llDimMin = -5*ones(1,llDim);      % Minimum bound accross dimensions
llDimMax = 10*ones(1,llDim);      % Maximum bound accross dimensions

ulStoppingCriteria = 1e-3;
llStoppingCriteria = 1e-5;

[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria)

save('smd6');
