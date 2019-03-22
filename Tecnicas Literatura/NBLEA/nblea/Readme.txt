Readme file for Nested Bilevel Evolutionary Algorithm (N-BLEA)



Quick instruction for execution:
------------------------------------------------------------------------
Code the upper level optimization task in ulExternalProblem.m
Code the lower level optimization task in llExternalProblem.m
xu and xl are the upper and lower level decision vectors respectively
Provide the problem and algorithm parameters in externalProblem.m
Execute it as: externalProblem()
------------------------------------------------------------------------



N-BLEA
------------------------------------------------------------------------
N-BLEA is a nested strategy that uses evolutionary algorithm at both levels to handle bilevel optimization problems. However, the N-BLEA approach is intelligent such that it switches to quadratic programming at lower level whenever possible. It is able to handle lower dimensional SMD-Suite as well as the TP-Suite successfully. However, it being a nested approach, the function evaluations are high, particularly for the lower level. This approach has been used in the following two papers that are in press.

Authors: Ankur Sinha, Pekka Malo, and Kalyanmoy Deb.
Paper Title: Test problem construction for single-objective bilevel optimization.
Journal: Evolutionary Computation Journal.
Year: 2013.

Authors: Ankur Sinha, Pekka Malo, Anton Frantsev, and Kalyanmoy Deb.
Paper Title: Finding optimal strategies in multi-period multi- leader-follower stackelberg games using an evolutionary framework.
Journal: Computers and Operations Research.
Year: 2013.

The algorithm provided in this package is almost similar to what has been discussed in the above two papers, but it differs in the following ways:
1. QP is implemented at the lower level along with a GA.
2. Utilizes different kind of termination criterion at both levels.
3. Lower level searches are performed with a reduced lower level population size in the inermediate upper level generations.

The above modifications help in faster execution of the algorithm. It increases its efficiency and in most of the cases produces better results than what is reported in the above two papers.
------------------------------------------------------------------------



Files in the package
------------------------------------------------------------------------
There are the following Matlab (.m and .p) files in this package:

ulSearch.p: Obscured version of ulSearch.m that performs search at the upper level.
llSearch.p: Obscured version of llSearch.m that performs search at the lower level.
quadApprox.p: Obscured version of quadApprox.m that is a supporting file for quadratic programming.

ulTestProblem.m: Source code for upper level SMD and TP suite.
llTestProblem.m: Source code for lower level SMD and TP suite.
ulExternalProblem.m: Source code for upper level optimization task for a user defined problem.
ulExternalProblem.m: Source code for lower level optimization task for a user defined problem

smd1.m - smd12.m: These files contain the problem and algorithm parameters for SMD-Suite.
tp1.m - tp10.m: These files contain the problem and algorithm parameters for TP-Suite.
externalProblem.m: It contains the problem and algorithm parameters for a user defined problem.
------------------------------------------------------------------------



Executing a user-defined problem
------------------------------------------------------------------------
To execute a user-defined problem, code the upper level optimization task in ulExternalProblem.m and the lower level optimization task in llExternalProblem.m. The functions inside the files contain arguments xu and xl, which represent the upper level decision vector and lower level decision vector respectively.

Provide the problem and algorithm parameters for user-defined problem in externalProblem.m and call the following command to execute the user-defined bilevel optimization task.
externalProblem()

The results of the execution are printed on the screen as well as stored in 'externalProblem.mat'. A sample bilevel optimization problem is already coded as an external problem.

Problem parameters to be defined in externalProblem.m
ulDim: Number of dimensions at upper level.
llDim: Number of dimensions at lower level.
ulDimMax: Vector defining maximum values for upper level variables.
ulDimMin: Vector defining minimum values for upper level variables.
llDimMax: Vector defining maximum values for lower level variables.
llDimMin: Vector defining minimum values for lower level variables.

Algorithm parameters to be defined in externalProblem.m
ulPopSize: Upper level population size.
llPopSize: Lower level population size.
ulMaxGens: Upper level maximum generations.
llMaxGens: Lower level maximum generations.
ulStoppingCriteria: Stopping parameter at upper level. Smaller the value higher the accuracy.
llStoppingCriteria: Stopping parameter at lower level. Smaller the value higher the accuracy.

There are other parameters in the algorithm. However, most of them are either adaptive or not necessary to be adjusted.

Output of the execution
ulEliteFunctionValue: Upper level function value for the elite member.
llEliteFunctionValue: Lower level function value for the elite member.
ulEliteIndiv: Upper level elite member.
llEliteIndiv: Lower level elite member.
ulFunctionEvaluations: Upper level function evaluations required during the exection.
llFunctionEvaluations: Lower level function evaluations required during the exection.
------------------------------------------------------------------------



Executing the SMD or TP suite
------------------------------------------------------------------------
To execute one of the SMD test problems (say SMD1), the following command needs to be called:
smd1()
This executes the SMD1 test problem with problem and algorithm parameters coded in smd1.m. It executes a 5 variable SMD1 problem with 2 upper level variables and 3 lower level variables. The results are printed on the screen as well as stored in 'smd1.mat'

To execute one of the problems in TP-Suite (say SMD1), the following command needs to be called:
tp1()
This executes tp1 with problem and algorithm parameters coded in tp1.m. The results are printed on the screen as well as stored in 'tp1.mat'
------------------------------------------------------------------------



Contact
------------------------------------------------------------------------
In case you have any questions, comments, suggestions, or you want to report any bugs, you can send an email to Ankur Sinha (Ankur.Sinha@aalto.fi)

Ankur Sinha, PhD
Aalto University School of Business
Helsinki Finland
------------------------------------------------------------------------