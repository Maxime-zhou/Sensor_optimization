function mutationChildren = mutationpower(parents,options,GenomeLength, ...
    ~,~,~,thisPopulation)
%MUTATIONPOWER Mutation operator for mixed-integer GA.
% 
%   MUTATIONCHILDREN = MUTATIONPOWER(PARENTS, OPTIONS, ...
%   GENOMELENGTH, FITNESSFCN, STATE, THISSCORE, THISPOPULATION) creates the
%   mutation children from the PARENTS in THISPOPULATION for mixed-integer
%   GA problems.
%
%   This function is set as the mutation function in the options structure
%   that is passed to galincon.
%
%   Power mutation mutates a parent, x, via the following. For each
%   component of the parent, the i-th component of the child is given by:
%
%   mutationChild(i) = x(i) - s(x(i) - lb(i)) if t < r
%                    = x(i) + s(ub(i) - x(i)) if t >= r

%
%   where t is the scaled distance of x(i) from the i-th component of the
%   lower bound, lb(i). s is a random variable drawn from a power
%   distribution and r is a random number drawn from a uniform
%   distribution.
% 
%   Note, that this function can handle lb(i) = ub(i). New children are
%   generated with the i-th component set to lb(i) (= ub(i)).
%
%   For more information on this mutation function see section 2.2 of the
%   following reference:
%
%   A real coded genetic algorithm for solving integer and mixed integer
%   optimization problems, Kusum Deep, Krishna Pratap Singh, M.L. Kansal,
%   C.Mohan, Applied Mathematics and Computation, 212 (2009), 505-518

%   Copyright 2011-2021 The MathWorks, Inc.

% Number of mutation children
numMutation = length(parents);

% Variable bounds and range. The mutation function requires finite bounds.
lb = options.LinearConstr.lb;
ub = options.LinearConstr.ub;
[lb, ub] = globaloptim.internal.pCreateCandidatePointBounds(lb(:)', ub(:)');
    
% Initialize mutation chldren
mutationChildren = zeros(numMutation, GenomeLength);

for i = 1:numMutation
    
    % Generate a random number from a power distribution to control the level
    % of mutation on the continuous variables.
    pCnt = 10;
    sCnt = (rand)^pCnt;
    
    % Generate a random number from a power distribution to control the level
    % of mutation on the integer variables.
    % pInt = 4;
	
	pInt = 2;
    sInt = (rand)^pInt;
    
    % Scaled distance of each variable to its lower bound
    parent = thisPopulation(parents(i),:);
    scaleDist = (parent - lb)./(ub - lb);

    % Generate mutated continuous variables    
    for j = options.LinearConstr.ContinuousVars
        if scaleDist(j) < rand
            mutationChildren(i, j) = parent(j) - sCnt*(parent(j) - lb(j));
        else
            mutationChildren(i, j) = parent(j) + sCnt*(ub(j) - parent(j));
        end
    end
        
    % Generate mutated integer variables
    for j = options.LinearConstr.IntegerVars
        if scaleDist(j) < rand
            mutationChildren(i, j) = parent(j) - sInt*(parent(j) - lb(j));
        else
            mutationChildren(i, j) = parent(j) + sInt*(ub(j) - parent(j));
        end
    end
        
    % Enforce integer restrictions, ensuring that rounded variables lie
    % within the bounds.
    mutationChildren(i, :) = gaminlpsatisfyintconstr(mutationChildren(i, :), ...
        lb, ub, options.LinearConstr.IntegerVars);    

end

