%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2012, Joseph T. Lizier
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

% = Example 5 - Multivariate transfer entropy on binary data =

% Multivariate transfer entropy (TE) calculation on binary data using the discrete TE calculator:

% Modified example with more than 2 columns

numObservations = 1000000;

base = 2; %base=2 for binary data
sourceColumns = 3;
destColumns = 4;

% Create or import data
sourceArray = (rand(numObservations,sourceColumns)>0.5)*1;
destArray = (rand(numObservations,destColumns)>0.5)*1;

% We need to construct the joint values of the dest and source before we pass them in,
%  and need to use the matrix conversion routine when calling from Matlab/Octave:
mUtils = javaObject('infodynamics.utils.MatrixUtils');
obs1 = mUtils.computeCombinedValues(octaveToJavaIntMatrix(sourceArray), base);
obs2 = mUtils.computeCombinedValues(octaveToJavaIntMatrix(destArray), base);

% Create a TE calculator and run it:
sourceSymbols = base^sourceColumns;
destSymbols = base^destColumns;
maxSymbols = max(sourceSymbols,destSymbols);
teCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', maxSymbols, 1);
teCalc.initialise();
teCalc.addObservations(obs1, obs2);
result2 = teCalc.computeAverageLocalOfObservations()

