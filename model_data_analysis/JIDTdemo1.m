% Change location of jar to match yours:
javaaddpath('../infodynamics-dist-1.3/infodynamics.jar');
% Generate some random binary data.
% Note that we need the *1 to make this a number not a Boolean,
% otherwise this will not work (as it cannot match the method signature)
sourceArray=(rand(100,1)>0.5)*1;
destArray = [zeros(50,1); sourceArray(1:50)];
sourceArray2=(rand(100,1)>0.5)*1;
% Create a TE calculator and run it:
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 2, 1);
teCalc.initialise();
% Since we have simple arrays of doubles, we can directly pass these in:
teCalc.addObservations(sourceArray, destArray);
fprintf('For copied source, result should be close to 1 bit : ');
result = teCalc.computeAverageLocalOfObservations()
teCalc.initialise();
teCalc.addObservations(sourceArray2, destArray);
fprintf('For random source, result should be close to 0 bits: ');
result2 = teCalc.computeAverageLocalOfObservations()