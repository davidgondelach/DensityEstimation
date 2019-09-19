% Test BC estimation performance
addpath( 'AstroFunctions' );
addpath( 'nrlmsise_matlab' );
i = 1;
    BCestimationPeriodInDays = 1;
    satrecIndex = find([objects(i).satrecs.jdsatepoch]>=objects(i).satrecs(1).jdsatepoch+BCestimationPeriodInDays,1,'first');
    [ estimate ] = estimateBallisticCoefficient( objects(i).satrecs(satrecIndex), objects(i).satrecs(1), 1 );
    BCestimates(i) = estimate(end);