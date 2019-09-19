function [ goodData, goodIndices ] = filterDataUsingMovingMedian( data, span, subsetSize, varargin )
%FILTERDATAUSINGMOVINGMEDIAN - find outliers in the data by sliding a
%window through the data and comparing the distances of the datapoints to
%the median in the window. Use robust statistics.
%
% The median is computed using a different window length than the actual
% data rejection. This enables larger or smaller subsets to be used to
% compute the dispersion of the sample to do outlier rejection. Median
% Absolute Deviation is used to set the outlier rejection threshold because
% this is a more robust metric of sample dispersion than, say, the standard
% deviation.
%
% Syntax: [goodData, goodIndices ] = filterDataUsingMovingMedian(data, span, subsetSize, noMAD)
%
% Inputs:
%    data - array with the data where the outleirs are to be found.
%    span - span of the window to use to detrend the data by computing the
%       moving median nd subtracting it from all the datapoints. Units are
%       the number of datapoints.
%    subsetSize - how many datapoints to use to compute the median absolute
%       deviation and reject outliers in the subset. If set to length(data)
%       MAD of the entire sample will be used.
%
% Outputs:
%    goodData - datapoints that are not outliers.
%    goodIndices - indices of the input data vector that correspond to
%       non-outlying datapoints.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% globals required: none
%
% Author: Aleksander Lidtke
% University of Southampton, UK
% email address: aleksander.lidtke@gmail.com 
% Website: http://www.aleksanderlidtke.com
% November 2015; Last revision: 25-Feb-2016


medianValues = movingMedian(data,span)';

goodIndices = find( data./medianValues < 1.5 & data./medianValues > 0.5 );

goodData = data(goodIndices);

end