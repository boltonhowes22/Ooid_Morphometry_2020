function statistical_flattening = flattening(long_axis, short_axis)
% This function calcluates the prolateness of an ellipse using the
% equation: (long_axis - short_axis) / long_axis
% INPUTS
% long_axis: long axis of ellipse
% short_axis: short axis of ellipse
%
% OUT
% statistical_flattening: the flattening value of the ellipsoe 
% Example: 
% long_axis = 100
% short_axis = 50% radius_3 = 100
% scaling_factor = 30
% statistical_flattening = flattening(long_axis, short_axis)
%
% Bolton Howes
% February 2019
% 
% ======================== Begin Function ==============================

statistical_flattening = (long_axis - short_axis) ./ long_axis;
