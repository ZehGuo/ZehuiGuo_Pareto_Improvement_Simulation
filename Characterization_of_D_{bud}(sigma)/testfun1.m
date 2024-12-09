function [outputy,outputb] = testfun1(inputv1,inputv2,inputx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
outputy=-(inputv1(1)*inputx+inputv1(3))/inputv1(2);
outputb=-inputv2(1)*inputx-inputv2(2)*outputy;
end