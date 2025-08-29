function [gray_level] = photoelectrons2graylevel(number_of_electrons,fw,bit_level)
gray_level=number_of_electrons*(2^bit_level-1)/fw;



