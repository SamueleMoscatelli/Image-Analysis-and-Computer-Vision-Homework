function [ lines ] = lineDetector( imInput , minLength)
% Identify the straight line segments with Hough Transform.

im = imInput;

[H,T,R]= hough(im); % Hough Transform

%P = houghpeaks(H,25,'threshold',ceil(0.1*max(H(:))),'NHoodSize',[73,11]);
P = houghpeaks(H,20,'threshold',ceil(0.01*max(H(:))),'NHoodSize',[7,7]);

lines = houghlines(im,T,R,P,'FillGap',11,'MinLength',minLength);
%lines = houghlines(im,T,R,P);

end

