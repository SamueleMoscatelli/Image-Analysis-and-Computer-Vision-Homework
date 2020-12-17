function [lines] = computeLines(edgs, minLength, cs)
switch(cs)
    case 1
        % Compute the Hough transform of the binary image returned by edge.
        [H,theta,rho] = hough(edgs);
        % Find the peaks in the Hough transform matrix, H, using the houghpeaks function.
        P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
        % Find lines in the image using the houghlines function.
        lines = houghlines(edgs,theta,rho,P,'FillGap',70,'MinLength',minLength);
        %lines = houghlines(edgs,theta,rho,P,'FillGap',70,'MinLength',10);
    case 2
        minLineLength_vert = 160;
        fillGap = 20;
        numPeaks = 40;
        NHoodSize = [101 51];
        vertical_angles = -90:0.5:89.8;

        % find lines vertical
        [H,theta,rho] = hough(edgs,'RhoResolution', 1, 'Theta', vertical_angles);
        % find peaks in hough transform
        P = houghpeaks(H,numPeaks,'threshold',ceil(0.05*max(H(:))), 'NHoodSize', NHoodSize);

        % find lines using houghlines
        lines_ = houghlines(edgs,theta,rho,P,'FillGap',fillGap,'MinLength',minLineLength_vert);

        lines = [lines_(1,[41,40,54])]; %41,40,54,11,

end