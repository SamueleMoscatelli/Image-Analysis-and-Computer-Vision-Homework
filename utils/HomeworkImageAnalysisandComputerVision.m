close all
clear
clc
%% Reading the image
y = im2double(imread('CastelloDiMiramare.jpg'));
sizev = size(y);
x = y;
z = y;
f = y;

x([1500:end],[1:end]) = imgaussfilt(x([1500:end],[1:end]),5);
z = imgaussfilt(z, 7);
w = histeq(y);
w(1:end,3000:end) = imgaussfilt(w(1:end,3000:end),8);
%% Converting image to grayscale
I = rgb2gray(y);
J = rgb2gray(x);
K = rgb2gray(z);
W = rgb2gray(w);
%% Edge detection
edgs1 = edge(I, 'canny', [0.09 0.1], 1.4142);
edgs2 = edge(J, 'canny', [0.09 0.1], 1.4142);
edgs3 = edge(K, 'canny', [0.09 0.1], 1.4142);
edgs4 = edge(W, 'canny', [0.09 0.5], 1.4142);
%edgs = edge(I, 'canny', [0.09 0.1], 1.4142);
%figure(2), imshow(edgs3)
%% Finding straight lines
lines1 = computeLines(edgs1, 350);
lines2 = computeLines(edgs2, 350);
lines3 = computeLines(edgs3, 300);
lines4 = computeLines(edgs4, 220);
lines = [lines1, lines2, lines3, lines4];

% Create a plot that displays the original image with the lines superimposed on it.
figure(4), imshow(I), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
%% Find corners
corners = detectFASTFeatures(I);
imshow(I); hold all;
plot(corners.selectStrongest(50));

%% 2D rectification - Affine rectification
% y = imread('CastelloDiMiramare.jpg');
% % select a first pair of segments (images of 2 parallel lines)
% segment1 = [lines(1).point1(1) lines(1).point1(2);lines(1).point2(1) lines(1).point2(2)];
% segment2 = [lines(3).point1(1) lines(3).point1(2);lines(3).point2(1) lines(3).point2(2)];
% 
% % select a second pair of segments (images of 2 parallel lines)
% segment3 = [lines(2).point1(1) lines(2).point1(2);lines(2).point2(1) lines(2).point2(2)];
% segment4 = [lines(4).point1(1) lines(4).point1(2);lines(4).point2(1) lines(4).point2(2)];
img = im2double(imread('CastelloDiMiramare.jpg'));
imgR = imresize(img, 0.5);
figure(20); imshow(img);

% select a first pair of segments (images of 2 parallel lines)
% segment1 = [lines(3).point1(1), lines(3).point1(2); lines(3).point2(1), lines(3).point2(2)];
% segment2 = [lines(7).point1(1), lines(7).point1(2); lines(7).point2(1), lines(7).point2(2)];
segment1 = drawline('Color','red');
segment2 = drawline('Color','red');

% select a second pair of segments (images of 2 parallel lines)
segment3 = [lines(1).point1(1), lines(1).point1(2); lines(1).point2(1), lines(1).point2(2)];
segment4 = [lines(4).point1(1), lines(4).point1(2); lines(4).point2(1), lines(4).point2(2)];
% segment3 = drawline('Color','red');
% segment4 = drawline('Color','red');


%compute the image of the line at the infinity
l1 = segToLine(segment1.Position);
l2 = segToLine(segment2.Position);

m1 = segToLine(segment3);
m2 = segToLine(segment4);

% compute the vanishing points
L = cross(l1,l2);
L = L./L(3);
M = cross(m1,m2);
M = M./M(3);

% compute the image of the line at infinity
imLinfty = cross(L,M);
imLinfty = imLinfty./(imLinfty(3));

% dispaly the selection
figure;
hold all;
% plot vanishing points
plot(L(1),L(2),'r.','MarkerSize',100);
plot(M(1),M(2),'b.','MarkerSize',100);
imshow(img);
% plot vanishing line
line([L(1),M(1)],[L(2),M(2)],'Color','Green','Linewidth',3);
% plot selected segments
% line([segment1(1,1),segment1(2,1)],[segment1(1,2),segment1(2,2)],'Color','red','Linewidth',3);
% line([segment2(1,1),segment2(2,1)],[segment2(1,2),segment2(2,2)],'Color','red','Linewidth',3);
% line([segment3(1,1),segment3(2,1)],[segment3(1,2),segment3(2,2)],'Color','blue','Linewidth',3);
% line([segment4(1,1),segment4(2,1)],[segment4(1,2),segment4(2,2)],'Color','blue','Linewidth',3);
hold off;
legend('Vanishing point 1', 'Vanishing point 2','Image of l_\infty');
set(gca,'FontSize',20)

% build the rectification matrix
H = [eye(2),zeros(2,1); imLinfty(:)'];
% we can check that H^-T* imLinfty is the line at infinity in its canonical
% form:

% compute the rectifying matrix
fprintf('The vanishing line is mapped to:\n');
disp(inv(H)'*imLinfty);

%rectify and show the results
tform = projective2d(H');
J = imwarp(imgR,tform);

figure;
imshow(J);
imwrite(J,'AffineRectified.jpg');

hold all;

%% 2D reconstruction - Affine to metric

imgAffRect = im2double(imread('AffineRectified.jpg'));

figure(11);
imshow(imgAffRect);
numConstraints = 2;
sizen = [1488,1984];
imageDim = [sizev(1),sizev(2)];


%normalize coordinates of the points and adapt them to the dimension of the
%new image
line1 = [lines(15).point1;lines(15).point2];
line1 = line1 ./ imageDim;
line1 = line1.* [sizen(1),sizen(2)];
p = [line1(1,1);line1(1,2);1];
p = H * p;
p = p./p(3);
p2 = [line1(2,1);line1(2,2);1];
p2 = H * p2;
p2 = p2./p2(3);
segment5 = [p(1), p(2);p2(1), p2(2)];

hold all;

line(segment5(:,1),segment5(:,2),'Color','red','Linewidth',3);


%normalize coordinates of the points and adapt them to the dimension of the
%new image
line1 = [lines(16).point1;lines(16).point2];
line1 = line1 ./ imageDim;
line1 = line1.* [sizen(1),sizen(2)];
p = [line1(1,1);line1(1,2);1];
p = H * p;
p = p./p(3);
p2 = [line1(2,1);line1(2,2);1];
p2 = H * p2;
p2 = p2./p2(3);
segment6 = [p(1), p(2);p2(1), p2(2)];

line(segment6(:,1),segment6(:,2),'Color','red','Linewidth',3);

%transform segments to lines
l5 = segToLine(segment5);
m6 = segToLine(segment6);
l = [l1, m1];
m = [m6, l5];
doAgain = 1;
count = 1;
A = zeros(numConstraints,3);

while (count <=numConstraints)
    % each pair of orthogonal lines gives rise to a constraint on s
    % [l(1)*m(1),l(1)*m(2)+l(2)*m(1), l(2)*m(2)]*s = 0
    % store the constraints in a matrix A
     A(count,:) = [l(count,1)*m(count,1),l(count,1)*m(count,2)+l(count,2)*m(count,1), l(count,2)*m(count,2)];

    count = count+1;
end

% solve the system
%S = [x(1) x(2); x(2) 1];
[~,~,v] = svd(A)
s = v(:,end); %[s11,s12,s22];
S = [s(1),s(2); s(2),s(3)]

% compute the rectifying homography
imDCCP = [S,zeros(2,1); zeros(1,3)]; % the image of the circular points
[U,D,V] = svd(S);
A = U*sqrt(D)*V';
H = eye(3);
H(1,1) = A(1,1);
H(1,2) = A(1,2);
H(2,1) = A(2,1);
H(2,2) = A(2,2);

Hrect = inv(H);
Cinfty = [eye(2),zeros(2,1);zeros(1,3)];

tform = projective2d(Hrect');
J = imwarp(imgAffRect,tform);

figure(12);
imshow(J);

%% Camera calibration
y = imread('CastelloDiMiramare.jpg');
figure(1), imshow(y),title('original image');
hold all;
% extraction of vertical vanishing point 
% select a first pair of segments (images of 2 parallel lines)
segment7 = [lines(14).point1(1), lines(14).point1(2); lines(14).point2(1), lines(14).point2(2)];
segment8 = [lines(9).point1(1), lines(9).point1(2); lines(9).point2(1), lines(9).point2(2)];

% segment7 = drawline('Color','red');
% segment8 = drawline('Color','red');
n1 = segToLine(segment7);
n2 = segToLine(segment8);

% compute the vanishing point
N = cross(n1,n2);
N = N./N(3);

%line([lines(14).point1(1),N(1)],[lines(14).point1(2),N(2)],'Color','Green','Linewidth',3);

% select the two necessary pair of parallel lines
segment9 = [lines(15).point1(1), lines(15).point1(2); lines(15).point2(1), lines(15).point2(2)];
%segment9 = drawline('Color','red');
o3 = segToLine(segment9);
%o4 = segToLine(segment4.Position);
%compute vanishing point
O = cross(o3,imLinfty);
O = O./O(3);
%line([lines(15).point1(1),O(1)],[lines(15).point1(2),O(2)],'Color','Green','Linewidth',3);


%second pair
segment10 = [lines(16).point1(1), lines(16).point1(2); lines(16).point2(1), lines(16).point2(2)];
%segment10 = drawline('Color','red');
q5 = segToLine(segment10);
%q6 = segToLine(segment5.Position);
%compute vanishing point
Q = cross(q5,imLinfty);
Q = Q./Q(3);
%line([lines(16).point1(1),Q(1)],[lines(16).point1(2),Q(2)],'Color','Green','Linewidth',3);

%IAC = get_IAC(imLinfty,

% dispaly the selection
% plot vanishing points (comprising the 
% ones taken with rectification)
plot(N(1),N(2),'r.','MarkerSize',100);
plot(L(1),L(2),'g.','MarkerSize',100);
plot(M(1),M(2),'b.','MarkerSize',100);
plot(O(1),O(2),'b.','MarkerSize',100);
plot(Q(1),Q(2),'b.','MarkerSize',100);

% constraints on the IAC
syms u1 u2 u3 u4;
IAC = [u1, 0, u2;0, 1, u3;u2, u3, u4];
% constraint 1: L' IAC O = 0
C1 = L' * IAC * O;
% constraint 2: M' IAC Q = 0
C2 = M' * IAC * Q;
% constraint 3: N' IAC L = 0
C3 = N' * IAC * L;
% constraint 4: N' IAC O = 0
C4 = N' * IAC * O;
% solving contraints to compute IAC
S = solve(C1==0, C2==0, C3==0, C4==0);
D = S;
S.u1 = double(S.u1);
S.u2 = double(S.u2);
S.u3 = double(S.u3);
S.u4 = double(S.u4);
IACS = [S.u1, 0, S.u2;0, 1, S.u3;S.u2, S.u3, S.u4];

IACS = inv(IACS);

K = chol(IACS, 'upper')


%% Functions
function [l] = segToLine(pts)
% convert the endpoints of a line segment to a line in homogeneous
% coordinates.
%
% pts are the endpoits of the segment: [x1 y1;
%                                       x2 y2]

% convert endpoints to cartesian coordinates
a = [pts(1,:)';1];
b = [pts(2,:)';1];
l = cross(a,b);
l = l./norm(l);
end

function [lines] = computeLines(edgs, minLength)
% Compute the Hough transform of the binary image returned by edge.
[H,theta,rho] = hough(edgs);
% Display the transform, H, returned by the hough function.
figure(3)
imshow(imadjust(rescale(H)),[],...
       'XData',theta,...
       'YData',rho,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on
colormap(gca,hot)
% Find the peaks in the Hough transform matrix, H, using the houghpeaks function.
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% Superimpose a plot on the image of the transform that identifies the peaks.
x = theta(P(:,2));
y = rho(P(:,1));
plot(x,y,'s','color','black');
% Find lines in the image using the houghlines function.
% GOOD: lines = houghlines(edgs,theta,rho,P,'FillGap',70,'MinLength',350);
lines = houghlines(edgs,theta,rho,P,'FillGap',70,'MinLength',minLength);
%lines = houghlines(edgs,theta,rho,P,'FillGap',70,'MinLength',10);
end

