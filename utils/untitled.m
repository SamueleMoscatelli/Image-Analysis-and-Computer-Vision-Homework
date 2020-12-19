clear;
close all;
img = imread('CastelloDiMiramare.jpg');
imgr = imresize(img, 0.5);
figure; imshow(img);
title('Draw two pairs of parallel segments and press enter')

fprintf('Draw parallel segments\n');

% select a first pair of segments (images of 2 parallel lines)
segment1 = drawline('Color','red');
segment2 = drawline('Color','red');


% select a second pair of segments (images of 2 parallel lines)
segment3 = drawline('Color','blue');
segment4 = drawline('Color','blue');

fprintf('Press enter to continue\n');
pause

l1 = segToLine(segment1.Position);
l2 = segToLine(segment2.Position);

m1 = segToLine(segment3.Position);
m2 = segToLine(segment4.Position);

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
line([segment1.Position(1,1),segment1.Position(2,1)],[segment1.Position(1,2),segment1.Position(2,2)],'Color','red','Linewidth',3);
line([segment2.Position(1,1),segment2.Position(2,1)],[segment2.Position(1,2),segment2.Position(2,2)],'Color','red','Linewidth',3);
line([segment3.Position(1,1),segment3.Position(2,1)],[segment3.Position(1,2),segment3.Position(2,2)],'Color','blue','Linewidth',3);
line([segment4.Position(1,1),segment4.Position(2,1)],[segment4.Position(1,2),segment4.Position(2,2)],'Color','blue','Linewidth',3);
hold off;
legend('Vanishing point 1', 'Vanishing point 2','Image of l_\infty');
set(gca,'FontSize',20)

H = [eye(2),zeros(2,1); imLinfty(:)'];
% we can check that H^-T* imLinfty is the line at infinity in its canonical
% form:

% compute the rectifying matrix
fprintf('The vanishing line is mapped to:\n');
disp(inv(H)'*imLinfty);

tform = projective2d(H');
J = imwarp(imgr,tform);

figure;
imshow(J);
%imwrite(J,'provarect.JPG');


% imgAffRect = imread('provarect.JPG');
% 
% 
% figure;
% imshow(imgAffRect);
% numConstraints = 2;
% 
% hold all;
% fprintf('Draw pairs of orthogonal segments\n');
% doAgain = 1;
% count = 1;
% A = zeros(numConstraints,3);
% % select pairs of orthogonal segments
% while (count <=numConstraints)
%     figure(gcf);
%     title('Select pairs of orthogonal segments')
%     col = 'rgbcmykwrgbcmykw';
%     segment1 = drawline('Color',col(count));
%     segment2 = drawline('Color',col(count));
% 
%     l = segToLine(segment1.Position);
%     m = segToLine(segment2.Position);
% 
%     % each pair of orthogonal lines gives rise to a constraint on s
%     % [l(1)*m(1),l(1)*m(2)+l(2)*m(1), l(2)*m(2)]*s = 0
%     % store the constraints in a matrix A
%      A(count,:) = [l(1)*m(1),l(1)*m(2)+l(2)*m(1), l(2)*m(2)];
% 
%     count = count+1;
% end
% 
% %S = [x(1) x(2); x(2) 1];
% [~,~,v] = svd(A);
% s = v(:,end); %[s11,s12,s22];
% S = [s(1),s(2); s(2),s(3)];
% imDCCP = [S,zeros(2,1); zeros(1,3)]; % the image of the circular points
% [U,D,V] = svd(S);
% A = U*sqrt(D)*V';
% H = eye(3);
% H(1,1) = A(1,1);
% H(1,2) = A(1,2);
% H(2,1) = A(2,1);
% H(2,2) = A(2,2);
% 
% Hrect = inv(H);
% Cinfty = [eye(2),zeros(2,1);zeros(1,3)];
% 
% tform = projective2d(Hrect');
% J = imwarp(imgAffRect,tform);
% 
% figure;
% imshow(J);

%% Camera calibration
y = imread('CastelloDiMiramare.jpg');
figure(10), imshow(y),title('original image');
hold all;
% extraction of vertical vanishing point 
% select a first pair of segments (images of 2 parallel lines)
segment7 = drawline('Color','red');
segment8 = drawline('Color','red');
n1 = segToLine(segment7.Position);
n2 = segToLine(segment8.Position);
% compute the vanishing point
N = cross(n1,n2);
N = N./N(3);

% select the two necessary pair of parallel lines
segment9 = drawline('Color','red');
%segment4 = drawline('Color','red');
o3 = segToLine(segment9.Position);
%o4 = segToLine(segment4.Position);
%compute vanishing point
O = cross(o3,imLinfty);
O = O./O(3);


%second pair
segment10 = drawline('Color','red');
%segment6 = drawline('Color','red');
q5 = segToLine(segment10.Position);
%q6 = segToLine(segment5.Position);
%compute vanishing point
Q = cross(q5,imLinfty);
Q = Q./Q(3);


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