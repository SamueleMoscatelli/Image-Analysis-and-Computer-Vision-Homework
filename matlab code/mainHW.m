%% Initialization
clc
clear all
close all
% true if wanted to display all the images
debug = false;
% if true constant selection of lines, if false manual selection
auto_selection = true; 
%% Load the image 
im_rgb = imread('CastelloDiMiramare.jpg');

IMG_MAX_SIZE = max(size(im_rgb));

if debug
    figure, imshow(im_rgb);
end
%% Extract the lines
% compute edges
edgs1 = compute_edges(im_rgb,1);
edgs2 = compute_edges(im_rgb,2);
edgs3 = compute_edges(im_rgb,3);
edgs4 = compute_edges(im_rgb,4);
edgs5 = compute_edges(im_rgb,5);
edgs6 = compute_edges(im_rgb,6);
edgs7 = compute_edges(im_rgb,7);
edgs8 = compute_edges(im_rgb,8);

% Finding straight lines
lines1 = computeLines(edgs1, 350, 1);
lines2 = computeLines(edgs2, 350, 1);
lines3 = computeLines(edgs3, 300, 1);
lines4 = computeLines(edgs4, 220, 1);
lines5 = computeLines(edgs5, 100, 1);
lines6 = computeLines(edgs6, 100, 1);
lines7 = computeLines(edgs7, 300, 1);
lines8 = computeLines(edgs8, 0, 2);
lines = [lines1, lines2, lines3, lines4, lines5, lines6(8), lines7(5), lines7(2), lines7(1), lines8];

% Create a plot that displays the original image with the 
% lines superimposed on it
if debug
    figure, imshow(im_rgb), hold on
    max_len = 0;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');
    end
end
%% Corner detection
% put image in gray scale
I = rgb2gray(im2double(im_rgb))./255;
% Harris detection
corners = detectHarrisFeatures(I,'MinQuality',0.1,'FilterSize', 3, 'ROI', [1,1000,3967,800]);

% create a plot that displays the original image with the 
% corners superimposed on it
figure, imshow(im_rgb); hold on;
plot(corners.selectStrongest(500));
hold off

%% Straight lines selection
% Extract parallel lines

% extraction of the lines parallel to the left part ot the plane 
[line_ind_l1, lines_l1] = select_lines(lines,im_rgb,"Select the group of long parallel lines on the left part of the plane",auto_selection,1);

% extraction of the lines perpendicular to the previous lines 
[line_ind_l2, lines_l2] = select_lines(lines,im_rgb,"Select another group of parallel lines on the left part of the plane (orthogonal to the one selected before)", auto_selection, 3);

% extraction of the lines parallel to the right part ot the plane 
[line_ind_r1, lines_r1] = select_lines(lines,im_rgb,"Select the group of long parallel lines on the right part of the plane (keep the order of the first selection)", auto_selection, 2);

% extraction of the lines perpendicular to the previous lines 
[line_ind_r2, lines_r2] = select_lines(lines,im_rgb,"Select another group of parallel lines on the right part of the plane (orthogonal to the one selected before)", auto_selection, 4);

% extraction of other lines parallel to the plane
[line_ind_other, lines_other] = select_lines(lines,im_rgb,"Select another group of parallel lines parallel to plane pi", auto_selection, 5);

% plot the selected lines
line_ind = [line_ind_l1, line_ind_l2, line_ind_r1, line_ind_r2];
plot_lines(lines(1,line_ind), im_rgb);

%% Vanishing points extraction and line at infinity fitting
% compute vanishing point of directions in the horizontal plane
% then fit the line at infinity through these points

% get vanishing points
vp_l1 = getVp(lines_l1);
vp_l2 = getVp(lines_l2);
vp_r1 = getVp(lines_r1);
vp_r2 = getVp(lines_r2);
vp_other = getVp(lines_other);

% fit the line at infinity through these points
l_inf_prime = fitLine([vp_l1 vp_l2 vp_r1 vp_r2 vp_other],false);

%% Affine rectification
% compute the affine reconstruction matrix
H_r_aff = [1 0 0; 0 1 0; l_inf_prime(1) l_inf_prime(2) l_inf_prime(3)];

% apply affine rectification matrix to the image and show it
img_affine = transform_and_show(H_r_aff, im_rgb, "Affine rectification");

%% Metric rectification (1)
% computation of two pairs of perpendicular lines
perpLines = [createLinePairsFromTwoSets(lines_l1, lines_l2), createLinePairsFromTwoSets(lines_r1, lines_r2)];

% transform lines according to the affine rectification matrix
perpLines = transformLines(H_r_aff, perpLines);

%% Metric rectification (2)
% computation of the matrix of the transformation from affine rectification
% to euclidean reconstruction through linear regression
ls = [];
ms = [];
index = 1;
for ii = 1:2:size(perpLines,2)
    ls(:, index) = perpLines(:, ii);
    ms(:, index) = perpLines(:, ii+1);
    index = index + 1;
end

% fit the transformation from affinity to euclidean
H_a_e = getH_from_affine(ls,ms);

% Application of the transformation and image showing
transform_and_show(H_a_e, img_affine, "Euclidean Reconstruction");

%% Metric rectification (3)
% Flip the image along the x axis
angle = 180;
R = rotx(deg2rad(180));

% computation of the composite transformation
H_r = R * H_a_e * H_r_aff;

% application of the transformation to the original image and 
% image showing
out_e = transform_and_show(H_r, im_rgb, "Euclidean Reconstruction");

%% Camera calibration (1)

% computation of the normalization matrix
H_scaling = diag([1/IMG_MAX_SIZE, 1/IMG_MAX_SIZE, 1]);

% normalization of the line at infinity
l_infs = H_scaling.' \ l_inf_prime;

% select the vertical lines
[line_ind_vert, lines_vertical] = select_lines(lines,im_rgb,"Select vertical parallel lines", auto_selection, 6);

% plot all the selected lines
if debug
    line_ind = [line_ind, line_ind_vert];
    plot_lines(lines(1,line_ind), im_rgb);
end

% get the vertical vanishing point
vp_vertical = H_scaling * getVp(lines_vertical);

% get the vanishing point of other lines on the plane
%vp_other = H_scaling * getVp(lines_other);

%% Camera calibration (2)

% computation of the image of the absolute conic 
% using l_inf, vertical vp and homography
IAC = get_IAC(l_infs, vp_vertical, [], [], H_scaling/H_r);

% get the intrinsic parameters
alfa = sqrt(IAC(1,1));
u0 = -IAC(1,3)/(alfa^2);
v0 = -IAC(2,3);
fy = sqrt(IAC(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy /alfa;

% computation of K
K = [fx 0 u0; 0 fy v0; 0 0 1];

% denormalization of K
K = H_scaling \ K;

% get intrinsic parameters after denormalization
fx = K(1,1);
fy = K(2,2);
u0 = K(1,3);
v0 = K(2,3);
alfa = fx/fy;

%% Metric measurements functional to camera localization
% get the liles correspondig to the four sides of the reference rectangle
line_hor_dl = cross([lines(13).point2(1);lines(13).point2(2);1],vp_l1);
line_hor_ul = cross([lines(15).point2(1);lines(15).point2(2);1],vp_l1);
line_obl_lr = cross([lines(15).point1(1);lines(15).point1(2);1],vp_l2);
line_obl_ll = cross([lines(13).point1(1);lines(13).point1(2);1],vp_l2);

% normalize them
line_hor_dl = line_hor_dl./line_hor_dl(3);
line_hor_ul = line_hor_ul./line_hor_ul(3);
line_obl_lr = line_obl_lr./line_obl_lr(3);
line_obl_ll = line_obl_ll./line_obl_ll(3);

% transform the lines according tp the metric rectification
l_left = H_r.' \ line_obl_ll;
l_up = H_r.' \ line_hor_ul;
l_down = H_r.' \ line_hor_dl;
l_right = H_r.' \ line_obl_lr;

% compute the four points of the reference rectangle
x_ul_left = cross(l_left,l_up);
x_dl_left = cross(l_left,l_down);
x_ur_left = cross(l_right,l_up);
x_dr_left = cross(l_down,l_right);

% normalization of the points
x_ul_left = x_ul_left ./ x_ul_left(3,1); 
x_dl_left = x_dl_left ./ x_dl_left(3,1);
x_ur_left = x_ur_left./ x_ur_left(3,1);
x_dr_left = x_dr_left./ x_dr_left(3,1);

% compute the lenght of the 
length_side1 = norm(x_ul_left - x_dl_left,2);
length_side2 = norm(x_ur_left - x_dr_left,2);

% do the average
length_side12_img = (length_side1 + length_side2) / 2;

% measure the length of the short side on the left part (needed for
% homography estimation
length_side3 = norm(x_ul_left - x_ur_left,2);
length_side4 = norm(x_dl_left - x_dr_left,2);

% do the average
length_side34_img = (length_side3 + length_side4) / 2;

% calculate aspect ratio, this is the aspect ratio on the real world since
% now the image is a similarity.
aspect_ratio_left = length_side12_img/length_side34_img;


%% Camera localization
% Determine the relative pose (i.e. position and orientation) 
% between the reference attached to the horizontal plane ∏ and 
% the camera reference.
% [i j o] are respectively the two planar axes ad the center of the
% reference with respect to the planar object.
% Using the camera calibration matrix and the homografy which map the image
% into the planar object it is possible to determine i, j and o, so
% localizing the planar object with reference to the camera and vice versa.

% build up the reference rectangle of the left face keeping the right
% proportion between the sides (it has been choosen 200 as measure of one
% side, because the real measures are not known)
Np = 200;
x_dl = [0 0];
x_ul = [0 Np];
x_dr = [Np/aspect_ratio_left 0];
x_ur = [Np/aspect_ratio_left Np];

% transform lines according to the metric reconstruction
l_left = H_r.' \ line_obl_ll;
l_up = H_r.' \ line_hor_ul;
l_down = H_r.' \ line_hor_dl;
l_right = H_r.' \ line_obl_lr;

% compute the four points of the reference rectangle
x_ul_left = cross(l_left,l_up);
x_dl_left = cross(l_left,l_down);
x_ur_left = cross(l_right,l_up);
x_dr_left = cross(l_down,l_right);

% normalization of the points
x_ul_left = x_ul_left ./ x_ul_left(3,1); 
x_dl_left = x_dl_left ./ x_dl_left(3,1);
x_ur_left = x_ur_left./ x_ur_left(3,1);
x_dr_left = x_dr_left./ x_dr_left(3,1);

% fit the homography from scene to image 
H_omog = fitgeotrans([x_ul; x_dl; x_ur; x_dr], [x_ul_left(1:2).'; x_dl_left(1:2).'; x_ur_left(1:2).'; x_dr_left(1:2).'], 'projective');
H_omog = H_omog.T.';

% extract columns
h1 = H_omog(:,1);
h2 = H_omog(:,2);
h3 = H_omog(:,3);

% normalization factor
lambda = 1 / norm(K \ h1);

% computation of the rotation matrix coloumns
r1 = (K \ h1) * lambda;
r2 = (K \ h2) * lambda;
r3 = cross(r1,r2);

% computation of the rotation matrix of the world with respect to the camera 
% where the world in this case is the left horizontal face
R = [r1, r2, r3];

% due to noise in the data R may be not a true rotation matrix.
% approximate it through svd, obtaining a orthogonal matrix
[U, ~, V] = svd(R);
R = U * V';

% Computation of the translation vector, which represents the position 
% of the plane with respect to the reference frame of the camera
T = (K \ (lambda * h3));

% computation of the camera rotation
cameraRotation = R.';

% since T is expressed in the camera reference frame computation of the 
% translation vector in the plane reference frame
cameraPosition = -R.'*T;
%% Display orientation and position of the camera with referenze to the left part of the plane

figure
plotCamera('Location', cameraPosition, 'Orientation', cameraRotation.', 'Size', 20);
hold on
pcshow([[x_ul; x_dl; x_ur; x_dr], zeros(size([x_ul; x_dl; x_ur; x_dr],1), 1)], ...
    'red','VerticalAxisDir', 'up', 'MarkerSize', 200);
xlabel('X')
ylabel('Y')
zlabel('Z')
%% Reconstruction of the left vertical facade
%computation of the projection matrix
P = K * [R,T]; 
% computation of the reconstruction matrix
H_vert_l_sr = inv([P(:,1), P(:,3), P(:,4)]);
R_y = roty(deg2rad(180));

H_vert_l_sr = R_y * H_vert_l_sr;

% application of the reconstruction matrix to the original image 
% and showing of the result
out_lv = transform_and_show(H_vert_l_sr, im_rgb, "Shape reconstruction of " + ...
                                        "the left vertical face");
