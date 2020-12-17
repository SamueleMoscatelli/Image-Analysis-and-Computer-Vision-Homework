function [edges] = compute_edges(im_rgb,useGaussfilt)
img = im2double(im_rgb);
switch(useGaussfilt)
    case 1
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', [0.09 0.1], 1.4142);
    case 2
        img([1500:end],[1:end]) = imgaussfilt(img([1500:end],[1:end]),5);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', [0.09 0.1], 1.4142);
    case 3
        img = imgaussfilt(img, 7);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', [0.09 0.1], 1.4142);
    case 4
        img = histeq(img);
        img(1:end,3000:end) = imgaussfilt(img(1:end,3000:end),8);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', [0.09 0.5], 1.4142);
    case 5 
        img([1300:end],[1:end]) = imgaussfilt(img([1300:end],[1:end]),4);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', [0.09 0.5], 1.4142);
    case 6
        img = histeq(img);
        img([1500:end],[1:end]) = imgaussfilt(img([1500:end],[1:end]),8);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', 0.09, 1.4142);
    case 7
        img = histeq(img);
        img = imgaussfilt(img,8);
        I = rgb2gray(img)./255;
        edges = edge(I, 'canny', 0.09, 0.1);
    case 8
        I = double(rgb2gray(im_rgb))./255;
        [BW, th] = edge(I,'canny');
        th_modified = th.*[1.5, 2.5];
        edges = edge(I,'canny', th_modified);

end
end