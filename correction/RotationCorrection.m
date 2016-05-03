function correctedImage = RotationCorrection(im)
% This function gets a rotated image and corrects its clockwise or
% counterclockwise rotation. It can be useful for enhancing the output
% images of desktop scanners. 
% 
% Usage: please refer to 'example.m' script.
%
% See also: EDGE, BWMORPH, RADON, IMROTATE 
%
% Original version by Amir Hossein Omidvarnia,  October 2007
% Email: aomidvar@ece.ut.ac.ir

%%%%% Convertion of input image into a gray scale image....
if( ~isgray(im)) 
    grayImage = rgb2gray(im);
else
    grayImage = im;
end
%%%%%

%%%%% Edge detection and edge linking....
binaryImage = edge(grayImage,'canny'); % 'Canny' edge detector
%imshow(binaryImage);

binaryImage = bwmorph(binaryImage,'thicken'); % A morphological operation for edge linking
% correctedImage = binaryImage;
% return;
%%%%%

%%%%% Radon transform projections along 180 degrees, from -90 to +89....
% R: Radon transform of the intensity image 'grayImage' for -90:89 degrees.
% In fact, each column of R shows the image profile along corresponding angle. 
% xp: a vector containing the radial coordinates corresponding to each row of 'R'.
% Negative angles correspond to clockwise directions, while positive angles
% correspond to counterclockwise directions around the center point (up-left corner).
% R1: A 1x180 vector in which, each element is equal the maximum value of Radon transform along each angle.
% This value reflects the maximum number of pixels along each direction. 
% r_max: A 1x180 vector, which includes corresponding radii of 'R1'.
theta = -90:89;
[R,xp] = radon(binaryImage,theta);
figure, imagesc(theta,xp,R);colormap(hot);colorbar;
return;
[R1,r_max] = max(R); 
display(r_max);
%%%%%

%%%%% Line detection....
% This section performs a Hough-like search. It finds maximum value of Radon
% transform over all radii and angles in angles greater than 50 or 
% less than -50. First detected angle indicates the slope of the upper bond of the image. 
theta_max = 90;
while(theta_max > 50 || theta_max<-50)
    [R2,theta_max] = max(R1); % R2: Maximum Radon transform value over all angles. 
                              % theta_max: Corresponding angle 
    R1(theta_max) = 0; % Remove element 'R2' from vector 'R1', so that other maximum values can be found.
    theta_max = theta_max - 91;
end
display(theta_max);
correctedImage = imrotate(im,-theta_max); % Rotation correction
correctedImage(correctedImage == 0) = 255; % Converts black resgions to white regions
