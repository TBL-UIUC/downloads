%% Cortical Bone Segmentation

% Author: Ida Ang and Mariana Kersh
% Email: ia267@cornell.edu and mkersh@illinois.edu
% Date Created: July 1 2017
% Date Last modified: June 12 2019

% Description: 
% This code scans a thickness map obtained in ImageJ/Fiji from the BoneJ 
% plugin in order to segment cortical/subchondral bone from trabecular.

%% Code Input and Output 

% Format for Required Input: 
% This code requires the raw computed tomography (CT) images to have been
% resliced in the plane of interest. These files can be formatted as tif
% or png files. 

% This code also requires a corresponding set of thickness maps which can
% be obtained using the BoneJ plugin of ImageJ. 

% Format of Output: 
% The corresponding segmentation for each image can be saved as either a
% MATLAB .mat file or a image file (png/tif). 

%% Important Details and Warnings
% This code requires an ancillary code, to read in ImageJ thickness maps.
% This code can be generated using MATLAB > Home (first tab) > Import Data 

% This code can stall if the incorrect bone border is found. Several
% modifications have been made using structuring element values to control
% erosion and dilation morphological operations. Sections 3, 6, and 7 can
% be removed if the user has a standard method for finding the correct
% image border

% Image sections can be uncommented for the user to visualize features such
% as the border, extrema points, and quadranting for grayscale range
% calculations. 

%% 1. Define path names where data is stored

% dataPath1 directs to the thickness map dataset from ImageJ/Fiji
dataPath1 = '/fill_here/fill_here/Thickness/';
% dataPath2 directs to the binarized, resliced dataset
dataPath2 = '/fill_here/fill_here/Binarized/';

% Define part of the file name used to store new images. This is the
% location that corticalSegmentation.m resides in
partFile = '/fill_here/fill_here/';

%% 2. User inputs and folder/file definition 
% User inputs

% File names for all data, thickness maps and binarized resliced dataset 
name = 'name_of_file';

rangeCheck = 1:10; % Number of files to access

% Folder created in directory the script resides
folder = sprintf('Change_Folder_Name%d', rangeCheck(1), rangeCheck(end));
mkdir(folder)

% End value defines how far the algorithm searches within the list of 
% borders. If the whole list of borders is checked, a false border
% including trabecular bone will likely be identified
endVal = 10;
    % Note: increasing this value increases code run-time
  
% Percentage of image to analyze at one time
% Changing the dimension depends on the ratio of bone to space in the image 
perDim = 0.1;       % 0.1 = 10 Percent
    % Note: increasing this value increases code run-time

% Structuring Element Value used in the final step to fill gaps 
strelVal = 2;
seD = strel('disk',strelVal); 

%% 3. Option: finding the bone border for segementation 
% The correct bone border is dependent on the values set in this section. 
% If the code stalls, set strelVal1 and strelVal2 to be larger values or 
% decrease strelVal3

% When running a large sequence of images, it is best to have as low a
% strelVal as possible. StrelVal2 is set in case strelVal1 is not
% sufficient to produce a continuous border

% An estimation of the pixel value used to erode/dilate original image 
strelVal1 = 10; 
strelVal2 = 20; % This setting should always be larger than strelVal1

% The results are usually better if the original border is eroded by 1-2 
% pixels, because this erosion accounts for noise on the periphery of the 
% image
actTempErode = 1; % To not use this function > 1-on 0-off
strelVal3 = 2;

% Disk operators for erosion and dilation
seD1 = strel('disk',strelVal1); 
seD2 = strel('disk',strelVal2); 
seD3 = strel('disk',strelVal3);

%% 4. Reading in data
% An ancillary script (here defined as importThicknessMap) has to be 
% created in order to read in ImageJ/FIJI BoneJ thickness maps. 

% Required Input 1: thimage is the thickness map produced from BoneJ 
% Required Input 2: origImage reads binarized dataset 

% Output 1: fileName sets up the directory and name of the image file 
% Output 2: matName sets up the name of the mat file if required

for i = rangeCheck % Analyzing the amount of images defined in rangeCheck

    if (i<10)
        num = '_000';
    elseif i >= 10 && i < 100
        num = '_00';
    elseif i >= 100 && i < 1000
        num = '_0';
    else
        num = '_';
    end
    
    % Edit num and name according to consistent naming scheme
    thImage = importThicknessMap([dataPath1 name num num2str(i) '.txt']);
    origImage = imread([dataPath2 name num num2str(i) '.png']);
    fileName = [partFile folder '/' name num num2str(i) '.png'];
    matName = [partFile folder '/' name num num2str(i) '.mat'];
    
    %% 5. Converting the original image and text image for use

    origImageBW = mat2gray(origImage); % Convert to a grayscale image
    
    % The text images are sometimes artificially padded in ImageJ - this 
    % cuts the thickness image to the dimensions of the original images
    dim = size(origImageBW); % Dimensions of the original image
    % Use the original image's dimensions to fix the text image's dimensions
    thImage = thImage(1:dim(1),1:dim(2));
 
    % Not a Number (NaN) in the thickness image is changed to 0 
    thImage(isnan(thImage)) = 0; 
    
    % Establish half of minimum non-zero thickness value as a threshold for
    % cortical bone segmentation
    nonZth = find(thImage > 0); % Find all non-zero thickness values 
    thMin = min(thImage(nonZth)); % Find the minimum non-zero value
    thMinTwo = thMin/2; % Used for cortical bone segementation 
    
    %% 6A. Option: Border
    % Use if using a different method to find an accurate border
%    outerBorder = bwboundaries(origImageBW);

    %% 6B. Option: Rectangular Fit for Border Estimation 
    % Implemented to estimate the size of the bone border for a given image
    
    % Extrema = [top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]
    Axis = regionprops(origImageBW,'Extrema');
    % Max vertical distance - top-left to bottom-left of the bone image
    Y = [Axis.Extrema(1,1), Axis.Extrema(1,2); Axis.Extrema(6,1), Axis.Extrema(6,2)];
    Y_axis = pdist(Y,'euclidean'); % Length in the y direction 
    % Max horizontal distance - left-top to right-bottom of the bone image
    X = [Axis.Extrema(4,1), Axis.Extrema(4,2); Axis.Extrema(8,1), Axis.Extrema(8,2)];
    X_axis = pdist(X, 'euclidean'); % Length in the x direction 
    estBorder = 2*(X_axis + Y_axis); % Square/rectangle estimate of the bone border 

    %% 7. Option: Find the correct outer border
    % If the border found is larger than the estimated border length,
    % erosion/dilation procedures will be done to connect gaps in the
    % border and the border length will be recalculated
      
    imDilate = imdilate(origImageBW,seD1);
    fixedImage = imerode(imDilate,seD1);
    if actTempErode == 1
         fixedImage = imerode(fixedImage,seD3);        
    end
    
    % outerBorder is the list of border values for any image. The first
    % border listed is the leftmost border value, and each following
    % value is further from the outermost border    
    outerBorder = bwboundaries(fixedImage);
    
    % If the border list is shorter than endVal (set in section 3), set 
    % endVal to the length of the border list
    if length(outerBorder) < endVal 
        endVal = length(outerBorder);
    end

    % Find the maximum border within outerBorder - limited by endVal 

    % maxBord is the size of the largest border identified
    % maxInd is the index of the largest border identified
    [maxBord, maxInd] = max(cellfun('size', outerBorder(1:endVal), 1));

    % If the maximum border is much larger or smaller than the
    % estimated border, use the second erosion/dilation value 
    if maxBord > estBorder || maxBord < estBorder/2
        % Above operations are repeated except using seD2
        imDilate = imdilate(origImageBW,seD2);
        imErode = imerode(imDilate,seD2);
        if actTempErode == 1 
             fixedImage = imerode(imErode,seD3); 
        end
        outerBorder = bwboundaries(fixedImage);

        [maxBord, maxInd] = max(cellfun('size', outerBorder, 1));        
    end

    % Stop code if after using strelVal2 the border is still incorrect
    if maxBord > estBorder
        sprintf('Increase strelVal1/strelVal2 or decrease strelVal3')
        return
    end
    
    %% 7.1 Image for Border Check
    % Image for visualization of the border and extrema points 
    
%     figure
%     imshow(thImage); hold on;
%     plot(outerBorder{maxInd}(:,2),outerBorder{maxInd}(:,1),'r','LineWidth',5);
%     hold on
%     line([Axis.Extrema(1,1),Axis.Extrema(6,1)],[Axis.Extrema(1,2),Axis.Extrema(6,2)], 'Color', 'blue', 'LineWidth', 5);
%     hold on
%     line([Axis.Extrema(8,1),Axis.Extrema(4,1)],[Axis.Extrema(8,2),Axis.Extrema(4,2)], 'Color', 'blue', 'LineWidth', 5);
%     hold off

    %% 8. Initalization and Definitions

    % Initialize the final image with zeros (black) 
    segCortical = zeros(size(origImageBW,1),size(origImageBW,2));

    % Length of accurate border
    lenBord = length(outerBorder{maxInd});

    % X and Y coordinates of every thickness value on border
    pixelInt = outerBorder{maxInd,1};
    
    %% 9. Optimization, Separation of Subchondral Coordinates of Interest
    % In order to decrease code runtime, create a list of the thickness
    % values along the border and skip any repeated thickness values that
    % occur directly after the thickness value being evaluated
    
    for p = 1:lenBord
        % pBorder = pixel on border containing thickness value 
        pBorder(p) = thImage(outerBorder{maxInd}(p,1),outerBorder{maxInd}(p,2));
    end
    
    % If the subsequent thickness value is different, border_cut will be
    % populated with the x and y coordinates of the pixel of interest
    val = 1;
    for p = 1:lenBord-1
        if pBorder(p+1) ~= pBorder(p)
            borderCut(val,1) = pixelInt(p,1);
            borderCut(val,2) = pixelInt(p,2);
            val = val + 1;
        end
    end
    
    %% 10. Quadranting border for greater range accuracy
    
    lenBordNew = length(borderCut); % Decreased border length
    
    % These values mark where the first, second... quadrant ends
    qOne = floor((1/4)*lenBordNew);
    qTwo = floor((1/2)*lenBordNew);
    qThr = floor((3/4)*lenBordNew);
    qFour = lenBordNew;
    
    % Each array holds the thickness value for each pixel along the border
    % sectioned by quadrants
    c1 = 1; c2 = 1; c3 = 1; c4 = 1; 
    for p = 1:lenBordNew
        pBorder(p) = thImage(borderCut(p,1), borderCut(p,2));
        if p <= qOne
            listOne(c1) = pBorder(p);
            c1 = c1 + 1; 
        end
        if p > qOne && p <= qTwo
            listTwo(c2) = pBorder(p);  
            c2 = c2 + 1; 
        end
        if p > qTwo && p <= qThr
            listThr(c3) = pBorder(p); 
            c3 = c3 + 1; 
        end
        if p > qThr && p <= qFour 
            listFour(c4) = pBorder(p);   
            c4 = c4 + 1; 
        end
    end
    
    % Change the thickness values that are zero to Not A Number (NaN)
    % Zeros skew the averages to lower values  
    listOne(listOne == 0) = NaN;
    listTwo(listTwo == 0) = NaN; 
    listThr(listThr == 0) = NaN; 
    listFour(listFour == 0) = NaN; 

    % Calculate the mean thickness value of each quadrant. This value will
    % be used to calculate the grayscale range of thickness values included
    % in each border pixel calculation. 
    meanFst = nanmean(listOne);
    meanSnd = nanmean(listTwo);
    meanTrd = nanmean(listThr);
    meanFth = nanmean(listFour);
    
    %% 10.1 Quadrant Figure
%     imshow(thImage);
%     hold on;
%     plot(borderCut(:,2),borderCut(:,1),'c','LineWidth',3);
%     hold on;
%     plot(borderCut(qOne,2),borderCut(qOne,1),'r.','MarkerSize',25)
%     hold on
%     plot(borderCut(qTwo,2),borderCut(qTwo,1),'r.','MarkerSize',25)
%     hold on
%     plot(borderCut(qThr,2),borderCut(qThr,1),'r.','MarkerSize',25)
%     hold on
%     plot(borderCut(qFour,2),borderCut(qFour,1),'r.','MarkerSize',25)

    %% 11. Setting up grayscale range for segmentation 
    % For every pixel along the length of the border, evaulate a grayscale
    % range of values that will be segmented out with the image. 
    
    for p = 1:length(borderCut)  % p = pixel location within border
        
        % New image to evaulate the grayscale range for each pixel 
        GSimage = zeros(size(origImageBW,1),size(origImageBW,2));
        
        % Set the grayscale range (gsRange) of images according to each
        % quadrant
        if p <= qOne
            gsRange = meanFst*thMinTwo;
        end
        if p > qOne && p <= qTwo
            gsRange = meanSnd*thMinTwo;
        end
        if p > qTwo && p <= qThr
            gsRange = meanTrd*thMinTwo;
        end
        if p > qThr && p <= qFour 
            gsRange = meanFth*thMinTwo;           
        end
        
        % Identify the lower and upper limits of the range
        lowRange = pBorder(p) - gsRange;
        highRange = pBorder(p) + gsRange;
        
        % If lowRange is negative change value to zero
        if lowRange < 0
            lowRange = 0;
        end
        
        %% 12. Code optimization 
        % This section decreased run time by using the find function on a
        % percentage value set by perDim instead of the complete image
        
        % Dimension 2 relates to the x axis/distance
        X10 = round(perDim*dim(2)); 
        % Dimension 1 relates to the y axis/distance
        Y10 = round(perDim*dim(1)); 
        
        % Define the x and y coordinates of a square around the pixel of
        % interest. The height and width of the square is perDim of the 
        % total image, unless the pixel is close to the edge of the image. 
        X10_l = borderCut(p,2) - X10; % X10_l lower/smaller x coordinate 
        X10_h = X10 + borderCut(p,2); % X10_h higher/larger x coordinate
        Y10_l = borderCut(p,1) - Y10; % Y10_l lower/smaller y coordinate 
        Y10_h = Y10 + borderCut(p,1); % Y10_h higher/largerr y coordinate 
        
        % If the pixel of interest is close to the edge of the image, and
        % the coordinate value is negative or larger than the largest
        % x-dimension or y-dimension then the pixel is replaced 
        if X10_l < 0 
            X10_l = 0;
        end
        if X10_h > dim(2)
            X10_h = dim(2);
        end
        if Y10_l < 0 
            Y10_l = 0;
        end 
        if Y10_h > dim(1)
            Y10_h = dim(1);
        end 
        
        % After cropping, pad the image to the size of the original images
        % Crop the area around the pixel of interest
        thCrop = imcrop(thImage,[X10_l+1 Y10_l+1 X10*2 Y10*2]);
        % Pad top and right sides
        thPad1 = padarray(thCrop,[Y10_l,X10_l],'pre'); 
        % Pad bottom and left sides 
        thPad2 = padarray(thPad1,[dim(1)-Y10_h,dim(2)-X10_h],'post');
        % Cut image to original image dimensions
        cutImage = thPad2(1:dim(1),1:dim(2)); 
        
        if pBorder(p) > 0 % Run if there is a border to process
            % Find pixels within the gsRange of outer border pixel of interest 
            [rTemp,cTemp] = find(cutImage > lowRange & cutImage < highRange);

            if ~isempty(rTemp) % Run if there is data to process
                % Load the thickness values into the empty image
                for d = 1:length(rTemp)
                    GSimage(rTemp(d),cTemp(d)) = pBorder(p);
                end
                
                % Stats is a list of the areas and the pixels in each area
                stats = regionprops(imbinarize(GSimage),'PixelList');
                
                %% 13. Find the cc that had the pixel of interest in it
                % cc = connected component
                
                % Looking for the pixel in border
                for s = 1:length(stats)
                    c = intersect(stats(s).PixelList,[borderCut(p,2) borderCut(p,1)],'rows');
                    if c % Break/stop once border pixel is found 
                        break
                    end
                end
                
                % Creation of final segmented image 
                for j = 1:size(stats(s).PixelList,1)
                    segCortical(stats(s).PixelList(j,2),stats(s).PixelList(j,1)) = 1;
                end
                
            end 
            clear stats GSimage c
        end
    end

    %% 14. Optional Smoothing Step
    % Dilation and erosion smooths the final image and fills gaps
    imDilFin = imdilate(segCortical,seD);
    segCortical = imerode(imDilFin,seD);   
    
    %% 15. Saving image
    imwrite(segCortical, fileName); % Save png/tif image file
    save(matName,'segCortical'); % If needed, save .mat file
end
