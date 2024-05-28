%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is meant to analyse the area of infected cells in cell 
% monolayers over time from live cell imaging data of the bacterial 
% fluorescence. The output is a binary mask of the infected cell area over 
% time and the area over time plotted and as .xlsx file. 
% Involved steps are: 
% 1. Binarize bacterial fluorescence images
% 2. Get border of infected cell area with alpha shapes
% Created: 15/11/2023 by Lara Hundsdorfer, with help from Julio César
% Sánchez Rendón and Marie Münkel.

% INPUT DATA: Imaging data of the bacterial fluorescence as multi-tif file. 
% Example: we provide an example, see "my_sample.tif" file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear
set(0, 'DefaultFigureVisible', 'on'); % Set 'off' to supress figures

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ADAPT Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timesteps   = 16;             % number of frames in multitif file
delta       = 60;             % images taken every "delta" (minutes)
fcal        = 0.28;           % calibration factor (µm/px)
Dir         = 'Example/';     % directory with imaging data

%%% Change these values for better binarization:
background  = 30;            % check background values e.g., in Fiji (influences strongly first frames, not so much towards the end)
erode       = 4;              % increase to get rid of single points/ noise
erode_end   = 4;              % increase to erode more towards later timepoints

alphaRadius = 50;             % radius used for alphaShape function
alphaHole   = 50000;          % size (px^2) of holes which will be filled for alphaShape function
Range_fluorescence = [background 200]; % For plotting of figures

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Binarize bacterial fluorescence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DirCheck    = [Dir '/CheckBinarization']; if ~exist(DirCheck, 'dir'); mkdir(DirCheck); end  %Check whether directory exists, if not, create directory
DirBin      = [Dir '/CorrectedBinary'];   if ~exist(DirBin, 'dir');   mkdir(DirBin);   end  %Check whether directory exists, if not, create directory

Binary      = cell(1,timesteps);      % prepare cell for binary
Binary_corr = cell(1,timesteps);      % prepare cell for corrected binary
Area        = cell(timesteps,3);      % prepare array for area


parfor k = 1:timesteps   % Install AddOn 'Parallel Computing Toolbox'
    disp(k)
    I1      = imread([Dir 'my_sample.tif'], k);
    I2      = I1-background; 
    I3      = double(imadjust(I2));
    level   = multithresh(I3,5);
    I4      = imquantize(I3,level(1));
    RGB     = label2rgb(I4); 
    I5      = rgb2gray(RGB);
    level5  = graythresh(I5);
    bw      = im2bw(I5,level5);
    bw2     = imerode(bw, strel('disk', round(erode/exp(k/timesteps/erode_end)))); 
    bw3     = imdilate(bw2, strel('disk', 4));
    [rows, columns] = size(bw3);               % size of the image

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Calculating area using alpha shape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [x,y]     = find(bw3);        % co-ordinates of binarized image
    shp       = alphaShape(y,-x,alphaRadius,'HoleThreshold',alphaHole); % forming an alphashape around flurescent signal
    Area{k,2} = area(shp);        % second column: area of alpha shape in (px^2)
    

    % Binarize the alpha shape:
    shp.Points(:,2)     = shp.Points(:,2)+rows;
    [qx, qy]            = meshgrid(1:rows, 1:columns);
    bw_alphaShape_noisy = inShape(shp,qx,qy);
    
    % Filter bw of alphaShape for only the largest area 
    Label_matrix        = labelmatrix(bwconncomp(bw_alphaShape_noisy));
    stats               = regionprops(Label_matrix,'FilledArea'); %Get labelled image, sizes in pixels, than filter by size
    label_index         = find([stats.FilledArea]==max([stats.FilledArea]));
    Label_matrix(Label_matrix~=label_index) = 0;
    Label_matrix(Label_matrix>0)            = 1;
    bw_alphaShape       = logical(Label_matrix);
    Binary{1,k}         = bw_alphaShape;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plotting the intermediate steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Boundary,~]        = bwboundaries(flip(bw_alphaShape),'noholes');
    Fig                 = figure;
    Fig.Position        = [100 100 1000 1000]; % [x,y,width,height]. With x and y defining the left lower border of the figure
    subplot(2,2,1); imagesc(I1);axis image; axis off; caxis(Range_fluorescence); colorbar; title('1. Input data'); 
    subplot(2,2,2); imagesc(bw3);axis image;axis off; title('2. Binarized image');
    subplot(2,2,4); plot(shp);axis([0 rows 0 columns]);title('3. Alpha shape');xlabel('pixel');ylabel('pixel');
    subplot(2,2,3); imagesc(I1);axis image; axis off; caxis(Range_fluorescence); title('4. Input data with boundary'); % colorbar;
    hold on
    for i = 1:length(Boundary)
       boundary = Boundary{i};
       plot(boundary(:,2), boundary(:,1), ':r', 'LineWidth', 3)
    end
    hold off
    sgtitle(['Timestep ' num2str(k)]);
    [ControlFigure, Map] = frame2im(getframe(gcf));
    imwrite(ControlFigure, [DirCheck '/CheckBinarization_' num2str(k) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Save area over time  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:timesteps
    Area{k,1} = k*delta/60;       % first column: time (h)
    Area{k,3} = Area{k,2}*fcal^2; % third column: area  of alpha shape in (µm^2)
end

titles        = {'Time (h)', 'Area (px^2)', 'Area (µm^2)'}; % Add column titles 
OutputFile    = [Dir 'Area.xlsx'];

% Write excel file
writecell(Area,   OutputFile,'Sheet','Area','Range','A2')
writecell(titles, OutputFile,'Sheet','Area','Range','A1')

% Plot figure
Fig = figure;
plot(horzcat(Area{1:(end-1),1}),horzcat(Area{1:(end-1),3}))
xlabel('Time (h)'); ylabel('Area (µm^2)');
saveas(gcf, [Dir 'Area.tif']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Improve binary by filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only if the pixel is also present in the area one timestep later, it is 
% part of the binary, to get rid of noise

Binary_corr{1,1}         = Binary{1,1}; 
Binary_corr{1,timesteps} = Binary{1,timesteps};
for k = 2:(timesteps-1)
    disp(k)
    Binary_corr{1,k}=Binary{1,k};
    for i = 1:length(Binary{1,k})
        for j = 1:length(Binary{1,k})
            if Binary{1,k}(i,j)==1 
                if Binary{1,k-1}(i,j)==0 && Binary{1,k+1}(i,j)==0
                    Binary_corr{1,k}(i,j)=0;
                else
                    Binary_corr{1,k}(i,j)=1;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plotting of improved binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Boundary,~] = bwboundaries(flip(Binary_corr{k}),'noholes');
    I1           = imread([Dir 'my_sample.tif'], k);
    Fig          = figure;
    imagesc(I1);axis image; axis off; caxis(Range_fluorescence); hold on
    for i = 1:length(Boundary) 
        boundary = Boundary{i}; 
        plot(boundary(:,2), boundary(:,1), ':r', 'LineWidth', 3); 
    end
    hold off
    title('Input data with improved boundary'); 
    Save_figure = getframe(Fig); imwrite(Save_figure.cdata, [DirBin '/CorrectedBinary_' num2str(k) '.png']);
end

save([Dir 'Binary.mat'], 'Binary_corr','-v7.3'); % Use version 7.3 because large file has to be saved!
close all
