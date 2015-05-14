%% Infarct Tissue Heterogeneity
function myocardial()
%% Callback: Browse
    function browseButton_Callback(hObject, eventdata, handles)
        folder_name = uigetdir();
        if(folder_name == 0)
            return
        end
        h = findobj('Tag', 'loadImagesEdit');
        set(h,'String',folder_name);
    end

%% Callback: Load Dicom Images
    function loadDicomImagesButton_Callback(hObject, eventdata, handles)
        loadImagesEdit = findobj('Tag', 'loadImagesEdit');
        imageListBox = findobj('Tag', 'imageListBox');
        images_dir = get(loadImagesEdit,'String');
        images = dir(strcat(images_dir,''));
        [sorted_names,~] = sortrows({images.name}');
        sorted_names = sorted_names(3:end);
        set(imageListBox,'String',sorted_names,'Value',1);
        slicessize = length(sorted_names);
        [imgs, infos] = readAllDICOMImages(images_dir,sorted_names);
        displayImage(imgs(:,:,1));
        displayImageInfo(infos(1));
    end

%% Callback: Contrast Images
    function contrastImagesButton_Callback(hObject, eventdata, handles)
        imcontrast(imageAxes);
        I = getimage(imgca);
        assignin('base','im', I);  
    end

%% Callback: Image List Box
    function imageListBox_Callback(hObject, eventdata, handles)
        [selected_image,selected_image_info,selected_image_index] = getSelectedImage();
        displayImage(selected_image);
        displayImageInfo(selected_image_info);
    end

%% Callback: Mark Epicardial
    function markEpicardialButton_Callback(hObject, eventdata, handles)
        currentEpicardialmark = clearEpicardial();
        if ~isempty(currentEpicardialmark)
            poly = impoly(imageAxes,[currentEpicardialmark(:,1),currentEpicardialmark(:,2)]);
        else
            poly = impoly(imageAxes,Epicardialpolyposition);
        end
        Epicardialpolyposition = wait(poly);
        delete(poly);
        
        Epicardialpolyposition = [Epicardialpolyposition; Epicardialpolyposition(1,:)];
        xi = Epicardialpolyposition(:,1); yi = Epicardialpolyposition(:,2);
        selected_image_index = getSelectedImageIndex();
        zi = ones(length(xi),1) * selected_image_index;
        Epicardial = [Epicardial; [xi, yi, zi]];
        
        hold on;
        Epicardialplot = plot(imageAxes,xi,yi,Epicardialcolor,'Linewidth', lineWidth);
        
        assignin('base', 'Epicardial', Epicardialpolyposition);
    end

%% Callback: Mark Endocardial
    function markEndocardialButton_Callback(hObject, eventdata, handles)
        currentEndocardialmark = clearEpicardial();
        if ~isempty(currentEndocardialmark)
            poly = impoly(imageAxes,[currentEndocardialmark(:,1),currentEndocardialmark(:,2)]);
        else
%             poly = impoly(imageAxes,Endocardialpolyposition);
            poly = impoly(imageAxes);
        end
        Endocardialpolyposition = wait(poly);
        delete(poly);
        
        Endocardialpolyposition = [Endocardialpolyposition; Endocardialpolyposition(1,:)];
        xi = Endocardialpolyposition(:,1); yi = Endocardialpolyposition(:,2);
        selected_image_index = getSelectedImageIndex();
        zi = ones(length(xi),1) * selected_image_index;
        Endocardial = [Endocardial; [xi, yi, zi]];
        
        hold on;
        Endocardialplot = plot(imageAxes,xi,yi,Endocardialcolor,'Linewidth', lineWidth);
        
        assignin('base', 'Endocardial', Endocardialpolyposition);
    end

%% Callback: Mark Remote
    function markRemoteButton_Callback(hObject, eventdata, handles)
        currentRemotemark = clearRemote();
        if ~isempty(currentRemotemark)
            poly = impoly(imageAxes,[currentRemotemark(:,1),currentRemotemark(:,2)]);
        else
%             poly = impoly(imageAxes,Remotepolyposition);
            poly = impoly(imageAxes);
        end
        Remotepolyposition = wait(poly);
        delete(poly);
        
        Remotepolyposition = [Remotepolyposition; Remotepolyposition(1,:)];
        xi = Remotepolyposition(:,1); yi = Remotepolyposition(:,2);
        selected_image_index = getSelectedImageIndex();
        zi = ones(length(xi),1) * selected_image_index;
        Remote = [Remote; [xi, yi, zi]];
        
        hold on;
        Remoteplot = plot(imageAxes,xi,yi,Remotecolor,'Linewidth', lineWidth);

        assignin('base', 'Remote', Remotepolyposition);
    end

%% Callback: Mark Hyperenhanced
    function markHyperenhancedButton_Callback(hObject, eventdata, handles)
        currentHyperenhanced = clearHyperenhanced();
        if ~isempty(currentHyperenhanced)
            poly = impoly(imageAxes,[currentHyperenhanced(:,1),currentHyperenhanced(:,2)]);
        else
%             poly = impoly(imageAxes,Hyperenhancedpolyposition);
            poly = impoly(imageAxes);
        end
        Hyperenhancedpolyposition = wait(poly);
        delete(poly);

        Hyperenhancedpolyposition = [Hyperenhancedpolyposition; Hyperenhancedpolyposition(1,:)];
        xi = Hyperenhancedpolyposition(:,1); yi = Hyperenhancedpolyposition(:,2);
        selected_image_index = getSelectedImageIndex();
        zi = ones(length(xi),1) * selected_image_index;
        Hyperenhanced = [Hyperenhanced; [xi, yi, zi]];
        
        hold on;
        Hyperenhancedplot = plot(imageAxes,xi,yi,Hyperenhancedcolor,'Linewidth', lineWidth);
        
        assignin('base', 'Hyperenhanced', Hyperenhancedpolyposition);
        
        displayHistogram([xi,yi],getAxesImage());
    end

%% 
    function displayHistogram(points, I)
        maskedImage=getMaskedImage(points,I);
        maskedVector = maskedImage(:);
        maskedVector = maskedVector(maskedVector>0);
        axes(histogramAxes);
        imhist(maskedVector); 
    end


%% get a masked Image of region of interst marked by user
    function maskedImage = getMaskedImage(points, I)
        [r, c] = size(I);
        maskIm = poly2mask(points(:,1), points(:,2), r, c);
        I = I.*cast(maskIm, class(I));
        maskedImage = im2uint8(I);
    end

%% Calculate Statistics
    function [maxSI, meanSI, stdSI] = calculateStat(points, I)
        maskedImage = getMaskedImage(points, I);
        maskedVector = maskedImage(:);
        maxSI = max(maskedVector);
        
        maskedVector = maskedVector(maskedVector>0); %removeing zeros from vector for stats
        assignin('base', 'maskedVector', maskedVector);
        
        maxSI = max(maskedVector);
        meanSI = mean(maskedVector);
        stdSI = std(double(maskedVector));
    end

    function thresholdButton_Callback(hObject, eventdata, handles)
%         mask1 = poly2mask(Epicardialpolyposition(:,1), Epicardialpolyposition(:,2), r, c);
%         mask2 = poly2mask(Endocardialpolyposition(:,1), Endocardialpolyposition(:,2), r, c);
%         mask = mask1 - mask2;
%         maskedRing = I.*cast(mask, class(I));
%         assignin('base', 'I', maskedRing);
%         V = maskedRing(:);
%         V = V(V>0); %removeing zeros from vector for stats
%         assignin('base', 'V', V);
%         maxSIRing = max(V);     
%         assignin('base', 'maxRemote', maxRemote);

        selected_image_index =  getSelectedImageIndex();
        I = getAxesImage();
%         figure, imshow(I);
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            xi = Remoteslice(:,1); yi = Remoteslice(:,2);
%             maskedRemote = getMaskedImage([xi,yi],getAxesImage());
%             maskedRemote = maskedRemote(:);
%             maskedRemote(maskedRemote>0);
            [maxRemote,meanRemote,stdRemote] = calculateStat([xi,yi],I);
            set(findobj('Tag', 'maxRemoteEdit'),'String',maxRemote);
            set(findobj('Tag', 'meanRemoteEdit'),'String',meanRemote);
            set(findobj('Tag', 'stdRemoteEdit'),'String',stdRemote);
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            xi = Hyperenhancedslice(:,1); yi = Hyperenhancedslice(:,2);
            maskedHyperenhanced = getMaskedImage([xi,yi],getAxesImage());
            [maxHyperenhanced,meanHyperenhanced,stdHyperenhanced] = calculateStat([xi,yi],I);
            set(findobj('Tag', 'maxHyperenhancedEdit'),'String',maxHyperenhanced);
            set(findobj('Tag', 'meanHyperenhancedEdit'),'String',meanHyperenhanced);
            set(findobj('Tag', 'stdHyperenhancedEdit'),'String',stdHyperenhanced);
        end
        
        grayZoneIm(:,:,1) = im2double(I);
        grayZoneIm(:,:,2) = im2double(I);
        grayZoneIm(:,:,3) = im2double(I);

        [r, c] = size(I);
        for ii = 1:r
            for jj=1:c
                if(maskedHyperenhanced(ii, jj) > maxRemote ...
                        && maskedHyperenhanced(ii, jj) < maxHyperenhanced/2.0)
                    grayZoneIm(ii, jj, 1) = 1.0;
                    grayZoneIm(ii, jj, 2) = 0.8;
                    grayZoneIm(ii, jj, 3) = 0;
                end
            end
        end
        figure, imshow(grayZoneIm,[]);
        assignin('base', 'mr', maxRemote);
        assignin('base', 'mh', maxHyperenhanced);
        assignin('base', 'hi', maskedHyperenhanced);
        assignin('base', 'g', grayZoneIm);
    end

%% Callback: Clear Epicardial
    function clearEpicardialButton_Callback(hObject, eventdata, handles)
        clearEpicardial();
    end

%% Callback: Clear Endocardial
    function clearEndocardialButton_Callback(hObject, eventdata, handles)
        clearEndocardial(); 
    end

%% Callback: Clear Remote
    function clearRemoteButton_Callback(hObject, eventdata, handles)
        clearRemote(); 
    end

%% Callback: Clear Hyperenhanced
    function clearHyperenhancedButton_Callback(hObject, eventdata, handles)
        clearHyperenhanced();
    end

%% Callback: Save Marks as .mat file for current slice
    function saveSliceMarksButton_Callback(hObject, eventdata, handles)
        selected_image_index = getSelectedImageIndex();
        [filename,pathname] = uiputfile('*.mat');
        filepath = strcat(pathname,filesep,filename);
        slice.Epicardial = [];
        slice.Endocardial = [];
        slice.Remote = [];
        slice.Hyperenhanced = [];
        if ~isempty(Epicardial)
            Epicardialslice = Epicardial(Epicardial(:,3) == selected_image_index,:)
            slice.Epicardial = Epicardialslice;
        end
        
        if ~isempty(Endocardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            slice.Endocardial = Endocardialslice;
        end
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            slice.Remote = Remoteslice;
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            slice.Hyperenhanced = Hyperenhancedslice;
        end
        save(filepath,'slice');
    end

%% Callback: Load Marks from .mat file for current slice
    function loadSliceMarksButton_Callback(hObject, eventdata, handles)
        [filename,pathname] = uigetfile('*.mat','Select Prostate Zone Marks for Current Slice');
        if(~isempty(filename))
            filepath = strcat(pathname,filesep,filename);
            clearEpicardial();
            clearEndocardial();
            clearRemote();
            clearHyperenhanced();
            data = load(filepath);
            Epicardial = [Epicardial; data.slice.Epicardial];
            Endocardial = [Endocardial; data.slice.Endocardial];
            Remote = [Remote; data.slice.Remote];
            Hyperenhanced = [Hyperenhanced; data.slice.Hyperenhanced];
            % Now Display Marks for current slice
            selected_image_index = getSelectedImageIndex();
            displayImageMarks(selected_image_index);
        end
    end

%% Callback: Save Marks as .mat file for all slices
    function saveAllMarksButton_Callback(hObject, eventdata, handles)
        [filename,pathname] = uiputfile('*.mat');
        if(~isempty(filename))
            filepath = strcat(pathname,filesep,filename);
            slices.Epicardial = Epicardial;
            slices.Endocardial = Endocardial;
            slices.Remote = Remote;
            slices.Hyperenhanced = Hyperenhanced;
            save(filepath,'slices');
        end
    end

%% Callback: Load Marks from .mat file for all slices
    function loadAllMarksButton_Callback(hObject, eventdata, handles)
        [filename,pathname] = uigetfile('*.mat','Select Prostate Zone Marks for Current Slice');
        if(~isempty(filename))
            filepath = strcat(pathname,filesep,filename);
            data = load(filepath);
            Epicardial = data.slices.Epicardial;
            Endocardial = data.slices.Endocardial;
            Remote = data.slices.Remote;
            Hyperenhanced = data.slices.Hyperenhanced;
            % Now Display Marks for current slice
            selected_image_index = getSelectedImageIndex();
            displayImageMarks(selected_image_index);
        end
    end

%% Clear Epicardial
    function currentEpicardialmark = clearEpicardial()
        currentEpicardialmark = [];
        if ~isempty(Epicardial)
            selected_image_index = getSelectedImageIndex();
            selected_index_rows = Epicardial(:,3) == selected_image_index;
            currentEpicardialmark = Epicardial(selected_index_rows,:);
            Epicardial(selected_index_rows,:) = [];
            if ishandle(Epicardialplot)
                delete(Epicardialplot);
            end
        end
    end

%% Clear Endocardial
    function currentEndocardialmark = clearEndocardial()
        currentEndocardialmark = [];
        if ~isempty(Endocardial)
            selected_image_index = getSelectedImageIndex();
            selected_index_rows = Endocardial(:,3) == selected_image_index;
            currentEndocardialmark = Endocardial(selected_index_rows,:);
            Endocardial(selected_index_rows,:) = [];
            if ishandle(Endocardialplot)
                delete(Endocardialplot);
            end
        end
    end

%% Clear Remote
    function currentRemotemark = clearRemote()
        currentRemotemark = [];
        if ~isempty(Remote)
            selected_image_index = getSelectedImageIndex();
            selected_index_rows = Remote(:,3) == selected_image_index;
            currentRemotemark = Remote(selected_index_rows,:);
            Remote(selected_index_rows,:) = [];
            if ishandle(Remoteplot)
                delete(Remoteplot);
            end
        end
    end

%% Clear Hyperenhanced
    function currentHyperenhancedmark = clearHyperenhanced()
        currentHyperenhancedmark = [];
        if ~isempty(Hyperenhanced)
            selected_image_index = getSelectedImageIndex();
            selected_index_rows = Hyperenhanced(:,3) == selected_image_index;
            currentHyperenhancedmark = Hyperenhanced(selected_index_rows,:);
            Hyperenhanced(selected_index_rows,:) = [];
            if ishandle(Hyperenhancedplot)
                delete(Hyperenhancedplot);
            end
        end
    end

%% Get Selected Image Index
    function selected_image_index = getSelectedImageIndex()
        imageListBox = findobj('Tag', 'imageListBox');
        selected_image_index = get(imageListBox,'Value');
    end

%% Get Selected Image Path
    function selected_image_path = getSelectedImagePath()
        loadImagesEdit = findobj('Tag', 'loadImagesEdit');
        imageListBox = findobj('Tag', 'imageListBox');
        images_names_list = get(imageListBox, 'String');
        selected_image_index = get(imageListBox,'Value');
        selected_image_name = images_names_list(selected_image_index);
        images_dir = get(loadImagesEdit,'string');
        selected_image_path = strcat(images_dir,filesep,selected_image_name);
        selected_image_path = selected_image_path{:};
    end

%% Get Selected Image with its Metadata and Index
    function [img, info, selected_image_index] = getSelectedImage()
        selected_image_index = 0;
        img = []; info = [];
        if ~isempty(imgs)
            imageListBox = findobj('Tag', 'imageListBox');
            selected_image_index = get(imageListBox,'Value');
            img = imgs(:,:,selected_image_index);
            info = infos(selected_image_index);
        end
    end

%% Get Current Image from Image Axes directly
    function img = getAxesImage()
        img = [];
        if ~isempty(imgs)
            axes(imageAxes);
            img = getimage(imgca);
        end
    end


%% Display Image
    function displayImage(selected_image)
        % Default behavior is showing in default axes of current figure
        % imageAxes must be cleared every time before showing new image
        % Otherwise, large memory leak
        cla(imageAxes);
%         imshow(imadjust(selected_image), [],'InitialMagnification',100);
        axes(imageAxes);
        imshow(selected_image, []);
    end

%% Display Image Info
    function displayImageInfo(metadata)
        set(findobj('Tag', 'familyNameEdit'),'String',metadata.PatientName.FamilyName);
        set(findobj('Tag', 'patientIDEdit'),'String',metadata.PatientID);
        patientBirthDate = datenum(metadata.PatientBirthDate,inDateFormat);
        patientBirthDate = datestr(patientBirthDate,outDateFormat);
        set(findobj('Tag', 'patientBirthDateEdit'),'String',patientBirthDate);
        set(findobj('Tag', 'studyIDEdit'),'String',metadata.StudyID);
        studyDate = datenum(metadata.StudyDate,inDateFormat);
        studyDate = datestr(studyDate,outDateFormat);
        set(findobj('Tag', 'studyDateEdit'),'String',studyDate);
        set(findobj('Tag', 'sliceLocationEdit'),'String',metadata.SliceLocation);
        set(findobj('Tag', 'instanceNumberEdit'),'String',metadata.InstanceNumber);
    end

%% Display Image Marks
    function displayImageMarks(selected_image_index)
        if ~isempty(Epicardial)
            Epicardialslice = Epicardial(Epicardial(:,3) == selected_image_index,:);
            xi = Epicardialslice(:,1); yi = Epicardialslice(:,2);
            hold on;
            Epicardialplot = plot(imageAxes,xi,yi,Epicardialcolor,'Linewidth', lineWidth);
        end
        
        if ~isempty(Endocardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            xi = Endocardialslice(:,1); yi = Endocardialslice(:,2);
            hold on;
            Endocardialplot = plot(imageAxes,xi,yi,Endocardialcolor,'Linewidth', lineWidth);
        end
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            xi = Remoteslice(:,1); yi = Remoteslice(:,2);
            hold on;
            Remoteplot = plot(imageAxes,xi,yi,Remotecolor,'Linewidth', lineWidth);
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            xi = Hyperenhancedslice(:,1); yi = Hyperenhancedslice(:,2);
            hold on;
            Hyperenhancedplot = plot(imageAxes,xi,yi,Hyperenhancedcolor,'Linewidth', lineWidth);
        end
        
    end

%% Read All Images
    function [imgs, infos] = readAllDICOMImages(images_dir,images_names_list)
        for i = 1:length(images_names_list)
            image_path = strcat(images_dir,filesep,images_names_list(i));
            image_path = image_path{:};
            img = dicomread(image_path);
            imgs(:,:,i) = imadjust(img);
%             imgs(:,:,i) = img;
            infos(i) = dicominfo(image_path);
        end
    end

%% All Globla Variables and GUI Construction

clear all;
close all;
clc;

database_path = '/media/student/58801F1B801EFEE6/University/Medical Image Analysis/Project/Patient_1';

%% Global Matricies
imgs = [];
infos = [];
Epicardial = [];
Endocardial = [];
Remote = [];
Hyperenhanced = [];

maskedHyperenhanced=[];

maxRemote=0; meanRemote=0; stdRemote=0;
maxHyperenhanced=0;

Epicardialplot = [];
Endocardialplot = [];
Remoteplot = [];
Hyperenhancedplot = [];

Epicardialpolyposition = [];
Endocardialpolyposition = [];
Remotepolyposition = [];
Hyperenhancedpolyposition = [];

inDateFormat = 'yyyymmdd';
outDateFormat = 'dd/mm/yyyy';
lineWidth = 3;
Epicardialcolor = 'r';
Endocardialcolor = 'g';
Remotecolor = 'b';
Hyperenhancedcolor = 'c';
slicessize = 0;

%% Figure
wMax = 100;
hMax = 100;

hFig = figure('Visible','off','Menu','none', 'Name',' Infarct Tissue Heterogeneity Analysis Tool',...
    'Resize','on', 'Position', [0 0 1000 600]);

imageAxes = axes('Parent',hFig,'Units','normalized',...
    'position',[22/wMax 10/hMax 60/wMax 80/hMax]);

histogramAxes = axes('Parent',hFig,'Units','normalized',...
    'position',[1/wMax 5/hMax 20/wMax 28/hMax]);

%% Image Directory
dirBG = uibuttongroup('Units','Normalized','Title','Images Directory',...
    'BackgroundColor',[1 0.5 0],'Position',[22/wMax 92/hMax 60/wMax 7/hMax]);

loadImagesEdit = uicontrol('Style','edit','Parent', dirBG,'Units','normalized',...
    'String',database_path,'Tag','loadImagesEdit',...
    'Position',[1/60 1/7 40/60 5/7]);

browseButton = uicontrol('Style','pushbutton','Parent',dirBG,'Units','normalized',...
    'String','Browse',...
    'Position',[42/60 1/7 6/60 5/7],'Callback',@browseButton_Callback);

loadDicomImagesButton = uicontrol('Style','pushbutton','Parent',dirBG,'Units','normalized',...
    'String','Load Images',...
    'Position',[49/60 1/7 10/60 5/7],'Callback',@loadDicomImagesButton_Callback);


%% Files Controllers
wFilesBG = 18;
hFilesBG = 7;

filesBG = uibuttongroup('Units','Normalized','Title','Files',...
    'BackgroundColor',[1 0.5 0],'Position',[1/wMax 92/hMax wFilesBG/wMax hFilesBG/hMax]);

contrastImagesButton = uicontrol('Style','pushbutton','Parent',filesBG,'Units','normalized',...
    'String','Contrast Images',...
    'Position',[1/wFilesBG 1/hFilesBG 16/wFilesBG 5/hFilesBG],'Callback',@contrastImagesButton_Callback);

%% Image ListBox
wListBG = 18;
hListBG = 20;

imagelistBG = uibuttongroup('Units','Normalized','Title','Image List',...
    'Position',[1/wMax 71/hMax wFilesBG/wMax hListBG/hMax]);

imageListBox = uicontrol('Style','listbox','Parent',imagelistBG,'Units','normalized',...
    'BackgroundColor','white','Tag','imageListBox',...
    'Position',[1/wListBG 1/hListBG 16/wListBG 19/hListBG],'Callback',@imageListBox_Callback);

%% Original Patient Parameters Handles
wOrgBG = 18;
hOrgBG = 35;

orgParameterBG = uibuttongroup('Units','Normalized','Title','Information',...
    'Position',[1/wMax 35/hMax wOrgBG/wMax hOrgBG/hMax]);

familyNameLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Family Name',...
    'Position',[1/wOrgBG 31/hOrgBG 7/wOrgBG 3/hOrgBG]);

familyNameEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','familyNameEdit','Enable','off',...
    'Position',[9/wOrgBG 31/hOrgBG 8/wOrgBG 3/hOrgBG]);

givenNameLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Given Name',...
    'Position',[1/wOrgBG 27/hOrgBG 7/wOrgBG 3/hOrgBG]);

givenNameEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','givenNameEdit','Enable','off',...
    'Position',[9/wOrgBG 27/hOrgBG 8/wOrgBG 3/hOrgBG]);

patientIDLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','ID',...
    'Position',[1/wOrgBG 23/hOrgBG 7/wOrgBG 3/hOrgBG]);

patientIDEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','patientIDEdit','Enable','off',...
    'Position',[9/wOrgBG 23/hOrgBG 8/wOrgBG 3/hOrgBG]);

patientBirthDateLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Birthdate',...
    'Position',[1/wOrgBG 19/hOrgBG 7/wOrgBG 3/hOrgBG]);

patientBirthDateEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','patientBirthDateEdit','Enable','off',...
    'Position',[9/wOrgBG 19/hOrgBG 8/wOrgBG 3/hOrgBG]);

studyIDLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Study ID',...
    'Position',[1/wOrgBG 15/hOrgBG 7/wOrgBG 3/hOrgBG]);

studyIDEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','studyIDEdit','Enable','off',...
    'Position',[9/wOrgBG 15/hOrgBG 8/wOrgBG 3/hOrgBG]);

studyDateLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Study Date',...
    'Position',[1/wOrgBG 11/hOrgBG 7/wOrgBG 3/hOrgBG]);

studyDateEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','studyDateEdit','Enable','off',...
    'Position',[9/wOrgBG 11/hOrgBG 8/wOrgBG 3/hOrgBG]);

sliceLocationLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Slice Location',...
    'Position',[1/wOrgBG 7/hOrgBG 7/wOrgBG 3/hOrgBG]);

sliceLocationEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','sliceLocationEdit','Enable','off',...
    'Position',[9/wOrgBG 7/hOrgBG 8/wOrgBG 3/hOrgBG]);

instanceNumberLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Instance Number',...
    'Position',[1/wOrgBG 1/hOrgBG 7/wOrgBG 5/hOrgBG]);

instanceNumberEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','instanceNumberEdit','Enable','off',...
    'Position',[9/wOrgBG 1/hOrgBG 8/wOrgBG 5/hOrgBG]);

%% Prostate Zones
wZoneBG = 16;
hZoneBG = 18;

%% Epicardial
EpicardialBG = uibuttongroup('Units','Normalized','Title','Epicardial Zone',...
    'BackgroundColor','red','Position',[83/wMax 81/hMax wZoneBG/wMax hZoneBG/hMax]);

markEpicardialButton = uicontrol('Style','pushbutton','Parent',EpicardialBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 14/hZoneBG 8/wZoneBG 4/hZoneBG],'Callback',@markEpicardialButton_Callback);

clearEpicardialButton = uicontrol('Style','pushbutton','Parent',EpicardialBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 14/hZoneBG 5/wZoneBG 4/hZoneBG],'Callback',@clearEpicardialButton_Callback);

maxEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
    'String','Max SI',...
    'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);

maxEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
    'Tag','maxEpicardialEdit','Enable','off',...
    'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);

meanEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
    'String','Mean SI',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

meanEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
    'Tag','meanEpicardialEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

stdEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
    'String','STD SI',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

stdEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
    'Tag','stdEpicardialEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Endocardial
EndocardialBG = uibuttongroup('Units','Normalized','Title','Endocardial Zone',...
    'BackgroundColor','green','Position',[83/wMax 62/hMax wZoneBG/wMax hZoneBG/hMax]);

markEndocardialButton = uicontrol('Style','pushbutton','Parent',EndocardialBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 14/hZoneBG 8/wZoneBG 4/hZoneBG],'Callback',@markEndocardialButton_Callback);

clearEndocardialButton = uicontrol('Style','pushbutton','Parent',EndocardialBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 14/hZoneBG 5/wZoneBG 4/hZoneBG],'Callback',@clearEndocardialButton_Callback);

maxEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
    'String','Max SI',...
    'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);

maxEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
    'Tag','maxEndocardialEdit','Enable','off',...
    'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);

meanEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
    'String','Mean SI',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

meanEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
    'Tag','meanEndocardialEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

stdEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
    'String','STD SI',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

stdEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
    'Tag','stdEndocardialEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Remote
RemoteBG = uibuttongroup('Units','Normalized','Title','Remote Zone',...
    'BackgroundColor','blue','Position',[83/wMax 43/hMax wZoneBG/wMax hZoneBG/hMax]);

markRemoteButton = uicontrol('Style','pushbutton','Parent',RemoteBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 14/hZoneBG 8/wZoneBG 4/hZoneBG],'Callback',@markRemoteButton_Callback);

clearRemoteButton = uicontrol('Style','pushbutton','Parent',RemoteBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 14/hZoneBG 5/wZoneBG 4/hZoneBG],'Callback',@clearRemoteButton_Callback);

maxRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Max SI',...
    'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);

maxRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','maxRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);

meanRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Mean SI',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

meanRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','meanRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

stdRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','STD SI',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

stdRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','stdRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Hyperenhanced
HyperenhancedBG = uibuttongroup('Units','Normalized','Title','Hyperenhanced Zone',...
    'BackgroundColor','cyan','Position',[83/wMax 24/hMax wZoneBG/wMax hZoneBG/hMax]);

markHyperenhancedButton = uicontrol('Style','pushbutton','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 14/hZoneBG 8/wZoneBG 4/hZoneBG],'Callback',@markHyperenhancedButton_Callback);

clearHyperenhancedButton = uicontrol('Style','pushbutton','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 14/hZoneBG 5/wZoneBG 4/hZoneBG],'Callback',@clearHyperenhancedButton_Callback);

maxHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Max SI',...
    'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);

maxHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','maxHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);

meanHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Mean SI',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

meanHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','meanHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

stdHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','STD SI',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

stdHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','stdHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Save Marking
wMarkBG = 16;
hMarkBG = 22;

markBG = uibuttongroup('Units','Normalized','Title','Marking',...
    'BackgroundColor','magenta','Position',[83/wMax 1/hMax wMarkBG/wMax hMarkBG/hMax]);

saveSliceMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Save for Current Slice',...
    'Position',[1/wMarkBG 18/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@saveSliceMarksButton_Callback);

loadSliceMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Load for Current Slice',...
    'Position',[1/wMarkBG 14/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@loadSliceMarksButton_Callback);

saveAllMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Save for All Slices',...
    'Position',[1/wMarkBG 10/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@saveAllMarksButton_Callback);

loadAllMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Load for All Slices',...
    'Position',[1/wMarkBG 6/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@loadAllMarksButton_Callback);

thresholdButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Threshold',...
    'Position',[1/wMarkBG 1/hMarkBG 14/wMarkBG 4/hMarkBG],'Callback',@thresholdButton_Callback);



set(hFig, 'Visible','on');

end