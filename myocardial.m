%% Infarct Tissue Heterogeneity
function myocardial()

%% Check for license
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
    message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
    reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
    if strcmpi(reply, 'No') 
        % Exit
        return;
    end
end

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
        displayImageMarks(selected_image_index);
    end

%% Callback: Mark Epicardial
    function markEpicardialButton_Callback(hObject, eventdata, handles)
        currentEpicardialmark = clearEpicardial();
        if ~isempty(currentEpicardialmark)
            poly = impoly(imageAxes,[currentEpicardialmark(:,1),currentEpicardialmark(:,2)]);
        else
            %             poly = impoly(imageAxes,Epicardialpolyposition);
            poly = impoly(imageAxes);
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
        
        showAreaAndVolume();
        assignin('base', 'Epicardial', Epicardialpolyposition);
    end

%% Callback: Mark Endocardial
    function markEndocardialButton_Callback(hObject, eventdata, handles)
        currentEndocardialmark = clearEndocardial();
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
        
        showAreaAndVolume();
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
        
        showAreaAndVolume();
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
        
        showAreaAndVolume();
        displayHistogram([xi,yi],getAxesImage(),'Hyperenhanced');
        assignin('base', 'Hyperenhanced', Hyperenhancedpolyposition);
    end

%% get a masked Image of region of interst marked by user
    function maskedImage = getMaskedImage(points, I)
        [r, c] = size(I);
        maskIm = poly2mask(points(:,1), points(:,2), r, c);
        I = I.*cast(maskIm, class(I));
        maskedImage = im2uint8(I);
    end

%% Calculate Statistics
    function [maxSI, meanSI, stdSI, modeSI] = calculateStat(points, I)
        maskedImage = getMaskedImage(points, I);
        maskedVector = maskedImage(:);
        % max should be calculated before removing zero elements
        maxSI = max(maskedVector);
        meanSI = mean(maskedVector);
        stdSI = std(double(maskedVector));
        modeSI = mode(maskedVector);
        maskedVector = maskedVector(maskedVector>0); %removeing zeros from vector for stats
        assignin('base', 'maskedVector', maskedVector);
        
        if ~isempty(maskedVector)
            maxSI = max(maskedVector);
            meanSI = mean(maskedVector);
            stdSI = std(double(maskedVector));
            modeSI = mode(maskedVector);
        end
    end
%% Callback: Calculate Infarct Core zone and Gray zone using Remote Max
    function computerGrayzoneButton_Callback(hObject, eventdata, handles)
        selected_image_index =  getSelectedImageIndex();
        I = getAxesImage();
        [r, c] = size(I);
        %                 figure, imshow(I,[]);
        
        if ~isempty(Epicardial)
            Epicardialslice = Epicardial(Epicardial(:,3) == selected_image_index,:);
            xi = Epicardialslice(:,1); yi = Epicardialslice(:,2);
            maskedEpicardial =  getMaskedImage([xi,yi], I);
            %             figure, imshow(maskedEpicardial,[]);
        end
        
        if ~isempty(Endocardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            xi = Endocardialslice(:,1); yi = Endocardialslice(:,2);
            maskedEndocardial =  getMaskedImage([xi,yi], I);
        end
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            xi = Remoteslice(:,1); yi = Remoteslice(:,2);
            %             maskedRemote = maskedRemote(:);
            %             maskedRemote(maskedRemote>0);
            %             maskedRemote = getMaskedImage([xi,yi],I);
            %             figure, imshow(maskedRemote,[]);
            [maxRemote,meanRemote,stdRemote, modeRemote] = calculateStat([xi,yi],I);          
            set(findobj('Tag', 'maxRemoteEdit'),'String',maxRemote);
            set(findobj('Tag', 'meanRemoteEdit'),'String',meanRemote);
            set(findobj('Tag', 'stdRemoteEdit'),'String',modeRemote);
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            xi = Hyperenhancedslice(:,1); yi = Hyperenhancedslice(:,2);
            %             maskedHyperenhanced = getMaskedImage([xi,yi],I);
            [maxHyperenhanced,meanHyperenhanced,stdHyperenhanced,~] = calculateStat([xi,yi],I);
            set(findobj('Tag', 'maxHyperenhancedEdit'),'String',maxHyperenhanced);
            set(findobj('Tag', 'meanHyperenhancedEdit'),'String',meanHyperenhanced);
            set(findobj('Tag', 'stdHyperenhancedEdit'),'String',stdHyperenhanced);
        end
        
        maskedRing = maskedEpicardial - maskedEndocardial;
        
        infarctCorePixelCount = 0;
        infarctGrayzonePixelCount = 0;
        remotezonePixelCount = 0;
        
        grayZoneIm(:,:,1) = im2double(I);
        grayZoneIm(:,:,2) = im2double(I);
        grayZoneIm(:,:,3) = im2double(I);
        for ii = 1:r
            for jj=1:c
                if(maskedRing(ii, jj) > maxHyperenhanced/2.0)
                    grayZoneIm(ii, jj, 1) = 1.0;
                    grayZoneIm(ii, jj, 2) = 0;
                    grayZoneIm(ii, jj, 3) = 0;
                    infarctCorePixelCount = infarctCorePixelCount + 1;                
                elseif(maskedRing(ii, jj) > maxRemote ...
                        && maskedRing(ii, jj) < maxHyperenhanced/2.0)
                    grayZoneIm(ii, jj, 1) = 1.0;
                    grayZoneIm(ii, jj, 2) = 0.8;
                    grayZoneIm(ii, jj, 3) = 0;
                    infarctGrayzonePixelCount = infarctGrayzonePixelCount + 1;
                else
                    remotezonePixelCount = remotezonePixelCount + 1;
                end
            end
        end
        figure, imshow(grayZoneIm,[]);
        showTotalAreaAndTotalVolume(infarctCorePixelCount, infarctGrayzonePixelCount, remotezonePixelCount);
    end

%% Callback: Calculate Infarct Core zone and Gray zone using Remote STD
    function computerGrayzoneWithRemoteSTDButton_Callback(hObject, eventdata, handles)
        remoteSTDFactorMin = str2double(get(findobj('Tag', 'remoteSTDFactorMinEdit'),'String'));
        remoteSTDFactorMax = str2double(get(findobj('Tag', 'remoteSTDFactorMaxEdit'),'String'));
        if ~(remoteSTDFactorMin >= 2 && remoteSTDFactorMax >= 3)
            msgbox('Remote STD Factor must be an integer value greater than 2');
            return
        end
        
        selected_image_index =  getSelectedImageIndex();
        I = getAxesImage();
        [r, c] = size(I);
        %                 figure, imshow(I,[]);
        
        if ~isempty(Epicardial)
            Epicardialslice = Epicardial(Epicardial(:,3) == selected_image_index,:);
            xi = Epicardialslice(:,1); yi = Epicardialslice(:,2);
            maskedEpicardial =  getMaskedImage([xi,yi], I);
            %             figure, imshow(maskedEpicardial,[]);
        end
        
        if ~isempty(Endocardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            xi = Endocardialslice(:,1); yi = Endocardialslice(:,2);
            maskedEndocardial =  getMaskedImage([xi,yi], I);
        end
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            xi = Remoteslice(:,1); yi = Remoteslice(:,2);
            %             maskedRemote = maskedRemote(:);
            %             maskedRemote(maskedRemote>0);
            %             maskedRemote = getMaskedImage([xi,yi],I);
            %             figure, imshow(maskedRemote,[]);
            [maxRemote,meanRemote,stdRemote, modeRemote] = calculateStat([xi,yi],I);          
            set(findobj('Tag', 'maxRemoteEdit'),'String',maxRemote);
            set(findobj('Tag', 'meanRemoteEdit'),'String',meanRemote);
            set(findobj('Tag', 'stdRemoteEdit'),'String',modeRemote);
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            xi = Hyperenhancedslice(:,1); yi = Hyperenhancedslice(:,2);
            %             maskedHyperenhanced = getMaskedImage([xi,yi],I);
            [maxHyperenhanced,meanHyperenhanced,stdHyperenhanced,~] = calculateStat([xi,yi],I);
            set(findobj('Tag', 'maxHyperenhancedEdit'),'String',maxHyperenhanced);
            set(findobj('Tag', 'meanHyperenhancedEdit'),'String',meanHyperenhanced);
            set(findobj('Tag', 'stdHyperenhancedEdit'),'String',stdHyperenhanced);
        end
        
        maskedRing = maskedEpicardial - maskedEndocardial;
        
        grayZoneIm(:,:,1) = im2double(I);
        grayZoneIm(:,:,2) = im2double(I);
        grayZoneIm(:,:,3) = im2double(I);
        
        infarctCorePixelCount = 0;
        infarctGrayzonePixelCount = 0;
        remotezonePixelCount = 0;
        
        grayZoneLowerValue = meanRemote + remoteSTDFactorMin*stdRemote;
        grayZoneUpperValue = meanRemote + remoteSTDFactorMax*stdRemote;
        
        for ii = 1:r
            for jj=1:c
                if(maskedRing(ii, jj) > grayZoneUpperValue)
                    % Infarct Core
                    grayZoneIm(ii, jj, 1) = 1.0;
                    grayZoneIm(ii, jj, 2) = 0;
                    grayZoneIm(ii, jj, 3) = 0;
                    infarctCorePixelCount = infarctCorePixelCount + 1;
                elseif(maskedRing(ii, jj) > grayZoneLowerValue ...
                        && maskedRing(ii, jj) < grayZoneUpperValue)
                     % Infarct Grayzone
                    grayZoneIm(ii, jj, 1) = 1.0;
                    grayZoneIm(ii, jj, 2) = 0.8;
                    grayZoneIm(ii, jj, 3) = 0;
                    infarctGrayzonePixelCount = infarctGrayzonePixelCount + 1;
                else
                    remotezonePixelCount = remotezonePixelCount + 1;
                end
            end
        end
        figure, imshow(grayZoneIm,[]);
        showTotalAreaAndTotalVolume(infarctCorePixelCount, infarctGrayzonePixelCount, remotezonePixelCount);
    end

%% Callback: Clear Epicardial
    function clearEpicardialButton_Callback(hObject, eventdata, handles)
        clearEpicardial();
        showAreaAndVolume();
    end

%% Callback: Clear Endocardial
    function clearEndocardialButton_Callback(hObject, eventdata, handles)
        clearEndocardial();
        showAreaAndVolume();
    end

%% Callback: Clear Remote
    function clearRemoteButton_Callback(hObject, eventdata, handles)
        clearRemote();
        showAreaAndVolume();
    end

%% Callback: Clear Hyperenhanced
    function clearHyperenhancedButton_Callback(hObject, eventdata, handles)
        clearHyperenhanced();
        showAreaAndVolume();
    end

%% Callback: Save Marks as .mat file for current slice
    function saveSliceMarksButton_Callback(hObject, eventdata, handles)
        selected_image_index = getSelectedImageIndex();
        [filename,pathname] = uiputfile('*.mat');
        
        if(filename == 0)
            return
        end
        
        filepath = strcat(pathname,filesep,filename);
        slice.Epicardial = [];
        slice.Endocardial = [];
        slice.Remote = [];
        slice.Hyperenhanced = [];
        if ~isempty(Epicardial)
            Epicardialslice = Epicardial(Epicardial(:,3) == selected_image_index,:);
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
        [filename,pathname] = uigetfile('*.mat','Select Myocardial Zone Marks for Current Slice');
        if(filename == 0)
            return
        end
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
        displayImageMarks(selected_image_index)
    end

%% Callback: Save Marks as .mat file for all slices
    function saveAllMarksButton_Callback(hObject, eventdata, handles)
        [filename,pathname] = uiputfile('*.mat');
        if(filename == 0)
            return
        end
        filepath = strcat(pathname,filesep,filename);
        slices.Epicardial = Epicardial;
        slices.Endocardial = Endocardial;
        slices.Remote = Remote;
        slices.Hyperenhanced = Hyperenhanced;
        save(filepath,'slices');
    end

%% Callback: Load Marks from .mat file for all slices
    function loadAllMarksButton_Callback(hObject, eventdata, handles)
        [filename,pathname] = uigetfile('*.mat','Select Myocardial Marks for Current Slice');
        if(filename == 0)
            return
        end
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

%%  Display Histogram of selected mask
    function displayHistogram(points, I,axesLegend)
        cla(histogramAxes,'reset');
        %         reset(histogramAxes);
        axes(histogramAxes);
        maskedImage=getMaskedImage(points,I);
        maskedVector = maskedImage(:);
        maskedVector = maskedVector(maskedVector>0);
        class(maskedVector)
        max(maskedVector)
        if ~isempty(maskedVector)
            hist(double(maskedVector),255);
            legend(histogramAxes,axesLegend,'Location','northeast');
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
            displayHistogram([xi,yi],getAxesImage(),'Hyperenhanced');
        end
        
        showAreaAndVolume();
    end

%% Read All Images
    function [imgs, infos] = readAllDICOMImages(images_dir,images_names_list)
        Epicardial = [];
        Endocardial = [];
        Remote = [];
        Hyperenhanced = [];
        for i = 1:length(images_names_list)
            image_path = strcat(images_dir,filesep,images_names_list(i));
            image_path = image_path{:};
            img = dicomread(image_path);
            imgs(:,:,i) = imadjust(img);
            %             imgs(:,:,i) = img;
            infos(i) = dicominfo(image_path);
        end
    end

%% Calculate and Show Area and Volume for automatic segmented zones
    function showTotalAreaAndTotalVolume(infarctCorePixelCount, infarctGrayzonePixelCount, remotezonePixelCount)
        [~, info, selected_image_index] = getSelectedImage();
        
        if selected_image_index == 0
            return
        end
        
        PixelSpacing = info.PixelSpacing; 
        SliceThickness = info.SliceThickness;
        
        infarctCorezoneTotalArea = infarctCorePixelCount * PixelSpacing(1) * PixelSpacing(2);
        infarctCorezoneTotalVolume = infarctCorezoneTotalArea * SliceThickness;
        set(findobj('Tag', 'areaTotalHyperenhancedEdit'),'String',infarctCorezoneTotalArea);
        set(findobj('Tag', 'volumeTotalHyperenhancedEdit'),'String',infarctCorezoneTotalVolume);
        
        infarctGrayzoneTotalArea = infarctGrayzonePixelCount * PixelSpacing(1) * PixelSpacing(2);
        infarctGrayzoneTotalVolume = infarctGrayzoneTotalArea * SliceThickness;
        set(findobj('Tag', 'areaGrayzoneEdit'),'String',infarctGrayzoneTotalArea);
        set(findobj('Tag', 'volumeGrayzoneEdit'),'String',infarctGrayzoneTotalVolume);
        
        remotezoneTotalArea = remotezonePixelCount * PixelSpacing(1) * PixelSpacing(2);
        remotezoneTotalVolume = remotezoneTotalArea * SliceThickness;
        set(findobj('Tag', 'areaTotalRemoteEdit'),'String',remotezoneTotalArea);
        set(findobj('Tag', 'volumeTotalRemoteEdit'),'String',remotezoneTotalVolume);
        
    end

%% Calculate and Show Area and Volume for manually selected zones
    function showAreaAndVolume()
        [~, info, selected_image_index] = getSelectedImage();
        
        if selected_image_index == 0
            return
        end
        
        PixelSpacing = info.PixelSpacing; 
        SliceThickness = info.SliceThickness;
        
        set(findobj('Tag', 'areaEndocardialEdit'),'String',0);
        set(findobj('Tag', 'volumeEndocardialEdit'),'String',0);
        set(findobj('Tag', 'areaEndocardialEdit'),'String',0);
        set(findobj('Tag', 'volumeEndocardialEdit'),'String',0);
        set(findobj('Tag', 'areaRemoteEdit'),'String',0);
        set(findobj('Tag', 'volumeRemoteEdit'),'String',0);
        set(findobj('Tag', 'areaHyperenhancedEdit'),'String',0);
        set(findobj('Tag', 'volumeHyperenhancedEdit'),'String',0);
        
        if ~isempty(Epicardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            xi = Endocardialslice(:,1); yi = Endocardialslice(:,2);
            areaEndocardial = polyarea(xi,yi) * PixelSpacing(1) * PixelSpacing(2);
            volumeEndocardial = areaEndocardial * SliceThickness;
            set(findobj('Tag', 'areaEpicardialEdit'),'String',areaEndocardial);
            set(findobj('Tag', 'volumeEpicardialEdit'),'String',volumeEndocardial);
        end
        
        if ~isempty(Endocardial)
            Endocardialslice = Endocardial(Endocardial(:,3) == selected_image_index,:);
            xi = Endocardialslice(:,1); yi = Endocardialslice(:,2);
            areaEndocardial = polyarea(xi,yi) * PixelSpacing(1) * PixelSpacing(2);
            volumeEndocardial = areaEndocardial * SliceThickness;
            set(findobj('Tag', 'areaEndocardialEdit'),'String',areaEndocardial);
            set(findobj('Tag', 'volumeEndocardialEdit'),'String',volumeEndocardial);
        end
        
        if ~isempty(Remote)
            Remoteslice = Remote(Remote(:,3) == selected_image_index,:);
            xi = Remoteslice(:,1); yi = Remoteslice(:,2);
            areaRemote = polyarea(xi,yi) * PixelSpacing(1) * PixelSpacing(2);
            volumeRemote = areaRemote * SliceThickness;
            set(findobj('Tag', 'areaRemoteEdit'),'String',areaRemote);
            set(findobj('Tag', 'volumeRemoteEdit'),'String',volumeRemote);
        end
        
        if ~isempty(Hyperenhanced)
            Hyperenhancedslice = Hyperenhanced(Hyperenhanced(:,3) == selected_image_index,:);
            xi = Hyperenhancedslice(:,1); yi = Hyperenhancedslice(:,2);
            areaHyperenhanced = polyarea(xi,yi) * PixelSpacing(1) * PixelSpacing(2);
            
            volumeHyperenhanced = areaHyperenhanced * SliceThickness;
            set(findobj('Tag', 'areaHyperenhancedEdit'),'String',areaHyperenhanced);
            set(findobj('Tag', 'volumeHyperenhancedEdit'),'String',volumeHyperenhanced);
        end
    end

%% All Global Variables and GUI Construction

clear all;
close all;
clc;

% database_path = '/media/student/58801F1B801EFEE6/University/Medical Image Analysis/Project/Patient_1';
database_path = 'D:\University\Medical Image Analysis\Project\Patient_1';

%% Global Matricies
imgs = [];
infos = [];
Epicardial = [];
Endocardial = [];
Remote = [];
Hyperenhanced = [];

maxRemote=0;
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
lineWidth = 1;
Epicardialcolor = 'c';
Endocardialcolor = 'g';
Remotecolor = 'b';
Hyperenhancedcolor = 'r';
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
hFilesBG = 20;
hFilesBGR = 18;

filesBG = uibuttongroup('Units','Normalized','Title','Files',...
    'BackgroundColor',[1 0.5 0],'Position',[1/wMax 82/hMax wFilesBG/wMax hFilesBGR/hMax]);

contrastImagesButton = uicontrol('Style','pushbutton','Parent',filesBG,'Units','normalized',...
    'String','Contrast Images','FontSize',9,...
    'Position',[1/wFilesBG 17/hFilesBG 16/wFilesBG 3/hFilesBG],'Callback',@contrastImagesButton_Callback);

computerGrayzoneButton = uicontrol('Style','pushbutton','Parent',filesBG,'Units','normalized',...
    'String','Grayzone with Remote Max','FontSize',9,...
    'Position',[1/wFilesBG 13/hFilesBG 16/wFilesBG 3/hFilesBG],'Callback',@computerGrayzoneButton_Callback);

computerGrayzoneWithRemoteSTDButton = uicontrol('Style','pushbutton','Parent',filesBG,'Units','normalized',...
    'String','Grayzone with Remote STD','FontSize',9,...
    'Position',[1/wFilesBG 9/hFilesBG 16/wFilesBG 3/hFilesBG],'Callback',@computerGrayzoneWithRemoteSTDButton_Callback);

remoteSTDFactorMinLabel = uicontrol('Style','text','Parent',filesBG,'Units','normalized',...
    'String','Remote STD Factor Min','FontSize',8,...
    'Position',[1/wFilesBG 5/hFilesBG 12.5/wFilesBG 3/hFilesBG]);

remoteSTDFactorMinEdit = uicontrol('Style','edit','Parent',filesBG,'Units','normalized',...
    'Tag','remoteSTDFactorMinEdit','FontSize',7,...
    'Position',[14/wFilesBG 5/hFilesBG 3/wFilesBG 3/hFilesBG]);

remoteSTDFactorMaxLabel = uicontrol('Style','text','Parent',filesBG,'Units','normalized',...
    'String','Remote STD Factor Max','FontSize',8,...
    'Position',[1/wFilesBG 1/hFilesBG 12.5/wFilesBG 3/hFilesBG]);

remoteSTDFactorMaxEdit = uicontrol('Style','edit','Parent',filesBG,'Units','normalized',...
    'Tag','remoteSTDFactorMaxEdit','FontSize',7,...
    'Position',[14/wFilesBG 1/hFilesBG 3/wFilesBG 3/hFilesBG]);

%% Image ListBox
wListBG = 18;
hListBG = 20;
hListBGR = 15;

imagelistBG = uibuttongroup('Units','Normalized','Title','Image List',...
    'Position',[1/wMax 67/hMax wFilesBG/wMax hListBGR/hMax]);

imageListBox = uicontrol('Style','listbox','Parent',imagelistBG,'Units','normalized',...
    'BackgroundColor','white','Tag','imageListBox',...
    'Position',[1/wListBG 1/hListBG 16/wListBG 19/hListBG],'Callback',@imageListBox_Callback);

%% Original Patient Parameters Handles
wOrgBG = 18;
hOrgBG = 32;

orgParameterBG = uibuttongroup('Units','Normalized','Title','Information',...
    'Position',[1/wMax 34/hMax wOrgBG/wMax hOrgBG/hMax]);

familyNameLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Family Name','FontSize',8,...
    'Position',[1/wOrgBG 29/hOrgBG 7/wOrgBG 3/hOrgBG]);

familyNameEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','familyNameEdit','Enable','off',...
    'Position',[9/wOrgBG 29/hOrgBG 8/wOrgBG 3/hOrgBG]);

givenNameLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Given Name','FontSize',8,...
    'Position',[1/wOrgBG 25/hOrgBG 7/wOrgBG 3/hOrgBG]);

givenNameEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','givenNameEdit','Enable','off',...
    'Position',[9/wOrgBG 25/hOrgBG 8/wOrgBG 3/hOrgBG]);

patientIDLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','ID',...
    'Position',[1/wOrgBG 21/hOrgBG 7/wOrgBG 3/hOrgBG]);

patientIDEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','patientIDEdit','Enable','off','FontSize',7,...
    'Position',[9/wOrgBG 21/hOrgBG 8/wOrgBG 3/hOrgBG]);

patientBirthDateLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Birthdate',...
    'Position',[1/wOrgBG 17/hOrgBG 7/wOrgBG 3/hOrgBG]);

patientBirthDateEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','patientBirthDateEdit','Enable','off','FontSize',9,...
    'Position',[9/wOrgBG 17/hOrgBG 8/wOrgBG 3/hOrgBG]);

studyIDLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Study ID',...
    'Position',[1/wOrgBG 13/hOrgBG 7/wOrgBG 3/hOrgBG]);

studyIDEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','studyIDEdit','Enable','off','FontSize',7,...
    'Position',[9/wOrgBG 13/hOrgBG 8/wOrgBG 3/hOrgBG]);

studyDateLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Study Date',...
    'Position',[1/wOrgBG 9/hOrgBG 7/wOrgBG 3/hOrgBG]);

studyDateEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','studyDateEdit','Enable','off','FontSize',9,...
    'Position',[9/wOrgBG 9/hOrgBG 8/wOrgBG 3/hOrgBG]);

sliceLocationLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Slice Location',...
    'Position',[1/wOrgBG 5/hOrgBG 7/wOrgBG 3/hOrgBG]);

sliceLocationEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','sliceLocationEdit','Enable','off',...
    'Position',[9/wOrgBG 5/hOrgBG 8/wOrgBG 3/hOrgBG]);

instanceNumberLabel = uicontrol('Style','text','Parent',orgParameterBG,'Units','normalized',...
    'String','Inst. Num',...
    'Position',[1/wOrgBG 1/hOrgBG 7/wOrgBG 3/hOrgBG]);

instanceNumberEdit = uicontrol('Style','edit','Parent',orgParameterBG,'Units','normalized',...
    'Tag','instanceNumberEdit','Enable','off',...
    'Position',[9/wOrgBG 1/hOrgBG 8/wOrgBG 3/hOrgBG]);

%% Myocardial Zones
%% Epicardial
wZoneBG = 16;
hZoneBG = 12;

EpicardialBG = uibuttongroup('Units','Normalized','Title','Epicardial Zone',...
    'BackgroundColor',Epicardialcolor,'Position',[83/wMax 88/hMax wZoneBG/wMax hZoneBG/hMax]);

markEpicardialButton = uicontrol('Style','pushbutton','Parent',EpicardialBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 9/hZoneBG 8/wZoneBG 3/hZoneBG],'Callback',@markEpicardialButton_Callback);

clearEpicardialButton = uicontrol('Style','pushbutton','Parent',EpicardialBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 9/hZoneBG 5/wZoneBG 3/hZoneBG],'Callback',@clearEpicardialButton_Callback);

areaEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
    'String','Area',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
    'Tag','areaEpicardialEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
    'String','Volume',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
    'Tag','volumeEpicardialEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

% maxEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
%     'String','Max SI',...
%     'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% maxEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
%     'Tag','maxEpicardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);
% 
% meanEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
%     'String','Mean SI',...
%     'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% meanEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
%     'Tag','meanEpicardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);
% 
% stdEpicardialLabel = uicontrol('Style','text','Parent',EpicardialBG,'Units','normalized',...
%     'String','STD SI',...
%     'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% stdEpicardialEdit = uicontrol('Style','edit','Parent',EpicardialBG,'Units','normalized',...
%     'Tag','stdEpicardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Endocardial
EndocardialBG = uibuttongroup('Units','Normalized','Title','Endocardial Zone',...
    'BackgroundColor',Endocardialcolor,'Position',[83/wMax 76/hMax wZoneBG/wMax hZoneBG/hMax]);

markEndocardialButton = uicontrol('Style','pushbutton','Parent',EndocardialBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 9/hZoneBG 8/wZoneBG 3/hZoneBG],'Callback',@markEndocardialButton_Callback);

clearEndocardialButton = uicontrol('Style','pushbutton','Parent',EndocardialBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 9/hZoneBG 5/wZoneBG 3/hZoneBG],'Callback',@clearEndocardialButton_Callback);

areaEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
    'String','Area',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
    'Tag','areaEndocardialEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
    'String','Volume',...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
    'Tag','volumeEndocardialEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

% maxEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
%     'String','Max SI',...
%     'Position',[1/wZoneBG 9/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% maxEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
%     'Tag','maxEndocardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 9/hZoneBG 7/wZoneBG 3/hZoneBG]);
% 
% meanEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
%     'String','Mean SI',...
%     'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% meanEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
%     'Tag','meanEndocardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);
% 
% stdEndocardialLabel = uicontrol('Style','text','Parent',EndocardialBG,'Units','normalized',...
%     'String','STD SI',...
%     'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);
% 
% stdEndocardialEdit = uicontrol('Style','edit','Parent',EndocardialBG,'Units','normalized',...
%     'Tag','stdEndocardialEdit','Enable','off',...
%     'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Remote
wZoneBG = 16;
hZoneBG = 33;

RemoteBG = uibuttongroup('Units','Normalized','Title','Remote Zone',...
    'BackgroundColor',Remotecolor,'Position',[83/wMax 51/hMax wZoneBG/wMax 25/hMax]);

markRemoteButton = uicontrol('Style','pushbutton','Parent',RemoteBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 29/hZoneBG 8/wZoneBG 3/hZoneBG],'Callback',@markRemoteButton_Callback);

clearRemoteButton = uicontrol('Style','pushbutton','Parent',RemoteBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 29/hZoneBG 5/wZoneBG 3/hZoneBG],'Callback',@clearRemoteButton_Callback);

areaRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Area',...
    'Position',[1/wZoneBG 25/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','areaRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 25/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Volume',...
    'Position',[1/wZoneBG 21/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','volumeRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 21/hZoneBG 7/wZoneBG 3/hZoneBG]);

areaTotalRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Total Area',...
    'Position',[1/wZoneBG 17/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaTotalRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','areaTotalRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 17/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeTotalRemoteLabel = uicontrol('Style','text','Parent',RemoteBG,'Units','normalized',...
    'String','Total Volume','FontSize',7,...
    'Position',[1/wZoneBG 13/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeTotalRemoteEdit = uicontrol('Style','edit','Parent',RemoteBG,'Units','normalized',...
    'Tag','volumeTotalRemoteEdit','Enable','off',...
    'Position',[8/wZoneBG 13/hZoneBG 7/wZoneBG 3/hZoneBG]);

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
    'BackgroundColor',Hyperenhancedcolor,'Position',[83/wMax 24/hMax wZoneBG/wMax 27/hMax]);

markHyperenhancedButton = uicontrol('Style','pushbutton','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Mark (Edit)',...
    'Position',[1/wZoneBG 29/hZoneBG 8/wZoneBG 3/hZoneBG],'Callback',@markHyperenhancedButton_Callback);

clearHyperenhancedButton = uicontrol('Style','pushbutton','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Clear',...
    'Position',[10/wZoneBG 29/hZoneBG 5/wZoneBG 3/hZoneBG],'Callback',@clearHyperenhancedButton_Callback);

areaHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Area',...
    'Position',[1/wZoneBG 25/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','areaHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 25/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Volume',...
    'Position',[1/wZoneBG 21/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','volumeHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 21/hZoneBG 7/wZoneBG 3/hZoneBG]);

areaTotalHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Total Area',...
    'Position',[1/wZoneBG 17/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaTotalHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','areaTotalHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 17/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeTotalHyperenhancedLabel = uicontrol('Style','text','Parent',HyperenhancedBG,'Units','normalized',...
    'String','Total Volume','FontSize',7,...
    'Position',[1/wZoneBG 13/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeTotalHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','volumeTotalHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 13/hZoneBG 7/wZoneBG 3/hZoneBG]);

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
    'String','STD SI','FontSize',8,...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

stdHyperenhancedEdit = uicontrol('Style','edit','Parent',HyperenhancedBG,'Units','normalized',...
    'Tag','stdHyperenhancedEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Grayzone
wZoneBG = 16;
hZoneBG = 8;

GrayzoneBG = uibuttongroup('Units','Normalized','Title','Grayzone',...
    'BackgroundColor',[1, 1, 0],'Position',[83/wMax 16/hMax wZoneBG/wMax hZoneBG/hMax]);

areaGrayzoneLabel = uicontrol('Style','text','Parent',GrayzoneBG,'Units','normalized',...
    'String','Total Area',...
    'Position',[1/wZoneBG 5/hZoneBG 6/wZoneBG 3/hZoneBG]);

areaGrayzoneEdit = uicontrol('Style','edit','Parent',GrayzoneBG,'Units','normalized',...
    'Tag','areaGrayzoneEdit','Enable','off',...
    'Position',[8/wZoneBG 5/hZoneBG 7/wZoneBG 3/hZoneBG]);

volumeGrayzoneLabel = uicontrol('Style','text','Parent',GrayzoneBG,'Units','normalized',...
    'String','Total Volume','FontSize',7,...
    'Position',[1/wZoneBG 1/hZoneBG 6/wZoneBG 3/hZoneBG]);

volumeGrayzoneEdit = uicontrol('Style','edit','Parent',GrayzoneBG,'Units','normalized',...
    'Tag','volumeGrayzoneEdit','Enable','off',...
    'Position',[8/wZoneBG 1/hZoneBG 7/wZoneBG 3/hZoneBG]);

%% Save Marking
wMarkBG = 16;
hMarkBG = 16;

markBG = uibuttongroup('Units','Normalized','Title','Marking',...
    'BackgroundColor','magenta','Position',[83/wMax 0/hMax wMarkBG/wMax hMarkBG/hMax]);

saveSliceMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Save for Current Slice','FontSize',9,...
    'Position',[1/wMarkBG 13/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@saveSliceMarksButton_Callback);

loadSliceMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Load for Current Slice','FontSize',9,...
    'Position',[1/wMarkBG 9/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@loadSliceMarksButton_Callback);

saveAllMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Save for All Slices','FontSize',9,...
    'Position',[1/wMarkBG 5/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@saveAllMarksButton_Callback);

loadAllMarkButton = uicontrol('Style','pushbutton','Parent',markBG,'Units','normalized',...
    'String','Load for All Slices','FontSize',9,...
    'Position',[1/wMarkBG 1/hMarkBG 14/wMarkBG 3/hMarkBG],'Callback',@loadAllMarksButton_Callback);


set(hFig, 'Visible','on');

end