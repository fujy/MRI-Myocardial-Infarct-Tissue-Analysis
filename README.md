# MRI-Myocardial-Infarct-Tissue-Analysis

A MATLAB implementation for publication titled [Infarct Tissue Heterogeneity by Magnetic Resonance Imaging Identifies Enhanced Cardiac Arrhythmia Susceptibility in Patients With Left Ventricular Dysfunction](http://circ.ahajournals.org/content/115/15/2006.short)

##The application has the following features:

- Load DICOM images
- Display patient information stored in the DICOM image (Family Name, Given
Name, ID, ... etc)
- Change the images contrast
- Manually segment Endocardial and Epicardial borders, remote (normal) and Hyperenhanced zone
- Save marked data in .mat file to be used or edited later 
- Calculate and highlight Remote, Core and Grayzone
- Calculate area and volume for all manual and automatic segmented zones

## Dependencies

- MATLAB R2013a
- MATLAB Image Processing Toolbox


## License

Licensed under the [MIT license](http://www.opensource.org/licenses/MIT)