% This function calculate the T1 map with B1 correction from two flip
% angles. Tested on Siemens Prisma vendor sequence. 
% Usage:
%         T1_B1cor=T1cal_nii;

% The output will be T1_B1cor.img and T1_B1cor.hdr
% B1 correction is done by smoothing FA3 for 20 times.
%--- written by Humberto Monsivais Jan. 12, 2022 -------------
%--- based on Dr. Chien-Lin's code for T1 mapping from Wabash welder study --------
% Separated the individual functions for easier reading and changed script to 
% accept .nii files from Siemens scanner %

%===please had SPM installed in MATLAB====

% User selects .nii files for the two flip angles used in increasing order
% and also the B1 map file that has already been coregistered to the T1w
% images
function T1_B1cor=T1cal_nii_Siemens

addpath '/Users/humbertomonsivais/Documents/MATLAB/spm12';
%% Store paths to load the files
[n1, p1]=uigetfile('*.nii','please select FA3 Nifti file');
[n2, p2]=uigetfile('*.nii','please select FA17 Nifti file');
[n3, p3]=uigetfile('*.nii','please select scaled FA map Nifti file from tfl_b1map sequence');
%-----------------------------------------------------------------------
% determine output directory path
outpath = '/Users/humbertomonsivais/Documents/MATLAB/W107A'; % input output path

%% Load header information for the T1w signal files
Sn.FA1 = spm_vol([p1, n1]);
Sn.FA2 = spm_vol([p2, n2]);

if isfield(Sn.FA1,'descrip')
     d=strfind(Sn.FA1.descrip,'TR=');
     par.TR=str2double(Sn.FA1.descrip(d+3:d+4)); %6 for GE data
     display(par.TR);
else
     [n3, p3]=uigetfile('*.*','please select a dicom file of FA3');
     info=dicominfo([p3,n3]);
     par.TR=info.RepetitionTime;
end

%% Load header information for the B1 map files
Sn_B1map_FA = spm_vol([p3, n3]);

%% Read volumes
Vol1 = spm_read_vols(Sn_B1map_FA);

%% Generates normalized B1 map

B1map = Vol1; 
B1map_norm = (1/900)*B1map; % normalizes to 90 degrees

%% Save everything in OUTPUT dir
sname = spm_file(Sn_B1map_FA.fname,'basename');
VB1 = Sn_B1map_FA;
VB1 = rmfield(VB1, 'pinfo');  % Updates scale of each plane
VB1.fname = fullfile(outpath, [sname '_B1map_norm.nii']);
VB1.descrip = 'Normalised (p.u.) B1 bias map - TFL B1map protocol';
spm_write_vol(VB1,B1map_norm);

% =====================================================================================

%% Read volumes for t1w images
Vol_Sn.FA1 = spm_read_vols(Sn.FA1);
Vol_Sn.FA2 = spm_read_vols(Sn.FA2);

par.FA1=3;  % Enter first FA
par.FA2=17; % Enter second FA

%correct FA values with B1 values for each voxel
FA1=par.FA1*B1map_norm;%corrected first flip angle 
FA2=par.FA2*B1map_norm;%corrected second flip angle

%Find sind(FA); FA input is in degrees 
SIN1=sind(FA1);%sin(FA) for each flip angle
SIN2=sind(FA2);%sin(FA) for each flip angle

%Find tand(FA); FA input is in degrees
TAN1=tand(FA1);%tan(FA) for each flip angle
TAN2=tand(FA2);%tan(FA) for each flip angle

%% Old code to calculate the T1 map

%{
% Creates T1 map content and corrects for B1
[SUM_X, SUM_Y,n,d,mask]=T1_mapcont(Vol_Sn,par,B1map_norm);

%% Creates out T1_B1cor map
T1_B1cor = T1_map(SUM_X, SUM_Y,n,d,mask,par);
%}

%create matrix to store output data
sz = size(Vol1);
T1 = zeros(sz);
m = zeros(sz);
%c = zeros(sz);
NCores = 4;
parpool('local',NCores); %Uses parallel computing toolbox to speed up calculations
%most outter loop (k): loops over each slice
%outter for loop (i): loops over each row
%inner for loop (j): loops over each element in a row
for k = 1:sz(3)
    for i = 1:sz(1)
        for j = 1:sz(1)
            warning('off');
            %temp-store data points for analysis
            S = [Vol_Sn.FA1(i,j,k); Vol_Sn.FA2(i,j,k)];
            TAN = [TAN1(i,j,k); TAN2(i,j,k)];
            SIN = [SIN1(i,j,k); SIN2(i,j,k)];
        
            X=S./TAN;%the x-axis values: signal/tan(FA)
            Y=S./SIN;%the y-axis values: signal/sin(FA)

            x=[ones(length(X),1) X];%including a y-intercept to improve the fit
  
            %b is the slope (m) and y-int (c)
            b=x\Y;%a [2*1] matrix includes [c,m]
            m(i,j,k)=b(2,1);%the slope of the graph   
            %c(i,j,k)=b(1,1);%the constant in y=mx+c
            %y=x*b;%calculated y-axis values
            %R2(i,j)=1 - sum((Y - y).^2)/sum((Y - mean(Y)).^2);%R squared value for the linearity of the graph    
            T1(i,j,k) =-par.TR/log(abs(m(i,j,k)));%T1 in ms for this ROI
        end
    end
end 

%%
%plotting the results
T1_B1cor = T1;
T1image=T1_B1cor(:,:,80);%for example slice 32 of R1
T1image = imrotate(T1image,90);
imagesc(T1image)
colorbar
colormap(parula)%changing the color map
axis image off%turn off the auto scale
caxis([500 3.5e3])%setting up a min and max for axis
%caxis([0 max(max(R1))])%or letting it scale to the max singnal


R1_B1cor = 1/T1_B1cor*1e3; % Takes the T1 map data (in ms) and converts it to R1 in (1/s)

%% Saves T1_B1cor map
fname = spm_file(Sn.FA1.fname,'basename');
VT1 = Sn.FA1; % Stores hdr information from first t1w image
%VT1 = rmfield(VT1, 'pinfo'); % Updates scale of each plane
VT1.fname = fullfile(outpath, [fname '_T1map_B1cor.nii']);
VT1.descrip = 'T1 map (in ms) with B1 bias map correction';
spm_write_vol(VT1,T1_B1cor);

%% Saves R1_B1cor map
fname = spm_file(Sn.FA1.fname,'basename');
VR1 = Sn.FA1; % Stores hdr information from first t1w image
VR1 = rmfield(VR1, 'pinfo'); % Updates scale of each plane
VR1.fname = fullfile(outpath, [fname '_R1map_B1cor.nii']);
VR1.descrip = 'R1 map (in 1/s) with B1 bias map correction';
spm_write_vol(VR1,R1_B1cor);

end
