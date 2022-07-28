clear all
% This demo loads

addpath(genpath(pwd))

% NX = 256;
% NY = 256; % If NX and NY are defined, then we are going to zero-filled the k-space (this is used for the deep learning data prep).


p.filename = ['/bmrNAS/people/yuxinh/highres_brain/Exam5908/30Apr19_Ex5908_Ser6']; 
% This variable is the path where the k-space data are saved. And it should be the "outputPath" if you are using our code to load k-space data.

p.ndir = 30; % This variable defines how many directions are acquired (not including the first non-diffusion-weighted image). Notice that
% each direction is saved into one individual file.

p.savepath = [p.filename, '/prep_sense_gcc']; % This defines the path to save the reconstructed results.
if (exist(p.savepath,'dir') ~= 7)
   mkdir(p.savepath);
end
    
p.lambda = 0.002; % Regularization parameter for the Locally low-rank reconstruction (here for initialiaztion we use SENSE).
p.gcc = true; % Whether to do geometric coil compression.
p.v = 8; % Number of virtual coil to be compresssed

p.reconmethod = 'SENSE'; % Reconstruction method for each slice and direction (options: SENSE and LLR).

p.b0 = []; % This is a list storing all the non-diffusion-weighted scans (not including the first b=0 scan). Usually it is empty.
% If using Qiyuan's tensor file for acquisition, set this to "1:15:p.ndir", since the first of every fifteen directions would be the non-diffusion-weighted scan.


for file_index = 1 : p.ndir
    p.dir = file_index;
    disp(['loading file from ',num2str(file_index),':',p.filename])
    SENSE_ini;
end

