clear all;
% This demo loads the kspace data from the Pfiles with modified Orchestra (Matlab, SDK 1.7-1).

addpath(genpath('/bmrNAS/people/yuxinh/bart/orchestra-sdk-1.7-1.matlab')) % Add installed Orchestra into the path. You probably need to change this.
addpath(genpath(pwd))
pdirname = ['/bmrNAS/scandata/2019/','02May19_Ex5916_Ser']; % where pfiles (and vrgf and ref files) are save
dirlist = [3]; % sometimes there might be multiple folders/scans. This (with the following for-loop) saves some code.

dirname = ['/bmrNAS/people/yuxinh/test_for_github/spa_llr/','02May19_Ex5916_Ser']; % where to save files 


for n = 1 : length(dirlist)
    dirname1 = [pdirname, int2str(dirlist(n))];
    filenames = dir(dirname1);
    p = findPfile(filenames);
    EpiDiffusionRecon_yuxin([dirname1,'/',p,'.7'], [dirname, int2str(dirlist(n))]);
end
