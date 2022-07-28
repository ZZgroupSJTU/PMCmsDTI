function p = findPfile(filenames)
% Given a list of file names, this function finds the one which saves k-space data.
% Notice that from the scan, all Pfiles are named as 'Pxxxxx.7'. The one with the smallest index (earliest) is what we are looking for.

% filenames: list of all filenames
% p: index of the strating pfilename

index = 99999;
for n = 1 : size(filenames,1)
    name = filenames(n).name;
    a = strsplit(name, '.');
    if (length(a) == 2) && strcmp(a{2}, '7')
        if str2num(a{1}(2:6)) < index
            index = str2num(a{1}(2:6));
            p = a{1};
        end
    end
end
        


end

