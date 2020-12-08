function PPCmat = WMPhaseLocking
PPCmat = [];
for ifolder = [1:4 7]%1:7
    if ifolder > 1
        PathName = sprintf('../SameCircuits%dItemsSimultaneousV2',ifolder);
        cd(PathName)
    end
    PPC
    dir_strut = dir('*PPC_LFPNo1.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
    for i = 1%:num_files
        load(files{i},'PPCvalue');
        PPCmat = [PPCmat PPCvalue];
    end
end
end