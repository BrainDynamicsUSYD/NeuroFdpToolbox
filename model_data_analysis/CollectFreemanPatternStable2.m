function [ApiceP,ApiceN] = CollectFreemanPatternStable2
% example on FreemanPatternStable Aone/two
% Due to heavy calculation, this is a simplied version to get data for
% plotting.
% Ref: Freeman, Walter J., and John M. Barrie. "Analysis of spatial patterns
%     of phase in neocortical gamma EEGs in rabbit." Journal of neurophysiology 84.3 (2000): 1266-1278.
tic;
% simplified = 100 300 
% tictoc     =207s 440s
% plotq      = 12  12
simplified = 100;
dir_strut = dir('0001-*NoFilterHz.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
step_tot = 1e5; % ms
seg_size = 1280; % 128 ms
gap = 20; % 2 ms
seg_num = floor((step_tot-seg_size)/gap) + 1;
Awi = cell(1,seg_num);
Awj = cell(1,seg_num);
AONE = []; % look like [fitr.x_c fitr.y_c Rwf -sign(fitr.a) i]
ATWO = [];
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    AONE = [AONE; R.Aone];
    ATWO = [ATWO; R.Aone];
end
clear dir_strut files R
for i = 1:seg_num
    Ind1 = find(AONE(:,5)==i);
    awi = AONE(Ind1,:);
    [~,Ind2] = sort(awi(:,3));
    Awi{i} = awi(Ind2,1:4);
    clear Ind1 Ind2
    Ind1 = find(ATWO(:,5)==i);
    awj = ATWO(Ind1,:);
    [~,Ind2] = sort(awj(:,3));
    Awj{i} = awj(Ind2,1:4);
    clear Ind1 Ind2
end
clear AONE ATWO awi awj
AwiP = cell(1,seg_num); % look like [fitr.x_c fitr.y_c Rwf -sign(fitr.a)]
AwjP = cell(1,seg_num);
AwiN = cell(1,seg_num);
AwjN = cell(1,seg_num);
HeavyP = [];
HeavyN = [];
for i = 1:seg_num
    Ind = find(Awi{i}(:,4)==1);
    AwiP{i} = Awi{i}(Ind,1:2);
    clear Ind
    Ind = find(Awi{i}(:,4)==-1);
    AwiN{i} = Awi{i}(Ind,1:2);
    clear Ind
    Ind = find(Awj{i}(:,4)==1);
    AwjP{i} = Awj{i}(Ind,1:2);
    clear Ind
    Ind = find(Awj{i}(:,4)==-1);
    AwjN{i} = Awj{i}(Ind,1:2);
    clear Ind
end
clear Awi Awj
ApiceP = cell(1,seg_num-26);
for i = 1:seg_num-26
    Mati = AwiP{i};
    q = 1;
    Matj = AwjP{i+q};
    Matk = AwjP{i+q+1};
    VVd2 = DistanceAmongMatrix(Matj,Matk);
    Ind2 = find(VVd2(:,5)*60e-3 < 0.5/7*0.6); % /7*0.6 % mm
    VVd2 = VVd2(Ind2,:);
    Candidate = [];
    for j = 1:length(Ind2)
        VVd1 = DistanceAmongMatrix(Mati,VVd2(j,1:2));
        for k = 1:size(Mati,1)
            if VVd1(k,5)*60e-3 < 1.2/7*0.6
                candidate = [Mati(k,:) VVd2(j,1:4)];
                Candidate = [Candidate;candidate];
            end            
        end
    end
    while q > 0
        mark = 0;
        Matk = AwjP{i+q+1};
        New = [];
        for k = 1:size(Matk,1)
            for j = 1:size(Candidate,1)
                VVd1 = DistanceAmongMatrix(Candidate(j,1:2),Candidate(j,end-1:end));
                VVd2 = DistanceAmongMatrix(Candidate(j,end-1:end),Matk(k,:));
                if VVd1(5)*60e-3 < 1.2/7*0.6 && VVd2(5)*60e-3 < 0.5/7*0.6
                    mark = 1;
                    new = [Candidate(j,:) Matk(k,:)];
                    New = [New;new];
                end
            end            
        end
        if size(New,1) > simplified
            r = randperm(size(New,1),simplified);
            New = New(r,:);
        end
        if mark == 1
            q = q + 1;
            Candidate = New;
        elseif q >= 25
            ApiceP{i} = Candidate;
            break
        elseif ~isempty(Mati)
            ApiceP{i} = Mati(1,:);
            break
        else
            break
        end
    end
end
clear AwiP AwjP
ApiceN = cell(1,seg_num-26);
for i = 1:seg_num-26
    Mati = AwiN{i};
    q = 1;
    Matj = AwjN{i+q};
    Matk = AwjN{i+q+1};
    VVd2 = DistanceAmongMatrix(Matj,Matk);
    Ind2 = find(VVd2(:,5)*60e-3 < 0.5/7*0.6); % /7*0.6 % mm
    VVd2 = VVd2(Ind2,:);
    Candidate = [];
    for j = 1:length(Ind2)
        VVd1 = DistanceAmongMatrix(Mati,VVd2(j,1:2));
        for k = 1:size(Mati,1)
            if VVd1(k,5)*60e-3 < 1.2/7*0.6
                candidate = [Mati(k,:) VVd2(j,1:4)];
                Candidate = [Candidate;candidate];
            end            
        end
    end
    while q > 1
        mark = 0;
        Matk = AwjN{i+q+1};
        New = [];
        for k = 1:size(Matk,1)
            for j = 1:size(Candidate,1)
                VVd1 = DistanceAmongMatrix(Candidate(j,1:2),Candidate(j,end-1:end));
                VVd2 = DistanceAmongMatrix(Candidate(j,end-1:end),Matk(k,:));
                if VVd1(5)*60e-3 < 1.2/7*0.6 && VVd2(5)*60e-3 < 0.5/7*0.6
                    mark = 1;
                    new = [Candidate(j,:) Matk(k,:)];
                    New = [New;new];
                end
            end            
        end
        if size(New,1) > simplified
            r = randperm(size(New,1),simplified);
            New = New(r,:);
        end
        if mark == 1
            q = q + 1;
            Candidate = New;
        elseif q >= 25
            ApiceN{i} = Candidate;
            break
        elseif ~isempty(Mati)
            ApiceN{i} = Mati(1,:);
            break
        else
            break
        end
    end
end
toc;
end