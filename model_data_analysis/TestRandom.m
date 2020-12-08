function TestRandom(varargin)
% test rand

% Loop number for PBS array job
loop_num = 0;
for i = 1:5
    
    loop_num = loop_num + 1;
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    date_now = datestr(now,'yyyymmddHHMM-SSFFF');
    scan_temp = textscan(date_now,'%s','Delimiter','-');
    rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
    alg='twister';
    rng(rand_seed,alg); %  Note that this effect is global!
    a = rand(5);
    disp(a)
end
end