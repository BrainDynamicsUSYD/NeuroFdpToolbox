% general post process
dir_strut = dir('Elec1Scales10Window*.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end
w = [];
f = [];
for i = 1:num_files
    fprintf('Loading RYG.mat file %s...\n', files{i});
    R = load(files{i});
    w = [w R.window];
    f = [f R.F];
end
A = polyfit(log(1e-3*w),log(f),1);
fitline = polyval(A,log(1e-3*w));
Alpha1 = A(1);
loglog(1e-3*w,f,'o');
xlabel('Window Size(s)')
ylabel('F(tau)')
hold on
loglog(1e-3*w,exp(fitline))
disp(Alpha1)