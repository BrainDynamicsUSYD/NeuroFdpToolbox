function job_convertVideo(array_ID)

switch array_ID
    case 1
        Name = 'Fullamp0.5_4Hz0.avi' ;
    case 2
        Name = 'Fullamp0.5_4Hz5000.avi' ;
    case 3
        Name = 'Fullamp4_8Hz0.avi' ;
    case 4
        Name = 'Fullamp4_8Hz5000.avi' ;
    case 5
        Name = 'Fullamp30_100Hz0.avi' ;
    case 6
        Name = 'Fullamp30_100Hz5000.avi' ;
    case 7
        Name = 'Fullphase0.5_4Hz0.avi' ;
    case 8
        Name = 'Fullphase0.5_4Hz5000.avi' ;
    case 9
        Name = 'Fullphase4_8Hz0.avi' ;
    case 10
        Name = 'Fullphase4_8Hz5000.avi' ;
    case 11
        Name = 'Fullphase30_100Hz0.avi' ;
    case 12
        Name = 'Fullphase30_100Hz5000.avi' ;
end
cd ..
cd([pwd,'/Results/GammaBurst'])
v = VideoReader(Name) ;
%%
v2 = VideoWriter(['new',Name,'.avi'],'Motion JPEG AVI');
v2.Quality = 50 ;
open(v2)
while hasFrame(v)
    video = readFrame(v);
    writeVideo(v2,video)
end
close(v2)