function animated_gif_creator(file_name, file_path, file_name2, file_path2, lps, delay, delay1, delay2)
switch lps
    case 'Forever'
        loops=65535;
    case 'None'
        loops=1;
    case 'Other'
        loops=inputdlg('Enter number of loops? (must be an integer between 1-65535)        .','Loops');
        loops=str2num(loops{1});
end

h = waitbar(0,['0% done'],'name','Progress') ;
for i=1:length(file_name)
    if strcmpi('gif',file_name{i}{1}(end-2:end))
        [M  c_map]=imread([file_path,file_name{i}{1}]);
    else
        a=imread([file_path,file_name{i}{1}]);
        [M  c_map]= rgb2ind(a,256);
    end
    if i==1
        imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay1)
    elseif i==length(file_name)
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay2)
    else
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
    end
    waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;
end
close(h);
end