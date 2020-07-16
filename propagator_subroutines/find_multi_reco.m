function [t,nargs]=find_multi_reco(stinput, fs_file, argstr,reco_num)

    
curdir_reco=strcat(stinput,'/',num2str(fs_file),'/pdata/',num2str(reco_num),'/reco');
reco=fopen(curdir_reco,'r');

tline = fgets(reco);                                   % read the first line. fgets automatically moves to the next line after a line is read.
while tline ~= -1                                       %when fgets hits the end of file, -1 is returned
    tstr = findstr(tline,argstr);                           % is argstr present?
    if isempty(tstr)                                        % not yet
        tline = fgets(reco);
    else                                                    % argstr is found
        tstr = findstr(tline,cat(2,'##$',argstr,'=( '));    % is it an array? (note:all new parameters start with ##$ and all parameters which are arrays have '=( ' and then the number of elemets in the array to follow.
%         nargs_tmp = textscan(tline(tstr+19:tstr+21), '%s');
        nargs_tmp = textscan(tline, '%s');
        nargs=str2num(nargs_tmp{1}{2});
        if ~isempty(tstr)
            tline = fgets(reco);                       % start reading elements
            tline_final ='';                                % initialize final string
            
            % while next string does not start with ## or $$, i.e. new
            % variable
            while isempty(findstr(tline,'##')) & isempty(findstr(tline,'$$'))
                tline_final = cat(2,tline_final,tline);     % cat final string with a new one
                tline = fgets(reco);                   % read the next line in the file
            end                                             % end while
%             t=strread(tline_final,'%f',nargs,'delimiter',' '); % convert string into an array
            t=textscan(tline_final, '%s');
            fclose(reco);
            return
        else                                                % scalar value
            tstr = findstr(tline,cat(2,'##$',argstr,'='));
            if ~isempty(tstr)
                tstr = cat(2,'##$',argstr,'=');             % string that includes argstr
                strlength = size(tstr,2);                   % get the length its length
                tline(strlength-1:end);                     % read only post value
                t=str2num(tline(strlength+1:end));          % convert it to number
                fclose(reco);
               return
            else
                tline = fgets(reco);
            end
        end
    end
end
t=0;
fclose(reco);