function xml = xml2matlab(fname,varargin)
xml = [];
fid = fopen(fname,'r');
eof = 0;
if (nargin==2)
    prnt = varargin{1};
else
    prnt = 0;
end
L = '';
count = 0;
while ~eof
    while isempty(strtrim(L))
        L = fgetl(fid);
        if ~ischar(L); 
            eof=1; return; 
        end
    end
    [xml,eof,L] = getnexttag(xml,L,fid,count,[],eof,prnt);
end

fclose(fid);

function [xml,eof,L] = getnexttag(xml,L,fid,count,curnt,eof,prnt)
    
    while 1
        if ~ischar(L); 
            eof=1; return; 
        end 

        while (1)
            [ncwrong]= regexp(L,['</[^ /<>\f\r\n\t\v]*>']);
            [ncl,ncla,nclb,nclc,ncld,ncle,nclf]= regexp(L,['</' curnt '>']);
            % Any non-white space character except "/", "<", and ">"
            [nel,nela,nelb,nelc,neld,nele,nelf]= regexp(L,'<([^/][^ /<>\f\r\n\t\v]*)');
            %[nel,nela,nelb,nelc,neld,nele,nelf]= regexp(L,'<([^/][a-zA-Z0-9]*)') 
            ext = 0;
            nn = length(L);
            n = nn;
            if ~isempty(ncl) 
                if isempty(nel);        
                    ext = 1; n = ncl(1)-1; % if there is no new element
                elseif ncl(1) < nel(1); 
                    ext = 1; n = ncl(1)-1; % if closing the current element happens before opening a new element
                end;
            elseif ~isempty(ncwrong)
                % If no element is created before the non-curnt one is
                % closed, then there is an errori.
                if isempty(nel)
                    error(['Missing close for element "' curnt '"'])
                elseif nel > ncwrong
                    error(['Missing close for element "' curnt '"'])
                end
            end
            if (nn==n && ~isempty(nel))
                n = nel(1)-1;
            end
            if ~isempty(strtrim(L(1:n)))
                xml = valassign(xml,L(1:n),curnt,prnt);
            end
            if (nn == n)
                L = ' ';
                while isempty(strtrim(L))
                    L = fgetl(fid);
                    if ~ischar(L); 
                        eof=1; return; 
                    end
                end
            else
                break;
            end
        end

        if (ext); 
            if (prnt==1)
                disp(['close "' curnt '"'])
            end
            % Before close, if this has a value element try to make it numeric
            if (isfield(xml,'Value'))
                if strcmpi(strtrim(xml.Value),'true')
                elseif strcmpi(strtrim(xml.Value),'false')
                else
                    a = str2num(xml.Value);
                    if ~isempty(a)
                        xml.Value = a;
                    end
                end
            %else
            %    xml.Value = [];
            end

            if (ncla(1)<length(L))
                L = L(ncla(1)+1:end);
            else
                L = '';
            end
            return; 
        end

        if ~isempty(nel)
            elm = neld{1}{1};
            [nen,nena,nenb,nenc,nend,nene,nenf]= regexp(nelf{2},['[/]*>']);
            if strcmpi(elm,'?xml')
                if (nela(1)+nena(1) < length(L))
                    L = L(nela(1)+nena(1)+1:end);
                else
                    L = fgetl(fid);
                    if ~ischar(L); 
                        eof=1; return; 
                    end
                end
                return
            elseif (strfind(elm,'!--')==1)
                % DEALING WITH COMMENTS ----------------------------
                comments = elm(4:end);
                L = nelf{2};
                while (1)
                    [nen,nena,nenb,nenc,nend,nene,nenf]= regexp(L,['-->']);
                    if ~isempty(nenf)
                        comments = [comments nenf{1}];
                    end

                    if isempty(nen)
                        L = fgetl(fid);
                        if ~ischar(L); 
                            eof=1; return; 
                        end
                        comments = [comments char(10)];
                    else
                        if (nena(1)==length(L))
                            L = '';
                        else
                            L = L(nenb+1:length(L));
                        end
                        break;
                    end

                end
                [xml,count] = addelement(xml,'Comments__',count,curnt,prnt);
                xml.('Comments__'){end} = comments;
                % DEALING WITH COMMENTS ----------------------------
            else 
                [xml,count] = addelement(xml,elm,count,curnt,prnt);
                if ~isempty(nenf)
                    fin = 0;
                    % Attribute can be any non-white space character except <>
                    [nat,nata,natb,natc,natd,nate,natf]= ...
                        regexp(nenf{1},'([^ /<>\f\r\n\t\v]*)="([^"]*)"');
                    %[nat,nata,natb,natc,natd,nate,natf]= ...
                    %    regexp(nenf{1},'([a-zA-Z0-9]*)="([a-zA-Z0-9_\-]*)"');
                    for ii = 1:length(natd)
                        xml = attassign(xml,elm,natd{ii}{1},natd{ii}{2},curnt,prnt);
                    end
                    if (nela(1)+nena(1) < length(L))
                        L = L(nela(1)+nena(1)+1:end);
                    else
                        L = fgetl(fid);
                        if ~ischar(L); 
                            eof=1; return; 
                        end
                    end
                    if (strcmp(nenc{1},'/>')); % This should close the element
                        L = ['</' elm '>' L];
                    else
                        while isempty(strtrim(L))
                            L = fgetl(fid);
                            if ~ischar(L); 
                                eof=1; return; 
                            end
                        end
                    end
                    %disp(['IN************' elm])
                    %if iscell(xml.(elm))
                        [xml.(elm){end},eof,L] = getnexttag(xml.(elm){end},L,fid,count,elm,eof,prnt);
                    %else
                    %    [xml.(elm),eof,L] = getnexttag(xml.(elm),L,fid,count,elm,eof,prnt);
                    %end
                    %disp(['OUT************' elm])
                    if ~ischar(L); 
                        eof=1; return; 
                    end
                    while isempty(strtrim(L))
                        L = fgetl(fid);
                        if ~ischar(L); 
                            eof=1; return; 
                        end
                    end
                else
                    error('huh?')
                end
            end
        else
            L = fgetl(fid);
            while isempty(strtrim(L))
                L = fgetl(fid);
                if ~ischar(L); 
                    eof=1; return; 
                end
            end
        end
    end
    
    
    function xml = attassign(xml,elm,a,v,curnt,prnt)
        %if iscell(xml.(elm))
            xml.(elm){end}.Attributes.(a) = v;
        %else
        %    xml.(elm).Attributes.(a) = v;
        %end
        if (prnt==1)
            disp(['Add attibute [' a '="' v '"] to element "' elm '" of "' curnt '"'])
        end
    function xml = valassign(xml,L,curnt,prnt)
        if (prnt==1)
            disp(['Add value "' L '" to "' curnt '"'])
        end
        if ~isfield(xml,'Value'); 
            xml.Value = L; 
        else
            xml.Value = [xml.Value char(10) L];
        end
        
    function [xml,count] = addelement(xml,elm,count,curnt,prnt)
        p = 1;
        count = count+1;
        if isfield(xml,elm)
            p = length(xml.(elm))+1;
%             if p==2
%                 tmp = xml.(elm);
%                 xml = rmfield(xml,elm);
%                 xml.(elm){1} = tmp;
%                 clear tmp
%             end
            xml.(elm){p} = [];
        else
            %xml.(elm) = [];
            xml.(elm){1} = [];
        end
        xml.(elm){end}.Count__ = count;
        if (prnt==1)
            disp(['Add ' num2str(p) 'th element "' elm '" to "' curnt '"'])
        end
        