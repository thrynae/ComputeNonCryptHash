function [timing,description]=aaa___ComputeNonCryptHash___performance_test
%This tester requires internet access to create files that will be loaded to persistent variables.
%
%The list of words is based on this list:
%http://web.archive.org/web/20061231134037/http://www.sitopreferito.it/html/all_english_words.html
%
%The images are from the Stanford Dog Dataset containing 20580 images (with 89 duplicates)
%         http://vision.stanford.edu/aditya86/ImageNetDogs/

n_tests=5;timing=zeros(4,n_tests);
description=cell(1,n_tests);
if nargout==0,k=0;K=numel(timing);h_wait=waitbar(0,'progress');end
for n=1:size(timing,1)
    [timing(n,1),description{1}]=test_1;
    if nargout==0,k=k+1;waitbar(k/K,h_wait);end
    [timing(n,2),description{2}]=test_2;
    if nargout==0,k=k+1;waitbar(k/K,h_wait);end
    [timing(n,3),description{3}]=test_3;
    if nargout==0,k=k+1;waitbar(k/K,h_wait);end
    [timing(n,4),description{4}]=test_4;
    if nargout==0,k=k+1;waitbar(k/K,h_wait);end
    [timing(n,5),description{5}]=test_5;
    if nargout==0,k=k+1;waitbar(k/K,h_wait);end
end
get_persistent([],'clear');
timing=median(timing(2:end,:),1);
if nargout==0
    delete(h_wait)
    fprintf('Timing results:\n')
    for k=1:numel(timing),fprintf('%.3f seconds (test: %s)\n',timing(k),description{k});end
    clear timing description
end
end

function [timing,description]=test_1
description='%d English words';
N=10000;
Words=get_persistent('Words','load',N);
then=now;for n=1:N,ComputeNonCryptHash(Words{n},128);end,timing=(now-then)*(24*60*60);
description=sprintf(description,N);
end
function [timing,description]=test_2
description='1 to 1e%d (in char format)';
E=4;
persistent Numbers N
if isempty(Numbers) || numel(Numbers)~=N
    N=10^E;
    Numbers=cell(1,N);
    for n=1:numel(Numbers),Numbers{n}=sprintf('%d',n);end
end
then=now;for n=1:N,ComputeNonCryptHash(Numbers{n},128);end,timing=(now-then)*(24*60*60);
description=sprintf(description,E);
end
function [timing,description]=test_3
description='1 to 1e%d (in double format)';
E=4;
persistent Numbers N
if isempty(Numbers) || numel(Numbers)~=N
    N=10^E;
    Numbers=num2cell(1:N);
end
then=now;for n=1:N,ComputeNonCryptHash(Numbers{n},128);end,timing=(now-then)*(24*60*60);
description=sprintf(description,E);
end
function [timing,description]=test_4
description='1 to 1e%d (in uint16 format)';
E=4;
persistent Numbers N
if isempty(Numbers) || numel(Numbers)~=N
    N=10^E;
    Numbers=num2cell(uint16(1:N));
end
then=now;for n=1:N,ComputeNonCryptHash(Numbers{n},128);end,timing=(now-then)*(24*60*60);
description=sprintf(description,E);
end
function [timing,description]=test_5
description='%d images';
N=25;
IM=get_persistent('IM','load',N);
then=now;for n=1:N,ComputeNonCryptHash(IM{n},128);end,timing=(now-then)*(24*60*60);
description=sprintf(description,N);
end
function out=get_persistent(var,operation,varargin)
%allow central loading and saving to prevent a memory hog
persistent Words IM
switch operation
    case 'load'
        if strcmp(var,'Words')
            if nargin>2,N=varargin{1};else,N=10000;end
            if isempty(Words),Words=get_Words(N);end
            out=Words;
        elseif strcmp(var,'IM')
            if nargin>2,N=varargin{1};else,N=25;end
            if isempty(IM),IM=download_example_images(N);end
            out=IM;
        else
            error('name not recognized')
        end
    case 'clear'
        Words=[];IM=[];
end
end
function Words=get_Words(N)
%store to a file in the temp folder
fn=fullfile(tempdir,'MATLAB','FileExchange','ComputeNonCryptHash','EnglishWordList.txt');
if ~exist(fn,'file') %load from internet
    txt=readfile(['http://web.archive.org/web/20201006125927id_/https://cdn-140.anonfiles.com/',...
        'B89ePdcdp9/1fa52173-1601989741/EnglishWordList.txt']);
    if ~exist(fileparts(fn),'dir'),mkdir(fileparts(fn));end
    txt=sprintf('%s\n',txt{:});txt(end)='';
    fid=fopen(fn,'w');fprintf(fid,'%s',txt);fclose(fid);
end
txt=readfile(fn);
for n=1:numel(txt),if strcmp(txt{n}(1),'%'),txt{n}='';end,end
txt(cellfun('isempty',txt))=[];
Words=txt;
if numel(Words)~=N
	Words=Words(unique(round(linspace(1,numel(Words),N))));%select a semi-random sample
end
end
function IM=download_example_images(N)
%download 25 images from the Stanford Dogs Dataset
%the links allow up to 5013 images (5043 if duplicates are not removed)
IM=cell(1,N);

p1='http://web.archive.org/web';
p2='http://vision.stanford.edu/aditya86/ImageNetDogs';
k=0;
%Affenpinscher up to Chesapeake Bay retriever
k=k+1;url{k}=[p1 '/20190429035422fw_/' p2 '/n02110627.html'];
k=k+1;url{k}=[p1 '/20190907200724fw_/' p2 '/n02088094.html'];
k=k+1;url{k}=[p1 '/20181211003822fw_/' p2 '/n02116738.html'];
k=k+1;url{k}=[p1 '/20190429035952fw_/' p2 '/n02096051.html'];
k=k+1;url{k}=[p1 '/20190429035720fw_/' p2 '/n02093428.html'];
k=k+1;url{k}=[p1 '/20190429041351fw_/' p2 '/n02107908.html'];
k=k+1;url{k}=[p1 '/20190429040002fw_/' p2 '/n02096294.html'];
k=k+1;url{k}=[p1 '/20190430212644fw_/' p2 '/n02110806.html'];
k=k+1;url{k}=[p1 '/20190430212609fw_/' p2 '/n02088238.html'];
k=k+1;url{k}=[p1 '/20190429035644fw_/' p2 '/n02088364.html'];
k=k+1;url{k}=[p1 '/20190430215444fw_/' p2 '/n02093647.html'];
k=k+1;url{k}=[p1 '/20190413043111fw_/' p2 '/n02107683.html'];
k=k+1;url{k}=[p1 '/20190429035654fw_/' p2 '/n02089078.html'];
k=k+1;url{k}=[p1 '/20190430212559fw_/' p2 '/n02086646.html'];
k=k+1;url{k}=[p1 '/20190429035649fw_/' p2 '/n02088466.html'];
k=k+1;url{k}=[p1 '/20181211003752fw_/' p2 '/n02088632.html'];
k=k+1;url{k}=[p1 '/20190430215535fw_/' p2 '/n02106166.html'];
k=k+1;url{k}=[p1 '/20190430215449fw_/' p2 '/n02093754.html'];
k=k+1;url{k}=[p1 '/20190429035659fw_/' p2 '/n02090622.html'];
k=k+1;url{k}=[p1 '/20181211003802fw_/' p2 '/n02096585.html'];
k=k+1;url{k}=[p1 '/20190429035830fw_/' p2 '/n02106382.html'];
k=k+1;url{k}=[p1 '/20190429035406fw_/' p2 '/n02108089.html'];
k=k+1;url{k}=[p1 '/20190429035850fw_/' p2 '/n02112706.html'];
k=k+1;url{k}=[p1 '/20190429035825fw_/' p2 '/n02105251.html'];
k=k+1;url{k}=[p1 '/20190429035810fw_/' p2 '/n02101388.html'];
k=k+1;url{k}=[p1 '/20190430215550fw_/' p2 '/n02108422.html'];
k=k+1;url{k}=[p1 '/20190429035957fw_/' p2 '/n02096177.html'];
k=k+1;url{k}=[p1 '/20190429041406fw_/' p2 '/n02113186.html'];
k=k+1;url{k}=[p1 '/20190430212624fw_/' p2 '/n02099849.html'];

%ignore/remove duplicate images
files_to_be_ignored={...
    'n02086646_172.jpg';'n02086646_280.jpg';'n02086646_422.jpg';'n02086910_103.jpg';
    'n02088466_8078.jpg';'n02088632_4584.jpg';'n02088632_600.jpg';'n02088632_982.jpg';
    'n02089867_3177.jpg';'n02090379_3300.jpg';'n02090379_855.jpg';'n02091831_6323.jpg';
    'n02092002_4745.jpg';'n02092002_5557.jpg';'n02092002_6503.jpg';'n02093428_1378.jpg';
    'n02093428_5635.jpg';'n02093428_8454.jpg';'n02093647_2687.jpg';'n02093647_3277.jpg';
    'n02093647_518.jpg';'n02093859_324.jpg';'n02093991_4670.jpg';'n02094258_1408.jpg';
    'n02094258_880.jpg';'n02095314_3227.jpg';'n02095314_80.jpg';'n02095570_3101.jpg';
    'n02095570_5404.jpg';'n02095570_6443.jpg';'n02095889_1582.jpg';'n02095889_2471.jpg';
    'n02095889_6458.jpg';'n02097130_5576.jpg';'n02097209_583.jpg';'n02097209_585.jpg';
    'n02097474_7300.jpg';'n02097658_8018.jpg';'n02098105_2456.jpg';'n02098105_2945.jpg';
    'n02098286_830.jpg';'n02099429_3513.jpg';'n02101556_5903.jpg';'n02102040_639.jpg';
    'n02102177_2388.jpg';'n02102480_3436.jpg';'n02102973_1714.jpg';'n02102973_4469.jpg';
    'n02104365_5670.jpg';'n02104365_7953.jpg';'n02105056_2194.jpg';'n02105056_785.jpg';
    'n02105162_6449.jpg';'n02105251_7647.jpg';'n02106166_1031.jpg';'n02106166_346.jpg';
    'n02106166_4107.jpg';'n02106166_6084.jpg';'n02106166_6512.jpg';'n02106166_6710.jpg';
    'n02107908_3531.jpg';'n02107908_47.jpg';'n02108000_331.jpg';'n02109525_7982.jpg';
    'n02110185_1614.jpg';'n02110185_4115.jpg';'n02112706_2070.jpg';'n02112706_2074.jpg';
    'n02112706_2161.jpg';'n02112706_2191.jpg';'n02112706_420.jpg';'n02112706_539.jpg';
    'n02112706_866.jpg';'n02113186_11400.jpg';'n02113712_1252.jpg';'n02113712_1546.jpg';
    'n02113712_1558.jpg';'n02113799_254.jpg';'n02113978_2572.jpg';'n02113978_3474.jpg';
    'n02113978_3832.jpg';'n02113978_530.jpg';'n02113978_737.jpg';'n02115913_2029.jpg';
    'n02115913_3270.jpg';'n02115913_4032.jpg';'n02115913_564.jpg';'n02115913_739.jpg';
    'n02115913_750.jpg'};

n=0;then=now;disp_ETA_flag=false;
parent_folder=fullfile(tempdir,'MATLAB','FileExchange','ComputeNonCryptHash','images');
if ~exist(parent_folder,'dir'),mkdir(parent_folder);end
for breed=1:numel(url)
    if sum(~cellfun('isempty',IM))==N,continue,end
    local_file=fullfile(parent_folder,sprintf('dog_breed_%02d.txt',breed));
    if ~exist(local_file,'file')
        str=readfile(url{breed});
        fid=fopen(local_file,'w');fprintf(fid,'%s\n',str{:});fclose(fid);
    end
    str=readfile(local_file);
    
    ind=find(ismember(str,{'<div class="highslide-gallery">'}));
    ind=ind+2;
    while strcmp(str{ind}(1:4),'<a c') && sum(~cellfun('isempty',IM))<N
        elem=str{ind};
        elem=strrep(elem,'<a class="highslide" href="','/');
        elem=strrep(elem,'" onclick="return hs.expand(this)">','');
        ind=ind+2;
        imurl=[p2 elem];
        [discard,fn,ext]=fileparts(elem);fn=[fn,ext]; %#ok<AGROW,ASGLU>
        if ismember(fn,files_to_be_ignored),continue,end
        fn=fullfile(parent_folder,fn);
        if ~exist(fn,'file')
            WBM(fn,imurl,'flag','id','tries',[3 1 2]);
            disp_ETA_flag=true;
        end
        n=n+1;
        try IM{n}=imread(fn);catch,end
        if disp_ETA_flag,ETA_disp(N,n,then);end
    end
end
IM(cellfun('isempty',IM))=[];
end

% The functions below were minified to make them more compact.
%
% Included functions (internal dependencies not listed):
%  - ETA_disp (version 1.0.1)
%  - readfile (version 3.0.0)
%  - WBM (version 2.0)

function v001=ETA_disp(v002,v003,v004),persistent v005,if isempty(v005)
v005.addtodate=f11('<',7,'Octave','<',4);end,v003=max(eps,v003);v006=(now-v004)*(24*60*60);v008=floor(v006/60);v007=round(v006-60*v008);
v010=v003/v002;v009=min(10^10,round((v006-v010*v006)/v010));if v005.addtodate,v011=datestr(now+(v009/(24*60*60)),'HH:MM:SS');
else,v011=datestr(addtodate(now,v009,'second'),'HH:MM:SS');end,v012=v006/v003;v014=floor(v012/60);v013=v012-60*v014;v001=cell(1,3);
v001{1}=sprintf('%05.1f%% done after a total time of %02d:%02d.',100*v003/v002,v008,v007);v001{2}=sprintf('estimated time of completion: %s',v011);
v001{3}=sprintf('(%02d:%05.2f per iteration, %d iterations left)',v014,v013,v002-v003);if nargout==0,clc,fprintf('%s\n',v001{:});clear v001;end,end
function v015=readfile(v016,varargin),if nargin<1,error('HJW:readfile:nargin','Incorrect number of input argument.')
end,if ~(nargout==0 || nargout==1),error('HJW:WBM:nargout','Incorrect number of output argument.')
end,[v017,v018,v019]=f20(v016,varargin{:});if ~v017,rethrow(v019),else,[v005,v020,v021]=deal(v018.legacy,v018.UseURLread,v018.err_on_ANSI);
v022=struct('con',v018.print_2_con,'fid',v018.print_2_fid,'obj',v018.print_2_obj);end,if isa(v016,'string'),v016=char(v016);end
if numel(v016)>=8 && ( strcmp(v016(1:7),'http://') || strcmp(v016(1:8),'https://') ),if ~v005.allows_https && strcmp(v016(1:8),'https://')
f30(v022,'HJW:readfile:httpsNotSupported',['This implementation of urlread probably doesn''t allow https requests.',char(10),'The next lines of code will probably result in an error.'])
end,v023=f19(v016,v020,v022);if isa(v023,'cell'),v015=v023;
else,v024=true;v023=f05(v023,v024);try [v025,v026,v027]=f28(v023);catch,v019=lasterror;if strcmp(v019.identifier,'HJW:UTF8_to_unicode:notUTF8')
v026=false;else,f06(v022,v019),end,end,if v026,v023=f25(v027);end,v015=f02(v023);end,else,v015=f18(v016,v022,v021);end,end
function v028=WBM(v016,v029,varargin),f35;if nargin<2,error('HJW:WBM:nargin','Incorrect number of input argument.'),end,if ~(nargout==0 || nargout==1)
error('HJW:WBM:nargout','Incorrect number of output argument.'),end,[v017,v018,v019]=f33(v016,v029,varargin{:});if ~v017,rethrow(v019),else
[v030,v031,v032,v033,v034,v035,v036,v037]=deal(v018.date_part,v018.tries,v018.response,v018.ignore,v018.verbose,v018.UseURLwrite,v018.flag,v018.err429);
v038=v031(2)>0;v022=struct('con',v018.print_2_con,'fid',v018.print_2_fid,'obj',v018.print_2_obj);v018.print_2=v022;
end,if ~f23(fileparts(v016)),f06(v022,'HJW:WBM:NoWriteFolder','The target folder doesn''t exist or Matlab doesn''t have write permission for it.')
end,v039=cellfun('length',v032(:,2));[v039,v040]=sort(v039);v040=v040(end:-1:1);v032=v032(v040,:);
v039=v039(end:-1:1);v041=1;v017=false;v042=[];v043=[];v044=0;while ~v017 &&sum(v031(1:2))>0 && v031(3)>=0,if v031(v041)<=0,v041=3-v041;
end,v045=v041;try if v045==1,v046=false;v031(v045)=v031(v045)-1;if v035,v028=urlwrite( ['http://web.archive.org/web/' v030 v036 '_/' v029],v016);
v028=f04(v016,v028);else,v028=websave(v016,['https://web.archive.org/web/' v030 v036 '_/' v029],weboptions('Timeout',10));
end,elseif v045==2,v046=true;v031(v045)=v031(v045)-1;if v035,v028=urlwrite( ['http://web.archive.org/save/' v029],v016);
v028=f04(v016,v028);else,v028=websave(v016,['https://web.archive.org/save/' v029],weboptions('Timeout',10));end,end,v017=true;v044=0;
if v038 && ~f03(v028,v018,v046),v017=false;v041=2;end,catch,v019=lasterror;v017=false;if ~f12,while ~f12,v047=datestr(now,'HH:MM:SS');
if v034>=1,f30(v022,'Internet connection down, retrying in %d seconds (@%s)',2^v044,v047),end,pause(2^v044),v044=min(1+v044,6);
end,continue,end,v044=0;v048=v019.identifier;v048=strrep(v048,':urlwrite:',':webservices:');if strcmp(v048,'MATLAB:webservices:Timeout')
v049=4080;v031(3)=v031(3)-1;else,v050=strrep(v048,'MATLAB:webservices:HTTP','');v050=strrep(v050,'StatusCodeError','');v050=str2double(v050);
if isnan(v050),v049=-1;if v034>=2,f30(v022,v019.identifier,v019.message),end,else,switch v019.message,case 'urlwrite: Couldn''t resolve host name'
v049=404;case ['urlwrite: Peer certificate cannot be ','authenticated with given CA certificates'],v049=403;
otherwise,v049=v050;end,end,end,if v034>=3,v051=sprintf('Error %d tries(%d,%d,%d) (download of %s)',double(v049),v031(1),v031(2),v031(3),v016);
if v022.fid.boolean,fprintf(v022.fid.fid,'%s\n',v051);
end,if v022.obj.boolean,set(v022.obj.obj,'String',v051),end,drawnow,end,if ~any(v049==v033),v042(end+1)=v049;v043(end+1)=v045;
for v052=1:size(v032,1),if length(v042) < v039(v052),continue,end,v053=v042( (end-v039(v052)+1):end);v054=v043( (end-v039(v052)+1):end);
v055=v032{v052,1}(2:2:end);v055=strrep(v055,'x',num2str(v045));v056=strcmp(v055,sprintf('%d',v054));if isequal( v032{v052,2},v053) && v056
switch v032{v052,3},case 'load',v041=1;case 'save',v041=2;case 'exit',v031=[0 0 -1];case 'pause_retry',if ~v037.CountsAsTry
v031(v041)=v031(v041)+1;end,if v034>=v037.PrintAtVerbosityLevel,v057=10;v058='Waiting a while until the server won''t block us anymore';
if v022.fid.boolean,fprintf(v022.fid.fid,v058);drawnow,end,for v003=1:v057,pause(v037.TimeToWait/v057),if v022.fid.boolean
fprintf(v022.fid.fid,'.');drawnow,end,if v022.obj.boolean,v058=[v058 '.'];set(v022.obj.obj,v058);drawnow,end,end,if v022.fid.boolean
fprintf(v022.fid.fid,'\nContinuing\n');drawnow,end,if v022.obj.boolean,set(v022.obj.obj,'Continuing');drawnow,end,else,pause(v037.TimeToWait),end
end,break,end,end,end,end,end,if ~v017 || ( ~v038 && ~f03(v028,v018,v046) ),if exist(v016,'file'),delete(v016);end,v059=[v016 '.html'];v079=dir(v059);
if numel(v079)==1 && v079.bytes==0 && abs(datenum(v079.date)-now)<=(1/24),delete(v059);end,v028=[];end,v059=[v016 '.html'];v079=dir(v059);
if numel(v079)==1 && v079.bytes==0 && abs(datenum(v079.date)-now)<=(1/24),try delete(v059);catch,end,end,if nargout==0,clear('outfilename');end,end
function v001=f01(v060,v061),try v001=v060+v061;catch,try v001=bsxfun(@plus,v060,v061);
catch,v063=size(v060); v062=size(v061);v060=repmat(v060,max(1,v062./v063)); v061=repmat(v061,max(1,v063./v062));v001=v060+v061;end,end,end
function v064=f02(v023),v065=[0 find(v023==10) numel(v023)+1];
v064=cell(numel(v065)-1,1);for v003=1:numel(v064),v066=(v065(v003 )+1);v067=(v065(v003+1)-1);v064{v003}=v023(v066:v067);end,end
function v068=f03(v028,v018,v046),if ~exist(v028,'file'),v068=false;return
end,[v069,v070,v022]=deal(v018.m_date_r,v018.date_bounds,v018.print_2);if ~v046,v071=f31;end,v072='<input type="hidden" name="date" value="';
v015=readfile(v028);v073=0;while v073<=numel(v015) && (v073==0 || ~strcmp(f21(v015{v073}),'<td class="u" colspan="2">'))
v073=v073+1;end,if numel(v015)>=(v073+1),v074=v015{v073+1};
v075=strfind(v074,v072);v075=v075+length(v072)-1;v076=str2double(v074(v075+(1:14)));v068= v070.double(1)<=v076 && v076<=v070.double(2);
return,end,v015=v015(:)';v015=cell2mat(v015);v075=strfind(v015,'/web/');if numel(v075)==0,if v069==0,v068=true;
return,elseif v069==1,f30(v022,'HJW:WBM:MissingDateWarning','No date found in file, unable to check date, assuming it is correct.'),v068=true;return
elseif v069==2,f06(v022,'HJW:WBM:MissingDateError',['Could not find date. This can mean there is an ','error in the save. Try saving manually.'])
end,end,v077=zeros(size(v075));v015=[v015 'abcdefghijklmnopqrstuvwxyz'];if exist('isstrprop','builtin')
for v003=1:length(v075),for v078=1:14,if ~isstrprop(v015(v075(v003)+4+v078),'digit'),break,end,end,v077(v003)=str2double(v015(v075(v003)+4+(1:v078)));
end,else,for v003=1:length(v075),for v078=1:14,if ~any(double(v015(v075(v003)+4+v078))==(48:57)),break
end,end,v077(v003)=str2double(v015(v075(v003)+4+(1:v078)));end,end,[v079,v080,v064]=unique(v077);try [v080,v081]=max(histc(v064,1:max(v064)));catch
[v080,v081]=max(accumarray(v064,1));end,v076=v079(v081);v068= v070.double(1)<=v076 && v076<=v070.double(2);if ~v046,if v076<1e4,if v069==0,v068=true;
return,elseif v069==1,f30(v022,'HJW:WBM:MissingDateWarning','No date found in file, unable to check date, assuming it is correct.'),v068=true;return
elseif v069==2,f06(v022,'HJW:WBM:MissingDateError',['Could not find date. This can mean there is an ','error in the save. Try saving manually.'])
end,end,v074=sprintf('%014d',v076);v074={v074(1:4),v074(5:6),v074(7:8),v074(9:10),v074(11:12),v074(13:14)};
v074=str2double(v074);timediff=(v071-datenum(v074))*24*60*60;if timediff<10,v068=false;elseif timediff<60
warning('HJW:WBM:LivePageStored',['The live page might have been saved instead of a capture.',char(10),'Check on the WBM if a capture exists.'])
end,end,end
function v028=f04(v016,v028),v083=[pwd filesep v016];if ~strcmp(v028,v083),if ~exist(v028,'file') &&exist(v083,'file'),v028=v083;end,end,end
function v023=f05(v023,v084)
persistent v085 v086,if isempty(v085),v087=[338 140;339 156;352 138;353 154;376 159;381 142;382 158;402 131;710 136;732 152;
8211 150;8212 151;8216 145;8217 146;8218 130;8220 147;8221 148;8222 132;8224 134;8225 135;
8226 149;8230 133;8240 137;8249 139;8250 155;8364 128;8482 153];v085=v087(:,2);v086=v087(:,1);end,if nargin>1 && v084
v089=v086;v088=v085;else,v089=v085;v088=v086;end,v023=uint32(v023);for v078=1:numel(v089),v023=f17(v023,v089(v078),v088(v078));end,end
function f06(v090,varargin)
if isempty(v090),v090=struct;end,if ~isfield(v090,'fid'),v090.fid.boolean=false;end,if ~isfield(v090,'obj'),v090.obj.boolean=false;end,if nargin==2
if isa(varargin{1},'struct') || isa(varargin{1},'MException'),v019=varargin{1};try v091=v019.stack;v092=f07(0,v091);catch,[v092,v091]=f07(2);
end,v093=v019.identifier;v094=v019.message;v095='Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback(';
if isa(v019,'struct') && numel(v094)>numel(v095) && strcmp(v095,v094(1:numel(v095)))
v094(1:find(v094==10,1))='';end,else,[v092,v091]=f07(2);[v093,v094]=deal('',varargin{1});end,else,[v092,v091]=f07(2);
if ~isempty(strfind(varargin{1},'%')),v093='';v096=varargin(2:end);v094=sprintf(varargin{1},v096{:});else,v093=varargin{1};v094=varargin{2};
if nargin>3,v096=varargin(3:end);v094=sprintf(v094,v096{:});end,end,end,v019=struct('identifier',v093,'message',v094,'stack',v091);if v090.obj.boolean
v097=v094;while v097(end)==10,v097(end)='';end,if any(v097==10),v097=regexp_outkeys(['Error: ' v097],char(10),'split');else,v097=['Error: ' v097];
end,set(v090.obj.obj,'String',v097),end,if v090.fid.boolean,fprintf(v090.fid.fid,'Error: %s\n%s',v094,v092);end,rethrow(v019),end
function [v023,v091]=f07(v098,v091),if nargin==0,v098=1;end
if nargin<2, v091=dbstack;end,v091(1:v098)=[];if ~isfield(v091,'file'),for v003=1:numel(v091),v099=v091(v003).name;if strcmp(v099(end),')')
v100=strfind(v099,'(');v101=v099( (v100(end)+1):(end-1) );v102=v099(1:(v100(end)-2));else,v102=v099;[v033,v101]=fileparts(v099);
end,[v033,v091(v003).file]=fileparts(v102);v091(v003).name=v101;end,end,persistent v103,if isempty(v103),v103=exist('OCTAVE_VERSION', 'builtin');end
if v103,for v003=1:numel(v091),[v033,v091(v003).file]=fileparts(v091(v003).file);end,end,v058=v091;v104='>';
v023=cell(1,numel(v058)-1);for v003=1:numel(v058),[v105,v058(v003).file,v106]=fileparts(v058(v003).file);if v003==numel(v058),v058(v003).file='';end
if strcmp(v058(v003).file,v058(v003).name),v058(v003).file='';end,if ~isempty(v058(v003).file),v058(v003).file=[v058(v003).file '>'];end
v023{v003}=sprintf('%c In %s%s (line %d)\n',v104,v058(v003).file,v058(v003).name,v058(v003).line);v104=' ';end,v023=horzcat(v023{:});end
function v107=f08(v108),if nargin==0,v109=f09;if isempty(v109),v109=f10;if isempty(v109)
error('HJW:getUTC:TimeReadFailed',['Both methods of retrieving the UTC timestamp failed.\nEnsure you ','have write access to the current folder and check your internet connection.'])
end,end,else,if v108==1,v109=f09;else,v109=f10;end,end,v110=v109/(24*60*60);v107=v110+datenum(1970,1,1);end
function v109=f09,persistent v111 v112 v113 v114 v115,if isempty(v111)
v113=['utc_time_' f16];v113=strrep(v113,['.' mexext],'');try v114=str2func(v113);catch,end,v112=fullfile(tempdir,'MATLAB','FileExchange','getUTC');try
if isempty(strfind([path ';'],[v112 ';'])),if ~exist(v112,'dir'),mkdir(v112);end,addpath(v112,'-end');end,catch,end,v115=5;v111={'#include "mex.h"';
'#include "time.h"';'';'/* Abraham Cohn,  3/17/2005 */';'/* Philips Medical Systems */';'';'void mexFunction(int nlhs, mxArray *plhs[], int nrhs,';
'                 const mxArray *prhs[])';'{';'  time_t utc;';'  ';'  if (nlhs > 1) {';'    mexErrMsgTxt("Too many output arguments");';'  }';'  ';
'  /* Here is a nice ref: www.cplusplus.com/ref/ctime/time.html */';'  time(&utc);';'  /* mexPrintf("UTC time in local zone: %s",ctime(&utc)); */';
'  /* mexPrintf("UTC time in GMT: %s",asctime(gmtime(&utc))); */';'  ';'  /* Create matrix for the return argument. */';
'  plhs[0] = mxCreateDoubleScalar((double)utc);';'   ';'}'};end,try v109=feval(v114);catch,if exist(['utc_time.' f16],'file')
v019=lasterror;rethrow(v019);end,v115=v115-1;if v115<0,v109=[];return,end,if f23(v112),v010=v112;else,v010=pwd;end,v116=cd(v010);
try if ~exist(fullfile(v010,[v113 '.c']),'file'),v117=fopen(fullfile(v010,[v113 '.c']),'w');for v074=1:numel(v111),fprintf(v117,'%s\n',v111{v074});
end,fclose(v117);end,try mex([v113 '.c']);catch,end,for v118={'c','o'},v102=fullfile(v010,[v113 '.' v118{1}]);if exist(v102,'file'),delete(v102),end
end,catch,end,cd(v116);if exist([v113 '.' mexext],'file'),v114=str2func(v113);v109=feval(v114);else,v109=[];end,end,end
function v109=f10,if ~f12,v109=[];return,end,for v031=1:3,try if exist('webread','file'),v015=webread('http://www.utctime.net/utc-timestamp');
else,v015=urlread('http://www.utctime.net/utc-timestamp');end,break,catch,end,end,try v015(v015==' ')='';v095='vartimestamp=';
v119=strfind(v015,v095)+numel(v095);v120=strfind(v015,';')-1;v120(v120<v119)=[];v109=str2double(v015(v119:v120(1)));catch,v109=[];end,end
function v121=f11(v122,v123,v124,v125,v126),persistent v127 v128 v129,if isempty(v127),v129=exist('OCTAVE_VERSION', 'builtin');v127=version;
v025=strfind(v127,'.');if numel(v025)~=1,v127(v025(2):end)='';v025=v025(1);end,v127=[str2double(v127(1:(v025-1))) str2double(v127((v025+1):end))];
v127=v127(1)+v127(2)/100;v127=round(100*v127);v128={ 'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;
'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;
'R2009a' 708;'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;
'R2012a' 714;'R2012b' 800;'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;
'R2015a' 805;'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;'R2020b',909};end,if v129,if nargin==2
warning('HJW:ifversion:NoOctaveTest',['No version test for Octave was provided.',char(10),'This function might return an unexpected outcome.'])
v130=ismember(v128(:,1),v123);
if sum(v130)~=1,warning('HJW:ifversion:NotInDict','The requested version is not in the hard-coded list.'),v121=NaN;return,else
v131=v128{v130,2};end,elseif nargin==4,[v122,v131]=deal(v124,v125);v131=0.1*v131+0.9*fix(v131);v131=round(100*v131);else,[v122,v131]=deal(v125,v126);
v131=0.1*v131+0.9*fix(v131);v131=round(100*v131);end,else,if isnumeric(v123),v131=0.1*v123+0.9*fix(v123);v131=round(100*v131);
else,v130=ismember(v128(:,1),v123);if sum(v130)~=1,warning('HJW:ifversion:NotInDict','The requested version is not in the hard-coded list.')
v121=NaN;return,else,v131=v128{v130,2};end,end,end,switch v122,case '==', v121= v127 == v131;
case '<' , v121= v127 < v131;case '<=', v121= v127 <= v131;case '>' , v121= v127 > v131;case '>=', v121= v127 >= v131;end,end
function [v132,v133]=f12,v121=f15;if isempty(v121),v132=0;v133=0;else,if v121,[v132,v133]=f13;else,[v132,v133]=f14;end,end,end
function [v132,v133]=f13,try v004=now;if exist('webread','file')
v023=webread('http://google.com');else,v023=urlread('http://google.com');end,v132=1;v133=(now-v004)*24*3600*1000;catch,v132=0;v133=0;end,end
function [v132,v133]=f14,if ispc,try [v080,v134]=system('ping -n 1 8.8.8.8');
v135=v134(strfind(v134,' = ')+3);v135=v135(1:3);if ~strcmp(v135,'110'),error('trigger error'),else,v132=1;[v119,v120]=regexp(v134,' [0-9]+ms');
v133=v134((v119(1)+1):(v120(1)-2));v133=str2double(v133);end,catch,v132=0;v133=0;end,elseif isunix,try [v080,v134]=system('ping -c 1 8.8.8.8');
v100=regexp(v134,', [01] ');if v134(v100+2)~='1',error('trigger error'),else,v132=1;[v119,v120]=regexp(v134,'=[0-9.]+ ms');
v133=v134((v119(1)+1):(v120(1)-2));v133=str2double(v133);end,catch,v132=0;v133=0;end,else,error('How did you even get Matlab to work?'),end,end
function [v121,v132,v133]=f15,persistent v136,if ~isempty(v136)
v121=v136;return,end,[v132,v133]=f14;if v132,v136=false;v121=false;return,end,[v132,v133]=f13;if v132,v136=true;v121=true;return,end,v121=[];end
function v118=f16
v131=version;v100=strfind(v131,'.');v131(v100(2):end)='';v131=['v' strrep(v131,'.','_')];if ~exist('OCTAVE_VERSION', 'builtin'),v045=computer;
else,v137=computer;v137=v137(1:(strfind(v137,'-')-1));if ispc,if strcmp(v137,'x86_64') ,v045= 'win_64';elseif strcmp(v137,'i686'),v045= 'win_i686';
elseif strcmp(v137,'x86') ,v045= 'win_x86';else ,v045=['win_' v137];end,elseif isunix && ~ismac,if strcmp(v137,'i686') ,v045= 'lnx_i686';
elseif strcmp(v137,'x86_64'),v045= 'lnx_64';else ,v045=['lnx_' v137];end,elseif ismac,if strcmp(v137,'x86_64'),v045= 'mac_64';
else ,v045=['mac_' v137];end,end,end,v045=strrep(strrep(v045,'.',''),'-','');v118=[v131 '_' v045 '.' mexext];end
function v001=f17(v138,v139,v140)
v001=v138(:)';if numel(v139)==0,v130=false(size(v138));elseif numel(v140)>numel(v139),error('not implemented (padding required)')
else,v130=true(size(v138));for v003=1:numel(v139),v141=find(v138==v139(v003));v141=v141-v003+1;v141(v141<1)=[];v142=false(size(v130));v142(v141)=true;
v130= v130 & v142;if ~any(v130),break,end,end,end,v141=find(v130);if ~isempty(v141),for v003=1:numel(v140),v001(v141+v003-1)=v140(v003);
end,if numel(v140)==0,v003=0;end,if numel(v139)>v003,v141=v141(:);v143=(v003+1):numel(v139);v075=f01(v141,v143-1);v001(v075(:))=[];end,end,end
function v023=f18(v016,v022,v021)
if nargin==1,v022=struct;v021=false;end,persistent v144,if isempty(v144),v144 = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
persistent v145,if isempty(v145),if v144,v146='Octave';else,v146='Matlab';end
v145=struct('identifier','HJW:readfile:ReadFail','message',sprintf('%s could not read the file %s.\nThe file doesn''t exist or is not readable.',v146,v016));
end,v117=fopen(v016,'rb');if v117<0,f06(v022,v145),end,v015=fread(v117,'uint8');
fclose(v117);v015=v015.';try v026=true;v027=f28(v015);catch,v019=lasterror;if strcmp(v019.identifier,'HJW:UTF8_to_unicode:notUTF8'),v026=false;
if v021,f06(v022,'HJW:readfile:notUTF8','The provided file "%s" is not a correctly encoded UTF-8 file.',v016),end,else,f06(v022,v019),end,end
if v144,if v026,v023=v027;else,try v023=fileread(v016);catch,f06(v022,v145),end,v023=f05(v023);end,else,if ispc,if v026,v023=v027;else,if f11('<',7)
try v023=fileread(v016);catch,f06(v022,v145),end,v023=f05(v023);else,try v023=fileread(v016);catch,f06(v022,v145),end,end,end,else,if v026,v023=v027;
else,v023=f05(v015);end,end,end,v023(v023==13)=[];if numel(v023)>=1 && double(v023(1))==65279,v023(1)=[];end,v023=f25(v023);v023=f02(v023);end
function v023=f19(v147,v020,v148)
try v149=false;v150=f24('readfile_from_URL_tmp_','.txt');try if v020,v150=urlwrite(v147,v150);else,v150=websave(v150,v147);end
v023=f18(v150);catch,v149=true;end,try if exist(v150,'file'),delete(v150);end,catch,end
if v149,error('revert to urlread'),end,catch,try if v020,v023=urlread(v147);else,v023=webread(v147);end,catch,f06(v148,lasterror),end,end,end
function [v017,v090,v019]=f20(v016,varargin),v017=false;v090=struct;v019=struct('identifier','','message','');
if ~( isa(v016,'char') || isa(v016,'string') ) ||( isa(v016,'string') && numel(v016)~=1 ) ||( isa(v016,'char') && numel(v016)==0 )
v019.identifier='HJW:readfile:IncorrectInput';v019.message='The file name must be a non-empty char or a scalar string.';
return,end,persistent v151,if isempty(v151),v005.split = f11('<','R2007b','Octave','<',4);
v005.allows_https=f11('>',0,'Octave','<',0);v151.legacy=v005;v151.UseURLread=isempty(which('webread')) || isempty(which('websave'));
v151.print_2_con=true;v151.print_2_fid.boolean=false;v151.print_2_obj.boolean=false;v151.err_on_ANSI=false;end,if nargin==1
v090=v151;v017=true;return,end,v152=nargin==2 && isa(varargin{1},'struct');v153=mod(nargin,2)==1 &&all(cellfun('isclass',varargin(1:2:end),'char'));
if ~( v152 || v153 ),v019.message=['The second input (options) is expected to be either a struct, ','or consist of Name,Value pairs.'];
v019.identifier='HJW:readfile:incorrect_input_options';
return,end,if v153,for v003=1:2:numel(varargin),try v090.(varargin{v003})=varargin{v003+1};catch,v019.message='Parsing of Name,Value pairs failed.';
v019.identifier='HJW:readfile:incorrect_input_NameValue';return,end,end,else,v090=varargin{1};end,v150=fieldnames(v090);for v141=1:numel(v150)
v154=v150{v141};v155=v090.(v154);v019.identifier=['HJW:readfile:incorrect_input_opt_' lower(v154)];switch v154,case 'UseURLread'
[v156,v155]=f22(v155);if ~v156,v019.message='UseURLread should be either true or false';return,end,v090.UseURLread=v155 || v151.UseURLread;
case 'err_on_ANSI',[v156,v155]=f22(v155);if ~v156,v019.message='err_on_ANSI should be either true or false';return,end,v090.show_UTF_err=v155;
case 'print_to_fid',if isempty(v155),v090.print_2_fid.boolean=false;else,v090.print_2_fid.boolean=true;try v157=ftell(v155);catch,v157=-1;end
if v155~=1 && v157==-1,v019.message=['Invalid print_to_fid parameter: ','should be a valid file identifier or 1'];
return,end,v090.print_2_fid.fid=v155;end,case 'print_to_obj',if isempty(v155)
v090.print_2_obj.boolean=false;else,v090.print_2_obj.boolean=true;try v051=get(v155,'String');set(v155,'String','');set(v155,'String',v051);
v090.print_2_obj.obj=v155;catch,v019.message=['Invalid print_to_obj parameter: ','should be a handle to an object with a writeable String property'];
return,end,end,case 'print_to_con',[v156,v090.print_2_con]=f22(v155);if ~v156,v019.message='print_to_con should be either true or false';
return,end,end,end,if ~isfield(v090,'print_2_con') &&( isfield(v090,'print_2_fid') || isfield(v090,'print_2_obj') ),v090.print_2_con=false;
end,v150=fieldnames(v151);for v141=1:numel(v150),if ~isfield(v090,v150(v141)),v090.(v150{v141})=v151.(v150{v141});end,end,v017=true;v019=[];end
function v023=f21(v023),if exist('strtrim','builtin'),v023=strtrim(v023);else,if numel(v023)==0,return,end,v130=isspace(v023);if v130(end)
v075=find(~v130);if isempty(v075),v023='';return,end,v023((v075(end)+1):end)='';end,if isempty(v023),return,end,if v130(1),v075=find(~v130);
v023(1:(v075(1)-1))='';end,end,v158=inf;while v158~=0,v159=length(v023);v023=strrep(v023,'  ',' ');v160=length(v023);v158=v159-v160;end,end
function [v161,v162]=f22(v162),persistent v163,if isempty(v163),v163={true,false;1,0;'on','off'};
try v163(end+1,:)=eval('{"on","off"}');catch,end,end,v161=true;try for v003=1:size(v163,1),for v078=1:2,if isequal(v162,v163{v003,v078})
v162=v163{1,v078};return,end,end,end,if isa(v162,'matlab.lang.OnOffSwitchState'),v162=logical(v162);return,end,catch,end,v161=false;end
function v121=f23(v010),if ~( isempty(v010) || exist(v010,'dir') ),v121=false;return
end,v150='';while isempty(v150) || exist(v150,'file'),[v033,v150]=fileparts(f24('write_permission_test_','.txt'));v150=fullfile(v010,v150);
end,try v117=fopen(v150,'w');fprintf(v117,'test');fclose(v117);delete(v150);v121=true;catch,try delete(v150);catch,end,v121=false;end,end
function v023=f24(v164,v118)
if nargin<1,v164='';end,if ~isempty(v164),v164=[v164 '_'];end,if nargin<2,v118='';else,if ~strcmp(v118(1),'.'),v118=['.' v118];end,end
v023=tempname;[v165,v010]=fileparts(v023);v023=fullfile(v165,[v164 v010 v118]);end
function v023=f25(v166,v167),persistent v144,if isempty(v144),v144 = exist('OCTAVE_VERSION', 'builtin') ~= 0;end,if nargin==1,v167=~v144;
end,if v167,if all(v166<65536),v023=uint16(v166);v023=reshape(v023,1,numel(v023));else,[v168,v033,v169]=unique(v166);v023=cell(1,numel(v166));
for v003=1:numel(v168),v170=f26(v168(v003));v170=uint16(v170);v023(v169==v003)={v170};end,v023=cell2mat(v023);end,if ~v144,v023=char(v023);
end,else,if all(v166<128),v023=char(v166);v023=reshape(v023,1,numel(v023));else,[v168,v033,v169]=unique(v166);v023=cell(1,numel(v166));
for v003=1:numel(v168),v170=f27(v168(v003));v170=uint8(v170);v023(v169==v003)={v170};end,v023=cell2mat(v023);v023=char(v023);end,end,end
function v023=f26(v166)
if v166<65536,v023=v166;return,end,v171=double(v166)-65536;v171=dec2bin(v171,20);v023=bin2dec(['110110' v171(1:10);'110111' v171(11:20)]).';end
function v023=f27(v166)
if v166<128,v023=v166;return,end,persistent v172,if isempty(v172),v172=struct;v172.limits.lower=hex2dec({'0000','0080','0800', '10000'});
v172.limits.upper=hex2dec({'007F','07FF','FFFF','10FFFF'});v172.scheme{2}='110xxxxx10xxxxxx';v172.scheme{2}=reshape(v172.scheme{2}.',8,2);
v172.scheme{3}='1110xxxx10xxxxxx10xxxxxx';v172.scheme{3}=reshape(v172.scheme{3}.',8,3);v172.scheme{4}='11110xxx10xxxxxx10xxxxxx10xxxxxx';
v172.scheme{4}=reshape(v172.scheme{4}.',8,4);for v134=2:4,v172.scheme_pos{v134}=find(v172.scheme{v134}=='x');
v172.bits(v134)=numel(v172.scheme_pos{v134});end,end,v173=find(v172.limits.lower<v166 & v166<v172.limits.upper);
v023=v172.scheme{v173};v174=v172.scheme_pos{v173};v134=dec2bin(v166,v172.bits(v173));v023(v174)=v134;v023=bin2dec(v023.').';end
function [v166,v026,v175]=f28(v176,v148),if nargin<2,v148=[];end,v177= nargout==1 ;v176=uint32(v176);[v175,v036,v019]=f29(v176,v177);
if strcmp(v036,'success'),v026=true;v166=v175;elseif strcmp(v036,'error'),v026=false;if v177,f06(v148,v019),end,v166=v176;end,end
function [v176,v036,v019]=f29(v176,v177),v036='success';v019=struct('identifier','HJW:UTF8_to_unicode:notUTF8','message','Input is not UTF-8.');
persistent v144,if isempty(v144),v144 = exist('OCTAVE_VERSION', 'builtin') ~= 0;end,if any(v176>255),v036='error';
if v177,return,end,elseif all(v176<128),return,end,for v173=4:-1:2,v162=bin2dec([repmat('1',1,v173) repmat('0',1,8-v173)]);v178=v176>=v162 & v176<256;
if any(v178),v178=find(v178);v178=v178(:).';if numel(v176)<(max(v178)+v173-1),v036='error';if v177,return,end,v178( (v178+v173-1)>numel(v176) )=[];end
if ~isempty(v178),v075=f01(v178 , (0:(v173-1)).' );v075=v075.';v178=v176(v075);end,else,v178=[];end,v179=[repmat('1',1,v173-1) repmat('10',1,v173)];
v180=unique([1:(v173+1) 1:8:(8*v173) 2:8:(8*v173)]);if numel(v178)>0,v178=unique(v178,'rows');v181=mat2cell(v178,ones(size(v178,1),1),v173);
for v003=1:numel(v181),v182=dec2bin(double(v181{v003}))';if ~strcmp(v179,v182(v180)),v036='error';if v177,return,end
continue,end,v182(v180)='';if ~v144,v183=uint32(bin2dec(v182 ));else,v183=uint32(bin2dec(v182.'));end,v176=f17(v176,v181{v003},v183);end,end,end,end
function f30(v090,varargin),if isempty(v090),v090=struct;end
if ~isfield(v090,'con'),v090.con=false;end,if ~isfield(v090,'fid'),v090.fid.boolean=false;end,if ~isfield(v090,'obj'),v090.obj.boolean=false;end
if nargin==2 || ~isempty(strfind(varargin{1},'%')),[v093,v094]=deal('',varargin{1});if nargin>3,v096=varargin(2:end);v094=sprintf(v094,v096{:});
end,else,[v093,v094]=deal(varargin{1},varargin{2});if nargin>3,v096=varargin(3:end);v094=sprintf(v094,v096{:});end,end,if v090.con,if ~isempty(v093)
warning(v093,'%s',v094),else,warning(v094),end,else,if ~isempty(v093),lastwarn(v094,v093);else,lastwarn(v094),end,end,if v090.obj.boolean
v097=v094;while v097(end)==10,v097(end)=[];end,if any(v097==10),v097=regexp_outkeys(['Warning: ' v097],char(10),'split');else,v097=['Warning: ' v097];
end,set(v090.obj.obj,'String',v097),end,if v090.fid.boolean,v098=2;v092=f07(v098);fprintf(v090.fid.fid,'Warning: %s\n%s',v094,v092);end,end
function v107=f31,v107=f08(2);if isempty(v107),v107=f08(1);if isempty(v107),v107=0;end,end,end
function [v184,v185]=f32(v030)
shift=[ 1000, 0, 0, 0, 0, 0;100 , 0, 0, 0, 0, 0;10 , 0, 0, 0, 0, 0;1 , 0, 0, 0, 0, 0;0 ,10, 0, 0, 0, 0;0 , 1, 0, 0, 0, 0;0 , 0,10, 0, 0, 0;0 , 0, 1, 0, 0, 0;0 , 0, 0,10, 0, 0;0 , 0, 0, 1, 0, 0;0 , 0, 0, 0,10, 0;0 , 0, 0, 0, 1, 0;0 , 0, 0, 0, 0,10;0 , 0, 0, 0, 0, 1];
v187=char(zeros(1,14)+'0');
v187(1:numel(v030))=v030;v187={v187(1:4),v187(5:6),v187(7:8),v187(9:10),v187(11:12),v187(13:14)};v187=str2double(v187);v187(1:3)=max(v187(1:3),1);
v188=num2cell(v187);v188=datenum(v188{:});v122=datestr(v188,'yyyymmddTHHMMSS');v122(9)='';if ~strcmp(v122(1:numel(v030)),v030)
error('incorrect date_part format'),end,v189=v187+shift(numel(v030),:);v190=[31 28 31 30 31 30 31 31 30 31 30 31];
v191=v189(1);if rem(v191,4)==0 && (rem(v191,100)~=0 || rem(v191,400)==0),v190(2) = 29;end,v190=v190(min(12,v189(2)));
v192=find(v189>[inf 12 v190 23 59 59])+1;if ~isempty(v192),v189(v192:end)=inf;else,v189(end)=-1;end,v189=min(v189,[inf 12 v190 23 59 59]);
v189=num2cell(v189);v189=datenum(v189{:});v185=(v188+v189)/2;v185=datestr(v185,'yyyymmddTHHMMSS');v185(9)='';v184=[v188 v189];end
function [v017,v090,v019]=f33(v016,v147,varargin),v017=false;v090=struct;
v019=struct('identifier','','message','');if ~ischar(v016) || numel(v016)==0,v019.message='The first input (filename) is not char and/or empty.';
v019.identifier='HJW:WBM:incorrect_input_filename';return,end,if ~ischar(v147) || numel(v147)==0
v019.message='The second input (url) is not char and/or empty.';v019.identifier='HJW:WBM:incorrect_input_url';return
end,persistent v151 v193,v194=false;if isempty(v151),v151.date_part='2';v151.tries=[5 4 4];v151.response={ 'tx',404,'load';'txtx',[404 404],'save';
'tx',403,'save';'t2t2',[403 403],'exit';'tx',429,'pause_retry' };v151.ignore=4080;v151.verbose=3;v151.m_date_r=1;v151.flag='*';
v151.UseLocalTime=false;v151.UseURLwrite=isempty(which('websave'));v151.err429=struct('CountsAsTry',false,'TimeToWait',15,'PrintAtVerbosityLevel',3);
v151.print_2_con=true;v151.print_2_fid.boolean=false;v151.print_2_obj.boolean=false;[v184,v185]=f32(v151.date_part);
v151.date_part=v185;v151.date_bounds.datenum=v184;v155={datestr(v184(1),'yyyymmddTHHMMSS'),datestr(v184(2),'yyyymmddTHHMMSS')};
v155{1}(9)='';v155{2}(9)='';v151.date_bounds.double=[str2double(v155{1}) str2double(v155{2})];v193=f35;end,if nargin==2
v090=v151;v017=true;return,end,v152=nargin==3 && isa(varargin{1},'struct');v153=mod(nargin,2)==0 &&all(cellfun('isclass',varargin(1:2:end),'char'));
if ~( v152 || v153 ),v019.message=['The third input (options) is expected to be either a struct, ','or consist of Name,Value pairs.'];
v019.identifier='HJW:WBM:incorrect_input_options';return,end,if v153,for v003=1:2:numel(varargin),try v090.(varargin{v003})=varargin{v003+1};
catch,v019.message='Parsing of Name,Value pairs failed.';v019.identifier='HJW:WBM:incorrect_input_NameValue';return,end,end,else,v090=varargin{1};
end,v150=fieldnames(v090);for v141=1:numel(v150),v154=v150{v141};v155=v090.(v154);v019.identifier=['HJW:WBM:incorrect_input_opt_' lower(v154)];
switch v154,case 'date_part',v195=true;if ~ischar(v155) || numel(v155)==0 || numel(v155)>14 || any(v155<48 & v155>57),v195=false;
end,try [v184,v185]=f32(v155);catch,v195=false;end,if ~v195,v019.message='The value of options.date_part is empty or not a valid numeric char.';
return,end,v090.date_part=v185;v090.date_bounds.datenum=v184;case 'tries',if ~isnumeric(v155) || numel(v155)~=3 || any(isnan(v155))
v019.message=['The value of options.tries has an incorrect format.',char(10),'The value should be a numeric vector with 3 integer elements.'];
return,end,case 'response',if f34(v155),v019.message='The value of options.response has an incorrect format.';
return,end,case 'ignore',if ~isnumeric(v155) || numel(v155)==0 || any(isnan(v155))
v019.message=['The value of options.ignore has an incorrect format.',char(10),'The value should be a numeric vector with HTML error codes.'];
return,end,case 'verbose',if ~isnumeric(v155) || numel(v155)~=1 || double(v155)~=round(double(v155))
v019.message='The value of options.verbose is not an integer scalar.';return,end,case 'm_date_r',if ~ischar(v155) || numel(v155)==0
v019.message='Options.m_date_r should be ''ignore'', ''warning'', or ''error''.';return,end,v194=true;switch lower(v155),case 'ignore'
v155=0;case 'warning',v155=1;case 'error',v155=2;otherwise,v019.message='Options.m_date_r should be ''ignore'', ''warning'', or ''error''.';
return,end,v090.m_date_r=v155;case 'flag',if ischar(v155)&&numel(v155)~=0&&~ismember({v155},{'*','id','js','cs','im','if','fw'})
v019.message='Invalid flag. Must be a char with *, id, js, cs, im, fw, or if.';
return,end,case 'UseLocalTime',[v156,v090.UseLocalTime]=f22(v155);if ~v156,v019.message='UseLocalTime should be either true or false';
return,end,case 'UseURLwrite',[v156,v155]=f22(v155);if ~v156,v019.message='UseURLwrite should be either true or false';
return,end,v090.UseURLwrite=v155 || v151.UseURLwrite;case 'err429',if ~isa(v155,'struct'),v019.message='The err429 parameter should be a struct.';
return,end,v090.err429=v151.err429;v196=fieldnames(v155);for v003=1:numel(v196),v197=v155.(v196{v003});switch lower(v196{v003})
case 'countsastry',[v156,v197]=f22(v197);if ~v156,v019.message=['Invalid field CountsAsTry in the err429 parameter: ','should be a logical scalar.'];
return,end,v090.err429.CountsAsTry=v197;case 'timetowait'
if ~isnumeric(v197) || numel(v197)~=1,v019.message=['Invalid field TimeToWait in the err429 parameter: ','should be a numeric scalar.'];
return,end,v090.err429.TimeToWait=double(v197);case 'printatverbositylevel',if ~isnumeric(v197) || numel(v197)~=1 ||double(v197)~=round(double(v197))
v019.message=['Invalid field PrintAtVerbosityLevel in the err429 ','parameter: should be a scalar double integer.'];
return,end,v090.err429.PrintAtVerbosityLevel=v197;
otherwise,warning('HJW:WBM:NameValue_not_found',['Name,Value pair not ','recognized during parsing of err429 parameter: %s'],v196{v003});
end,end,case 'print_to_fid',if isempty(v155),v090.print_2_fid.boolean=false;else,v090.print_2_fid.boolean=true;try v157=ftell(v155);catch,v157=-1;end
if v155~=1 && v157==-1,v019.message=['Invalid print_to_fid parameter: ','should be a valid file identifier or 1.'];
return,end,v090.print_2_fid.fid=v155;end,case 'print_to_obj',if isempty(v155)
v090.print_2_obj.boolean=false;else,v090.print_2_obj.boolean=true;try v051=get(v155,'String');set(v155,'String','');set(v155,'String',v051);
v090.print_2_obj.obj=v155;catch,v019.message=['Invalid print_to_obj parameter: ','should be a handle to an object with a writeable String property'];
return,end,end,case 'print_to_con',[v156,v090.print_2_con]=f22(v155);if ~v156,v019.message='print_to_con should be either true or false';
return,end,otherwise,v019.message=sprintf('Name,Value pair not recognized: %s',v154);v019.identifier='HJW:WBM:incorrect_input_NameValue';return,end
end,if ~isfield(v090,'print_2_con') &&( isfield(v090,'print_2_fid') || isfield(v090,'print_2_obj') ),v090.print_2_con=false;end,v150=fieldnames(v151);
for v141=1:numel(v150),if ~isfield(v090,v150(v141)),v090.(v150{v141})=v151.(v150{v141});end,end,try v198=f31;catch,errordlg(v193.errordlg_text{:})
v019=v193.ME;return,end,if v090.UseLocalTime,v199=v198-now;v199=round(v199*24*4)/(24*4);v090.date_bounds.datenum=v090.date_bounds.datenum-v199;
v184=v090.date_bounds.datenum;else,v184=v090.date_bounds.datenum;end,v155={datestr(v184(1),'yyyymmddTHHMMSS'),datestr(v184(2),'yyyymmddTHHMMSS')};
v155{1}(9)='';v155{2}(9)='';v090.date_bounds.double=[str2double(v155{1}) str2double(v155{2})];
if ~( v184(1)<=v198 && v198<=v184(2) ),v090.tries(2)=0;end,if v090.m_date_r==2 && ~strcmp(v090.flag,'*')
v019.message=['m_date_r set to ''error'' and the flag set to something other than ''*'' will',' by definition',char(10),'cause an error, as the downloaded pages will not contain',' any dates.',char(10),'See the help text for a work-around for images.'];
v019.identifier='HJW:WBM:IncompatibleInputs';return,end,if ~v194 && ~strcmp(v090.flag,'*'),v090.m_date_r=0;end,v017=true;v019=[];end
function v200=f34(v032),v200=false;if ~iscell(v032) || isempty(v032) || size(v032,2)~=3,v200=true;return
end,for v201=1:size(v032,1),v155=v032{v201,1};v202=numel(v155(2:2:end));if ~ischar(v155) || numel(v155)==0 || ~all(ismember(v155(2:2:end),'12x'))
v200=true;return,end,v155=v032{v201,2};if ~isa(v155,'double') || numel(v155)~=v202
v200=true;return,end,v155=v032{v201,3};if ~ischar(v155) || ~ismember({v155},{'load','save','exit','pause_retry'}),v200=true;return,end,end,end
function v203=f35,persistent v193,if isempty(v193) || v193.boolean
v193=struct;v193.boolean=true;v193.errordlg_text={{['The WBM function relies on a compiled c function to get the current time in UTC '.';
'when offline or for Matlab 6.5. The automatic compiling mechanism appears to have '.';
'failed. Please manually make sure utc_time.c is compiled for your system.'.'].'},'utc_time.c compile failed'};
	v193.ME.message='Retrieval of UTC failed.';v193.ME.identifier='HJW:WBM:UTC_missing';try if abs(f31-now)>(14.1/24)
error('HJW:WBM:getUTCfail',['The getUTC function failed,',char(10),'or the ','difference between now() and the UTC time is more than a day.',char(10),'The ','getUTC function requires either internet access or folder write access.'])
end,catch,errordlg(v193.errordlg_text{:}),error(v193.ME.identifier,v193.ME.message),end,v193.boolean=false;end,v203=v193;end