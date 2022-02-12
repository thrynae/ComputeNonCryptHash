function [hash,varargout]=ComputeNonCryptHash(data,varargin)
%Compute a non-cryptographic hash
%
% This function is intended to be fast, but without requiring a Java or mex implementation to do
% the actual hashing. It was *not* checked for any security flaws and is therefore probably
% vulnerable to most attacks.
% Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like
% storing passwords.
%
% Syntax:
%   hash=ComputeNonCryptHash(data)
%   hash=ComputeNonCryptHash(___,HashLength)
%   hash=ComputeNonCryptHash(___,VersionFlag)
%   hash=ComputeNonCryptHash(___,options)
%   hash=ComputeNonCryptHash(___,Name,Value)
%
% Input/output arguments:
% hash:
%   The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
% data:
%   The data to be hashed. Most common data types are allowed: uint*, int*, char, cell, struct,
%   double, or single (string is cast to a cell array of chars). The contents of the nested data
%   types (i.e. cell and struct) must also be one of the mentioned data types.
% HashLength/VersionFlag:
%   See below.
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%   Parameters related to warning/error redirection will be parsed first.
%
% Name,Value parameters:
%   HashLength:
%      The length of the hash (the number of bits). This value must be a multiple of 16. The
%      default is 256 bits. Depending on your input 64 bits might have some collisions, but 64 bits
%      and higher should be safe.
%   VersionFlag:
%      Either '-v1', '-v2'. This is provided for backwards compatibility. Version 1 of this
%      function has many hash collisions for scalar doubles and attempts to cast strings to chars,
%      instead of casting to a cell array of chars. Version 2 also decodes the UTF-8 chars from
%      Octave and re-encodes them with UTF-16. That way the output is stable for the Unicode code
%      points.
%   print_to_con:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;] if either print_to_fid, print_to_obj, or
%      print_to_fcn is specified then [default=false]
%   print_to_fid:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_obj:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.2.0                                                         |%
%|  Date:    2022-02-12                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - On releases prior to R2010b int64 and uint64 are converted to double as an intermediate step.
%   That means values larger than flintmax (i.e. 2^52) may be rounded, changing the hash. Version 1
%   applies this conversion on all releases, making it stable in this respect.
% - Attributes like the sparsity are ignored and gather will be called on tall/gpuArray objects.

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
end
if nargout>=2
    varargout=cell(nargout-1,1);
end

% To allow fast processing of an internal call, skip the input parsing.
if nargin==2 && isa(varargin{1},'struct') && varargin{1}.SkipInputParse
    opts=varargin{1};
else
    [success,opts,ME]=ComputeNonCryptHash_parse_inputs(varargin{:});
    if success
        opts=ComputeNonCryptHash_DefaultsByVersion(opts);
        % Undocumented shortcut: return the parsed options for faster future calls.
        if nargout>=2,hash='';opts.SkipInputParse=true;varargout{1}=opts;return,end
    else
        % If the parsing of print_to failed (which is tried first), the default will be used.
        error_(opts.print_to,ME)
    end
end
HashLength=opts.HashLength;
Version=opts.Version;

try ME=[]; %#ok<NASGU>
    % Convert the input to an uint16 array (Nx1).
    data=cast_to_uint16_vector(data,opts);
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    if strcmp(ME.identifier,'MATLAB:nomem')
        % Rethrow memory error.
        error_(opts.print_to,ME)
    else
        if isfield(opts,'debug') && opts.debug
            msg=sprintf('\n[original error: %s %s]',ME.identifier,ME.message);
        else
            msg='';
        end
        error_(opts.print_to,'HJW:ComputeNonCryptHash:UnwindFailed',...
            ['The nested input contains an unsupported data type.' msg])
    end
end

% Extend to a multiple of HashLength bits. Padding with zeros is generally not advised, and the
% performance penalty for this extension (compared to padding with zeros) should be negligible.
if mod(numel(data),HashLength/16)
    extra=uint16(1:HashLength/16).'; extra(1:mod(numel(data),HashLength/16))=[];
    data=[data;extra];
end

% Add perturbation to the data and convert to 16xN logical. Then further perturb the intermediate
% result by circshifting every column (by col_index-1 positions) and doing an XOR in blocks with
% itself (by reshaping and transposing).
if Version==1
    data=ComputeNonCryptHash_shuffle_uint16(data);
    data=ComputeNonCryptHash_uint16_to_logical(data);
    data=xor(data,reshape(data,[],16).');
else
    data=ComputeNonCryptHash_shuffle_uint16(data);
    data=ComputeNonCryptHash_uint16_to_logical(data);
    data=circshift_by_col(data);
    %data=circshift_by_col(data.').';
    %data=xor(data,not(reshape(data,[],16).'));
    %data=checker_shift(data);
end

% Reshape to HashLength cols and collapse the key size down to the hash length by counting the
% number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if opts.isSaltCall
    hash=data; % Return uint16 for the salting.
    return
end

% Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
% few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data,opts);

% Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=circshift_by_col(data)
% Circshift every column by col_index-1 positions.
persistent LUT
sz=size(data);
if isempty(LUT) || any(size(LUT)<sz) || isempty(LUT{sz(1),sz(2)})
    % Keep a record of ind, which speeds up similar sizes.
    [x,y]=meshgrid(1:size(data,2),1:size(data,1));
    z=mod(x+y-2,size(data,1))+1;
    ind=sub2ind(size(data),z,x);
    if prod(sz)<=1000 % To prevent a memory-hog, only keep ind for small sizes.
        LUT{sz(1),sz(2)}=ind;
    end
else
    ind=LUT{sz(1),sz(2)};
end
data=data(ind);
end
function data=ComputeNonCryptHash_add_salt(data,options)
% Apply a few transformations to the hash to increase the spread.
% If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
% A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
SaltOptions=options;
% Overwrite two options:
SaltOptions.Version=1;         SaltOptions.HashLength=saltHashLength;
% Avoid an infinite recursion by using two undocumented switches:
SaltOptions.SkipInputParse=1;  SaltOptions.isSaltCall=1;
% Do a recursive call to perturb the data.
salt=ComputeNonCryptHash(data,SaltOptions);
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
if options.Version>1
    salt=salt(end:-1:1);
end
data=mod(double(data).*double(salt),1+2^16);
data=uint16(data);
end
function hash=ComputeNonCryptHash_dec2hex(data)
% Look up the precomputed dec2hex for faster conversion.
persistent LUT
if isempty(LUT)
    LUT=upper(dec2hex(0:(-1+2^16),4)); % The default should already upper case.
end
data=double(data)+1;% plus() is not defined for uint16 on older Matlab releases.
hash=LUT(data,:);
end
function data=ComputeNonCryptHash_logical_to_uint16(data)
if mod(numel(data),16) % Pad to 16 bits with 0.
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch % No implicit expansion.
    %times() is not defined for uint16 on older Matlab releases.
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';
end
function data=ComputeNonCryptHash_shuffle_uint16(data)
% Input should be uint16.
base=65537;%base=(1+(2^16));
key=479001600;%key=1*2*3*4*5*6*7*8*9*10*11*12;
data = uint16(mod(double(data) * key , base));
end
function data=ComputeNonCryptHash_shuffle_uint16_inv(data)
base=65537;%base=(1+(2^16));
%key=1*2*3*4*5*6*7*8*9*10*11*12;
% %Solution suggested by John D'Errico, https://www.mathworks.com/matlabcentral/answers/81859
% [G,C]=gcd(key,base);invKey=mod(C,base);
invKey=1919;
data=uint16(mod(double(data) * invKey,base));
end
function data=ComputeNonCryptHash_uint16_to_logical(data)
% uint16 Nx1 vector in, logical 16xN array out
persistent LUT
if isempty(LUT)
    LUT=dec2bin(0:(-1+2^16))=='1';
    LUT=LUT.';
end
data=double(data)+1;% The plus() function is not defined for uint16 on old Matlab releases.
data=LUT(:,data);
end
function out=bsxfun_plus(in1,in2)
%Implicit expansion for plus(), but without any input validation.
try
    out=in1+in2;
catch
    try
        out=bsxfun(@plus,in1,in2);
    catch
        sz1=size(in1);                    sz2=size(in2);
        in1=repmat(in1,max(1,sz2./sz1));  in2=repmat(in2,max(1,sz1./sz2));
        out=in1+in2;
    end
end
end
function data=cast_to_uint16_vector(data,options)
%Linearize the input data and convert it to a uint16 vector
if isa(data,'uint16')
    % Append the array size and type to make it influence the hash.
    c='uint16';sz=size(data).';
    data=reshape(data,[],1);% Linearize the data.
    data=[data;uint16(c.');uint16(mod(sz,2^16))];
    return
end
data=cast_to_uint16_vector__cell({data},options);
data([end-1 end])=[];% Remove the [1 1] that is added because of the wrapping in a cell.
end
function data=cast_to_uint16_vector__cell(data,options)
sz=size(data).';data=data(:);
for n=1:numel(data)
    if numel(data{n})==0
        c=double(class(data{n})');
        data{n}=uint16([0;c;size(data{n})']);
        continue
    end
    switch class(data{n})
        case {'double','single'}
            data{n}=cast_to_uint16_vector__floats(data{n});
        case 'logical'
            data{n}=cast_to_uint16_vector__logical(data{n});
        case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
            data{n}=cast_to_uint16_vector__integer(data{n},options);
        case 'char'
            data{n}=cast_to_uint16_vector__char(data{n},options);
        case 'string'
            data{n}=cast_to_uint16_vector__string(data{n},options);
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n},options);
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n},options);
        case {'gpuArray','tall'}
            data{n}=cast_to_uint16_vector__cell({gather(data{n})},options);
        otherwise
            error_(options.print_to,'HJW:cast_to_uint16_vector:nosupport',...
                'Unsupported data type in nested variable')
    end
end
data=cell2mat(data); % Merge all cell contents.
data=[data;uint16(mod(sz,2^16))]; % Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__floats(data)
sz=size(data).';c=class(data); % The rest of the function treats singles as double.

% Convert to a uint64, separate it into 4 words and everything merge into a vector.
data=reshape(data,size(data,1),[]); % Reshape ND to 2D.
[bit64,bit4]=typecast_double_uint64(double(data));
bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;bit4=bit4.';
bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;bit3=bit3.';
bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;bit2=bit2.';
bit1      =mod(bit64,2^16);                                        bit1=bit1.';
data=[bit1;bit2;bit3;bit4];
data=uint16(data(:));

% Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz,2^16))];
end
function data=cast_to_uint16_vector__logical(data)
sz=size(data).';data=data(:);

if mod(numel(data),16) % Pad to 16 bits with 0.
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch % No implicit expansion.
    % Cast to double, because times() is not defined for uint16 on old Matlab releases.
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';

data=[data;uint16(mod(sz,2^16))]; % Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__integer(data,options)
% Large values (>2^52) will not have integer precision due to a conversion to double.
% This conversion is done pre-R2010b, because of limited availability of operations.
sz=size(data).';data=data(:);

persistent has_int64_arith % 64-bit arithmetic was added in R2010b
if isempty(has_int64_arith),has_int64_arith=ifversion('>=','R2010b','Octave','>',0);end

c=class(data);
cast_int64_to_double=~options.cast_int64_double && has_int64_arith && c(end)=='4';
if ~cast_int64_to_double
    % As an undocumented hack, cast_int64_double can be set to false in the options without
    % changing the hash function version.
    if any(abs(double(data(:)))>2^52)
        warning_(options,'HJW:ComputeNonCryptHash:int64rounding',...
            ['int64 and uint64 will be rounded pre-R2010b, resulting in rounding.',char(10),...
            'This will result in a hash that is different from newer releases.']) %#ok<CHARTEN>
    end
end
if cast_int64_to_double
    % Ensure the data is uint64 instead of converting to double.
    if c(1)~='u'
        % Shift int64 values to the uint64 data range.
        L=data>0;x=-int64(-inf);
        data_=uint64(data+x+1);               % This clips positive data.
        data_(L)=uint64(data(L))+uint64(x)+1; % This clips negative data.
        data=data_;
    end
elseif c(1)~='u'
    % Shift int* values to the uint* data range.
    data=double(data)-double(eval([c '(-inf)']));
else
    data=double(data);
end
switch c(end)
    case '8'
        % Append a 0 for odd length and merge pairs.
        if mod(numel(data),2),data(end+1)=0;end
        data=reshape(data,[],2);
        data=data(:,1)*255+data(:,2);
        data=uint16(data);
    case '6'
        data=uint16(data);
    case '2'
        % Split to 2 words.
        bit1=floor(data/2^16);bit1=bit1.';
        bit2=mod(data,2^16);  bit2=bit2.';
        data=[bit1;bit2];
        data=uint16(data(:));
    case '4'
        % Split to 4 words. Note that bit4 contains a rounding error for data >2^52 pre-R2010b.
        bit64=data;
        bit4=mod(bit64,2^16);bit64=bit64-bit4;bit64=bit64/2^16;bit4=bit4.';
        bit3=mod(bit64,2^16);bit64=bit64-bit3;bit64=bit64/2^16;bit3=bit3.';
        bit2=mod(bit64,2^16);bit64=bit64-bit2;bit64=bit64/2^16;bit2=bit2.';
        bit1=mod(bit64,2^16);                                  bit1=bit1.';
        
        data=[bit1;bit2;bit3;bit4];
        data=uint16(data(:));
end
% Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz,2^16))];
end
function data=cast_to_uint16_vector__char(data,options)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if isOctave && options.re_encode_char
    isColVector = size(data,1)==numel(data);
    if isColVector,data=data.';end
    data=cellstr(data);
    % Decode from UTF-8 and encode with UTF-16 (the output will be uint16).
    for n=1:numel(data)
        data{n}=unicode_to_char(UTF8_to_unicode(data{n},options.print_to),true);
    end
    % Pad with spaces (but in uint16).
    cell_length=cellfun('length',data);longest_cell=max(cell_length);
    for n=find(cell_length<longest_cell)
        data{n}( (numel(data{n})+1) : longest_cell)=uint16(' ');
    end
    data=cell2mat(data);
    if isColVector,data=data.';end
end
sz=size(data).';data=data(:);
data=uint16(data); % Chars are 16 bit in Matlab, as they are encoded with UTF-16.
data=[data;uint16(mod(sz,2^16))]; % Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__string(data,options)
if options.string_to_cellstr
    data=cellstr(data);
    data=cast_to_uint16_vector__cell(data,options);
else
    data=char(data); % Cast to char instead of a cell array of chars.
    data=cast_to_uint16_vector__char(data,options);
end
end
function data=cast_to_uint16_vector__struct(data,options)
sz=size(data).';data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output,options);
data=[data;uint16(mod(sz,2^16))]; % Append the array size to make it influence the hash.
end
function c=char2cellstr(str,LineEnding)
%Split char or uint32 vector to cell (1 cell element per line). Default splits are for CRLF/CR/LF.
%The input data type is preserved.
%
%Since the largest valid Unicode codepoint is 0x10FFFF (i.e. 21 bits), all values will fit in an
%int32 as well. This is used internally to deal with different newline conventions.
%
%The second input is a cellstr containing patterns that will be considered as newline encodings.
%This will not be checked for any overlap and will be processed sequentially.

returnChar=isa(str,'char');
str=int32(str);%convert to signed, this should not crop any valid Unicode codepoints.

if nargin<2
    %Replace CRLF, CR, and LF with -10 (in that order). That makes sure that all valid encodings of
    %newlines are replaced with the same value. This should even handle most cases of files that
    %mix the different styles, even though such mixing should never occur in a properly encoded
    %file. This considers LFCR as two line endings.
    if any(str==13)
        str=PatternReplace(str,int32([13 10]),int32(-10));
        str(str==13)=-10;
    end
    str(str==10)=-10;
else
    for n=1:numel(LineEnding)
        str=PatternReplace(str,int32(LineEnding{n}),int32(-10));
    end
end

%Split over newlines.
newlineidx=[0 find(str==-10) numel(str)+1];
c=cell(numel(newlineidx)-1,1);
for n=1:numel(c)
    s1=(newlineidx(n  )+1);
    s2=(newlineidx(n+1)-1);
    c{n}=str(s1:s2);
end

%Return to the original data type.
if returnChar
    for n=1:numel(c),c{n}=  char(c{n});end
else
    for n=1:numel(c),c{n}=uint32(c{n});end
end
end
function tf=CharIsUTF8
%This provides a single place to determine if the runtime uses UTF-8 or UTF-16 to encode chars.
%The advantage is that there is only 1 function that needs to change if and when Octave switches to
%UTF-16. This is unlikely, but not impossible.
persistent persistent_tf
if isempty(persistent_tf)
    if exist('OCTAVE_VERSION','builtin')~=0 %Octave
        %Test if Octave has switched to UTF-16 by looking if the Euro symbol is losslessly encoded
        %with char.
        w=struct('w',warning('off','all'));[w.msg,w.ID]=lastwarn;
        persistent_tf=~isequal(8364,double(char(8364)));
        warning(w.w);lastwarn(w.msg,w.ID); % Reset warning state.
    else
        persistent_tf=false;
    end
end
tf=persistent_tf;
end
function opts=ComputeNonCryptHash_DefaultsByVersion(opts)
%These defaults are based on the version number.

opts.Version=str2double(opts.VersionFlag(3:end));

if ~isfield(opts,'re_encode_char_on_Octave')
    %Chars will be decoded from UTF-8 and encoded with UTF-16.
    opts.re_encode_char=opts.Version>=2;
end

if ~isfield(opts,'string_to_cellstr')
    %Cast strings to a cell of chars, instead of a char.
    opts.string_to_cellstr=opts.Version>=2;
    
    if ~isfield(opts,'cast_int64_double')
        %Always cast int64 and uint64 (results in rounding >2^52).
        opts.cast_int64_double=opts.Version==1;
    end
end
end
function [success,options,ME]=ComputeNonCryptHash_parse_inputs(varargin)
%Parse the inputs to a struct.
%
%  hash=ComputeNonCryptHash(data)
%  hash=ComputeNonCryptHash(___,HashLength)
%  hash=ComputeNonCryptHash(___,VersionFlag)
%  hash=ComputeNonCryptHash(___,options)
%  hash=ComputeNonCryptHash(___,Name,Value)

%Assign default outputs.
success=false;
options=ComputeNonCryptHash_parse_inputs_default;

[ArgIn,ME]=ParseArgs(varargin{:});
if ~isempty(ME),options=default;return;end

% Attempt to match the struct to the available options. This will return a struct with the same
% fields as the default option struct.
try ME_=[];[options,replaced]=parse_NameValue(options,ArgIn{:}); %#ok<NASGU>
catch ME_;if isempty(ME_),ME_=lasterror;end,ME=ME_;return,end %#ok<LERR>

% Attempt to parse the error redirection options (this generates an ME struct on fail).
[options.print_to,ME]=validate_print_to__options(options,ME);
if isempty(options.print_to),return,end

if numel(replaced)==0,success=true;ME=[];return,end % No default values were changed.

%Perform the actual input validation.
for k=1:numel(replaced)
    if strcmp(replaced{k},'VersionFlag')
        %Test VersionFlag and set Version
        try
            Version=str2double(options.VersionFlag(3:end));
            if isnan(Version) || round(Version)~=Version || Version>2
                error('trigger');
            end
        catch
            ME.identifier='HJW:ComputeNonCryptHash:InputIncorrect';
            ME.message='VersionFlag input incorrect. Must be ''-v1'', ''-v2''.';
            return
        end
    elseif strcmp(replaced{k},'HashLength')
        %Test HashLength
        HashLength=options.HashLength;
        if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
            ME.identifier='HJW:ComputeNonCryptHash:InputIncorrect';
            ME.message='HashLength input must be a multiple of 16.';
            return
        end
    end
end

success=true;ME=[];
end

function [ArgIn,ME]=ParseArgs(varargin)
ME=[];
ArgIn=varargin;
extracted=[false false];
for count=1:nargin
    % Pop the two optional parameters from the varargin.
    [flag,arg,ArgIn]=ComputeNonCryptHash_parse_inputs__pop_queue(ArgIn);
    switch flag
        case 0,break
        case 1,VersionFlag=arg; extracted(1)=true;
        case 2,HashLength=arg;  extracted(2)=true;
        case 3
            ME.identifier='HJW:ComputeNonCryptHash:InputIncorrect';
            ME.message='Unable to determine a valid syntax.';
            return
    end
end

if nargin==sum(extracted)
    % Only VersionFlag and HashLength were present in the argument list.
    options=struct;
    if extracted(1),options.VersionFlag=VersionFlag;end
    if extracted(2),options.HashLength =HashLength ;end
    ArgIn={options};
else
    % VersionFlag and HashLength are extracted from the argument list, so now were left with the
    % normal optionstruct/NameValue in ArgIn.
    if isa(ArgIn{1},'stuct')
        % Add to the struct.
        options=ArgIn{1};
        if extracted(1),options.VersionFlag=VersionFlag;end
        if extracted(2),options.HashLength =HashLength ;end
        ArgIn{1}=options;
    else
        % Add to the Name,Value list.
        x=cell(1,0);
        if extracted(1),x=[x {'VersionFlag',VersionFlag}];end
        if extracted(2),x=[x {'HashLength' ,HashLength }];end
        ArgIn=[x ArgIn];
    end
end
end

function [flag,arg,argout]=ComputeNonCryptHash_parse_inputs__pop_queue(argin)
% The output flag encodes a VersionFlag (1), a HashLength (2), or an invalid input (3).
% If the output flag is 0, a Name,Value or struct was found.
arg=argin{1};
if ismember(class(arg),{'char','string'})
    % Is the parameter VersionFlag or the first from a NameValue pair?
    if isa(arg,'char')&&numel(arg)>=1 && strcmp('-',arg(1))
        flag=1;
    else
        flag=0;
    end
else
    if isa(arg,'struct')
        flag=0;
    else
        flag=2;
    end
end
argout=argin;if flag~=0,argout(1)=[];end
end
function opts=ComputeNonCryptHash_parse_inputs_default
persistent opts_
if isempty(opts_)
    opts_=struct;
    opts_.HashLength=256;
    opts_.VersionFlag='-v2';
    opts_.SkipInputParse=false;%undocumented shortcut
    opts_.isSaltCall=false;%undocumented: for use in internal call
    
    %Set defaults for the error redirection.
    opts_.print_to_con=[];
    opts_.print_to_fid=[];
    opts_.print_to_obj=[];
    opts_.print_to_fcn=[];
    opts_.print_to=validate_print_to__options(opts_);
end
opts=opts_;
end
function error_(options,varargin)
%Print an error to the command window, a file and/or the String property of an object.
%The error will first be written to the file and object before being actually thrown.
%
%Apart from controlling the way an error is written, you can also run a specific function. The
%'fcn' field of the options must be a struct (scalar or array) with two fields: 'h' with a function
%handle, and 'data' with arbitrary data passed as third input. These functions will be run with
%'error' as first input. The second input is a struct with identifier, message, and stack as
%fields. This function will be run with feval (meaning the function handles can be replaced with
%inline functions or anonymous functions).
%
%The intention is to allow replacement of every error(___) call with error_(options,___).
%
% NB: the error trace that is written to a file or object may differ from the trace displayed by
% calling the builtin error function. This was only observed when evaluating code sections.
%
%options.boolean.con: if true throw error with rethrow()
%options.fid:         file identifier for fprintf (array input will be indexed)
%options.boolean.fid: if true print error to file
%options.obj:         handle to object with String property (array input will be indexed)
%options.boolean.obj: if true print error to object (options.obj)
%options.fcn          struct (array input will be indexed)
%options.fcn.h:       handle of function to be run
%options.fcn.data:    data passed as third input to function to be run (optional)
%options.boolean.fnc: if true the function(s) will be run
%
%syntax:
%  error_(options,msg)
%  error_(options,msg,A1,...,An)
%  error_(options,id,msg)
%  error_(options,id,msg,A1,...,An)
%  error_(options,ME)               %equivalent to rethrow(ME)
%
%examples options struct:
%  % Write to a log file:
%  opts=struct;opts.fid=fopen('log.txt','wt');
%  % Display to a status window and bypass the command window:
%  opts=struct;opts.boolean.con=false;opts.obj=uicontrol_object_handle;
%  % Write to 2 log files:
%  opts=struct;opts.fid=[fopen('log2.txt','wt') fopen('log.txt','wt')];

persistent this_fun
if isempty(this_fun),this_fun=func2str(@error_);end

%Parse options struct, allowing an empty input to revert to default.
if isempty(options),options=validate_print_to__options(struct);end
options                   =parse_warning_error_redirect_options(  options  );
[id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs( varargin{:});
if no_op,return,end
ME=struct('identifier',id,'message',msg,'stack',stack);

%Print to object.
if options.boolean.obj
    msg_=msg;while msg_(end)==10,msg_(end)='';end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend 'Error: '.
        msg_=char2cellstr(['Error: ' msg_]);
    else              % Only prepend 'Error: '.
        msg_=['Error: ' msg_];
    end
    for OBJ=options.obj(:).'
        try set(OBJ,'String',msg_);catch,end
    end
end

%Print to file.
if options.boolean.fid
    T=datestr(now,31); % Print the time of the error to the log as well.
    for FID=options.fid(:).'
        try fprintf(FID,'[%s] Error: %s\n%s',T,msg,trace);catch,end
    end
end

%Run function.
if options.boolean.fcn
    if ismember(this_fun,{stack.name})
        %To prevent an infinite loop, trigger an error.
        error('prevent recursion')
    end
    for FCN=options.fcn(:).'
        if isfield(FCN,'data')
            try feval(FCN.h,'error',ME,FCN.data);catch,end
        else
            try feval(FCN.h,'error',ME);catch,end
        end
    end
end

%Actually throw the error.
rethrow(ME)
end
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers=1;end
if nargin<2, stack=dbstack;end
stack(1:skip_layers)=[];

%Parse the ML6.5 style of dbstack (the name field includes full file location).
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp=stack(n).name;
        if strcmp(tmp(end),')')
            %Internal function.
            ind=strfind(tmp,'(');
            name=tmp( (ind(end)+1):(end-1) );
            file=tmp(1:(ind(end)-2));
        else
            file=tmp;
            [ignore,name]=fileparts(tmp); %#ok<ASGLU>
        end
        [ignore,stack(n).file]=fileparts(file); %#ok<ASGLU>
        stack(n).name=name;
    end
end

%Parse Octave style of dbstack (the file field includes full file location).
persistent IsOctave,if isempty(IsOctave),IsOctave=exist('OCTAVE_VERSION','builtin');end
if IsOctave
    for n=1:numel(stack)
        [ignore,stack(n).file]=fileparts(stack(n).file); %#ok<ASGLU>
    end
end

%Create the char array with a (potentially) modified stack.
s=stack;
c1='>';
str=cell(1,numel(s)-1);
for n=1:numel(s)
    [ignore_path,s(n).file,ignore_ext]=fileparts(s(n).file); %#ok<ASGLU>
    if n==numel(s),s(n).file='';end
    if strcmp(s(n).file,s(n).name),s(n).file='';end
    if ~isempty(s(n).file),s(n).file=[s(n).file '>'];end
    str{n}=sprintf('%c In %s%s (line %d)\n',c1,s(n).file,s(n).name,s(n).line);
    c1=' ';
end
str=horzcat(str{:});
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
%   tf=ifversion(test,Rxxxxab)
%   tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Input/output arguments:
% tf:
%   If the current version satisfies the test this returns true. This works similar to verLessThan.
% Rxxxxab:
%   A char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the numeric
%   version (e.g. 6.5, 7, or 9.6).
% test:
%   A char array containing a logical test. The interpretation of this is equivalent to
%   eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.1.0                                                         |%
%|  Date:    2021-09-26                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910;'R2021b' 911};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
        else
            L=ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf=NaN;return
            else
                v=v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, ignoring incomplete names and case
%
% The first output is a struct with the same fields as the first input, with field contents
% replaced according to the supplied options struct or Name,Value pairs.
% The second output is a cell containing the field names that have been set.
%
% If this fails to find a match, this will throw an error with the offending name as the message.
%
% If there are multiple occurences of a Name, only the last Value will be returned. This is the
% same as Matlab internal functions like plot. GNU Octave also has this behavior.
%
% If a struct array is provided, only the first element will be used. An empty struct array will
% trigger an error.

opts=default;replaced={};
if nargin==1,return,end

%Unwind an input struct to Name,Value pairs.
try
    struct_input=numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input=mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        % If this is a struct array, explicitly select the first element.
        Names=fieldnames(varargin{1}(1));
        Values=struct2cell(varargin{1}(1));
    else
        %Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names=cellstr(varargin(1:2:end));
        Values=varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    %If this block errors, that is either because a missing Value with the Name,Value syntax, or
    %because the struct input is not a struct, or because an attempt was made to mix the two
    %styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as struct.')
end

%Convert the real fieldnames to a char matrix (with and without removing underscores).
d_Names1=fieldnames(default);
d_Names2=lower(d_Names1);
len=cellfun('prodofsize',d_Names2);maxlen=max(len);
for n=find(len<maxlen).' %Pad with spaces where needed
    d_Names2{n}((end+1):maxlen)=' ';
end
d_Names2=vertcat(d_Names2{:});
if any(any(d_Names2=='_'))
    d_Names3=strrep(lower(d_Names1),'_','');
    len=cellfun('prodofsize',d_Names3);maxlen=max(len);
    for n=find(len<maxlen).' %Pad with spaces where needed
        d_Names3{n}((end+1):maxlen)=' ';
    end
    d_Names3=vertcat(d_Names3{:});
else
    d_Names3=d_Names2;
end

%Attempt to match the names.
replaced=false(size(d_Names1));
for n=1:numel(Names)
    name=lower(Names{n});
    
    tmp=d_Names2(:,1:min(end,numel(name)));
    non_matching=numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
    match_idx=find(non_matching==0);
    if numel(match_idx)~=1
        % Retry without underscores.
        tmp=d_Names3(:,1:min(end,numel(name)));
        non_matching=numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
        match_idx=find(non_matching==0);
        
        if numel(match_idx)~=1
            error('parse_NameValue:NonUniqueMatch',Names{n})
        end
    end
    
    %Store the Value in the output struct.
    opts.(d_Names1{match_idx})=Values{n};
    
    replaced(match_idx)=true;
end
replaced=d_Names1(replaced);
end
function [id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs(varargin)
no_op=false;
if nargin==1
    %  error_(options,msg)
    %  error_(options,ME)
    if isa(varargin{1},'struct') || isa(varargin{1},'MException')
        ME=varargin{1};
        if numel(ME)==0
            no_op=true;
            [id,msg,stack,trace]=deal('');
            return
        end
        try
            stack=ME.stack;%Use the original call stack if possible.
            trace=get_trace(0,stack);
        catch
            [trace,stack]=get_trace(3);
        end
        id=ME.identifier;
        msg=ME.message;
        pat='Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback(';
        %This pattern may occur when using try error(id,msg),catch,ME=lasterror;end instead of
        %catching the MException with try error(id,msg),catch ME,end.
        %This behavior is not stable enough to robustly check for it, but it only occurs with
        %lasterror, so we can use that.
        if isa(ME,'struct') && numel(msg)>numel(pat) && strcmp(pat,msg(1:numel(pat)))
            %Strip the first line (which states 'error in function (line)', instead of only msg).
            msg(1:find(msg==10,1))='';
        end
    else
        [trace,stack]=get_trace(3);
        [id,msg]=deal('',varargin{1});
    end
else
    [trace,stack]=get_trace(3);
    if ~isempty(strfind(varargin{1},'%')) %The id can't contain a percent symbol.
        %  error_(options,msg,A1,...,An)
        id='';
        A1_An=varargin(2:end);
        msg=sprintf(varargin{1},A1_An{:});
    else
        %  error_(options,id,msg)
        %  error_(options,id,msg,A1,...,An)
        id=varargin{1};
        msg=varargin{2};
        if nargin>2
            A1_An=varargin(3:end);
            msg=sprintf(msg,A1_An{:});
        end
    end
end
end
function options=parse_warning_error_redirect_options(options)
%Fill the struct:
%options.boolean.con (this field is ignored in error_)
%options.boolean.fid
%options.boolean.obj
%options.boolean.fcn
if ~isfield(options,'boolean'),options.boolean=struct;end
if ~isfield(options.boolean,'con') || isempty(options.boolean.con)
    options.boolean.con=false;
end
if ~isfield(options.boolean,'fid') || isempty(options.boolean.fid)
    options.boolean.fid=isfield(options,'fid');
end
if ~isfield(options.boolean,'obj') || isempty(options.boolean.obj)
    options.boolean.obj=isfield(options,'obj');
end
if ~isfield(options.boolean,'fcn') || isempty(options.boolean.fcn)
    options.boolean.fcn=isfield(options,'fcn');
end
end
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
out=in(:)';
if numel(pattern)==0
    L=false(size(in));
elseif numel(rep)>numel(pattern)
    error('not implemented (padding required)')
else
    L=true(size(in));
    for n=1:numel(pattern)
        k=find(in==pattern(n));
        k=k-n+1;k(k<1)=[];
        %Now k contains the indices of the beginning of each match.
        L2=false(size(L));L2(k)=true;
        L= L & L2;
        if ~any(L),break,end
    end
end
k=find(L);
if ~isempty(k)
    for n=1:numel(rep)
        out(k+n-1)=rep(n);
    end
    if numel(rep)==0,n=0;end
    if numel(pattern)>n
        k=k(:);%Enforce direction.
        remove=(n+1):numel(pattern);
        idx=bsxfun_plus(k,remove-1);
        out(idx(:))=[];
    end
end
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%The char and string test are not case sensitive.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
%- 'enable' or 'disable'
%- 'enabled' or 'disabled'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    try
        states(end+1,:)=eval('{"on","off"}'); %#ok<EVLCS>
    catch
    end
end
isLogical=true;
try
    if isa(val,'char') || isa(val,'string')
        try val=lower(val);catch,end
    end
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end
function [bit64,last16bits]=typecast_double_uint64(FP)
%Turn a double into a uint64 with the same binary representation.
%This is similar to typecast(FP,'uint64'); the difference being that this function ignores
%endianness, supports array inputs, and is slower.
%
%Because of missing support for some operations in ML6.5, this function returns the uint64 as a
%double. Because this may cause rounding errors, the last 16 bits are returned separately.

[M,E]=log2(FP);
signBit =-floor(sign(FP)/2-0.5);
exponent=E+1022;
mantissa=abs(M)*2-1;

%no plus() for integer types in ML6.5, so we need to use double, instead of uint64
bit64=zeros(size(FP));
bit64=bit64+(signBit*2^63);
bit64=bit64+(exponent*2^52);
bit64=bit64+(mantissa*2^52);
last16bits=mod(mantissa*2^52,2^16);

%correct the 0 and hard-code the special cases
L=isinf(FP);
bit64(FP==0)=0;
bit64(isnan(FP))=18444492273895866368;
bit64(L & FP>0)=9218868437227405312;%positive inf
bit64(L & FP<0)=18442240474082181120;%negative inf
last16bits(FP==0)=0;
last16bits(isnan(FP))=0;
last16bits(L)=0;
end
function str=unicode_to_char(unicode,encode_as_UTF16)
%Encode Unicode code points with UTF-16 on Matlab and UTF-8 on Octave.
%Input is either implicitly or explicitly converted to a row-vector.

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if nargin==1
    encode_as_UTF16=~CharIsUTF8;
end
if encode_as_UTF16
    if all(unicode<65536)
        str=uint16(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-16.
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF16(char_list(n));
            str_element=uint16(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
    end
    if ~isOctave
        str=char(str);%Conversion to char could trigger a conversion range error in Octave.
    end
else
    if all(unicode<128)
        str=char(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-8
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF8(char_list(n));
            str_element=uint8(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
        str=char(str);
    end
end
end
function str=unicode_to_UTF16(unicode)
%Convert a single character to UTF-16 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
if unicode<65536
    str=unicode;return
end
U=double(unicode)-65536;%Convert to double for ML6.5.
U=dec2bin(U,20);
str=bin2dec(['110110' U(1:10);'110111' U(11:20)]).';
end
function str=unicode_to_UTF8(unicode)
%Convert a single character to UTF-8 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-8
if numel(unicode)>1,error('this should only be used for single characters'),end
if unicode<128
    str=unicode;return
end
persistent pers
if isempty(pers)
    pers=struct;
    pers.limits.lower=hex2dec({'0000','0080','0800', '10000'});
    pers.limits.upper=hex2dec({'007F','07FF','FFFF','10FFFF'});
    pers.scheme{2}='110xxxxx10xxxxxx';
    pers.scheme{2}=reshape(pers.scheme{2}.',8,2);
    pers.scheme{3}='1110xxxx10xxxxxx10xxxxxx';
    pers.scheme{3}=reshape(pers.scheme{3}.',8,3);
    pers.scheme{4}='11110xxx10xxxxxx10xxxxxx10xxxxxx';
    pers.scheme{4}=reshape(pers.scheme{4}.',8,4);
    for b=2:4
        pers.scheme_pos{b}=find(pers.scheme{b}=='x');
        pers.bits(b)=numel(pers.scheme_pos{b});
    end
end
bytes=find(pers.limits.lower<=unicode & unicode<=pers.limits.upper);
str=pers.scheme{bytes};
scheme_pos=pers.scheme_pos{bytes};
b=dec2bin(unicode,pers.bits(bytes));
str(scheme_pos)=b;
str=bin2dec(str.').';
end
function [unicode,isUTF8,assumed_UTF8]=UTF8_to_unicode(UTF8,print_to)
%Convert UTF-8 to the code points stored as uint32
% Plane 16 goes up to 10FFFF, so anything larger than uint16 will be able to hold every code point.
%
% If there a second output argument, this function will not return an error if there are encoding
% error. The second output will contain the attempted conversion, while the first output will
% contain the original input converted to uint32.
%
% The second input can be used to also print the error to a GUI element or to a text file.
if nargin<2,print_to=[];end
return_on_error= nargout==1 ;

UTF8=uint32(reshape(UTF8,1,[]));% Force row vector.
[assumed_UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error);
if strcmp(flag,'success')
    isUTF8=true;
    unicode=assumed_UTF8;
elseif strcmp(flag,'error')
    isUTF8=false;
    if return_on_error
        error_(print_to,ME)
    end
    unicode=UTF8; % Return input unchanged (apart from casting to uint32).
end
end
function [UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error)
flag='success';
ME=struct('identifier','HJW:UTF8_to_unicode:notUTF8','message','Input is not UTF-8.');

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end

if any(UTF8>255)
    flag='error';
    if return_on_error,return,end
elseif all(UTF8<128)
    return
end

for bytes=4:-1:2
    val=bin2dec([repmat('1',1,bytes) repmat('0',1,8-bytes)]);
    multibyte=UTF8>=val & UTF8<256;%Exclude the already converted chars.
    if any(multibyte)
        multibyte=find(multibyte);multibyte=multibyte(:).';
        if numel(UTF8)<(max(multibyte)+bytes-1)
            flag='error';
            if return_on_error,return,end
            multibyte( (multibyte+bytes-1)>numel(UTF8) )=[];
        end
        if ~isempty(multibyte)
            idx=bsxfun_plus(multibyte , (0:(bytes-1)).' );
            idx=idx.';
            multibyte=UTF8(idx);
        end
    else
        multibyte=[];
    end
    header_bits=[repmat('1',1,bytes-1) repmat('10',1,bytes)];
    header_locs=unique([1:(bytes+1) 1:8:(8*bytes) 2:8:(8*bytes)]);
    if numel(multibyte)>0
        multibyte=unique(multibyte,'rows');
        S2=mat2cell(multibyte,ones(size(multibyte,1),1),bytes);
        for n=1:numel(S2)
            bin=dec2bin(double(S2{n}))';
            % To view the binary data, you can use this: bin=bin(:)';
            % Remove binary header (3 byte example):
            % 1110xxxx10xxxxxx10xxxxxx
            %     xxxx  xxxxxx  xxxxxx
            if ~strcmp(header_bits,bin(header_locs))
                % Check if the byte headers match the UTF-8 standard.
                flag='error';
                if return_on_error,return,end
                continue %leave unencoded
            end
            bin(header_locs)='';
            if ~isOctave
                S3=uint32(bin2dec(bin  ));
            else
                S3=uint32(bin2dec(bin.'));% Octave needs an extra transpose.
            end
            % Perform actual replacement.
            UTF8=PatternReplace(UTF8,S2{n},S3);
        end
    end
end
end
function [opts,ME]=validate_print_to__options(opts_in,ME)
% If any input is invalid, this returns an empty array and sets ME.message and ME.identifier.
%
% Input struct:
% options.print_to_con=true;   % or false
% options.print_to_fid=fid;    % or []
% options.print_to_obj=h_obj;  % or []
% options.print_to_fcn=struct; % or []
%
% Output struct:
% options.fid
% options.obj
% options.fcn.h
% options.fcn.data
% options.boolean.con
% options.boolean.fid
% options.boolean.obj
% options.boolean.fcn

% Set defaults.
if nargin<2,ME=struct;end
if ~isfield(opts_in,'print_to_con'),opts_in.print_to_con=[];end
if ~isfield(opts_in,'print_to_fid'),opts_in.print_to_fid=[];end
if ~isfield(opts_in,'print_to_obj'),opts_in.print_to_obj=[];end
if ~isfield(opts_in,'print_to_fcn'),opts_in.print_to_fcn=[];end
print_to_con_default=true; % Unless a valid fid, obj, or fcn is specified.

% Initialize output.
opts=struct;

% Parse the fid. We can use ftell to determine if fprintf is going to fail.
item=opts_in.print_to_fid;
if isempty(item)
    opts.boolean.fid=false;
    opts.fid=[];
else
    print_to_con_default=false;
    opts.boolean.fid=true;
    opts.fid=item;
    for n=1:numel(item)
        try position=ftell(item(n));catch,position=-1;end
        if item(n)~=1 && position==-1
            ME.message=['Invalid print_to_fid parameter:',char(10),...
                'should be a valid file identifier or 1.']; %#ok<CHARTEN>
            ME.identifier='HJW:print_to:ValidationFailed';
            opts=[];return
        end
    end
end

% Parse the object handle. Retrieving from multiple objects at once works, but writing that output
% back to multiple objects doesn't work if Strings are dissimilar.
item=opts_in.print_to_obj;
if isempty(item)
    opts.boolean.obj=false;
    opts.obj=[];
else
    print_to_con_default=false;
    opts.boolean.obj=true;
    opts.obj=item;
    for n=1:numel(item)
        try
            txt=get(item(n),'String'    ); % See if this triggers an error.
            set(    item(n),'String','' ); % Test if property is writable.
            set(    item(n),'String',txt); % Restore original content.
        catch
            ME.message=['Invalid print_to_obj parameter:',char(10),...
                'should be a handle to an object with a writeable String property.']; %#ok<CHARTEN>
            ME.identifier='HJW:print_to:ValidationFailed';
            opts=[];return
        end
    end
end

% Parse the function handles.
item=opts_in.print_to_fcn;
if isempty(item)
    opts.boolean.fcn=false;
    opts.fcn=[];
else
    print_to_con_default=false;
    try
        for n=1:numel(item)
            if ~ismember(class(item(n).h),{'function_handle','inline'}) || numel(item(n).h)~=1
                error('trigger error')
            end
        end
    catch
        ME.message=['Invalid print_to_fcn parameter:',char(10),...
            'should be a struct with the h field containing a function handle,',char(10),...
            'anonymous function or inline function.']; %#ok<CHARTEN>
        ME.identifier='HJW:print_to:ValidationFailed';
        opts=[];return
    end
end

% Parse the logical that determines if a warning will be printed to the command window.
% This is true by default, unless an fid, obj, or fcn is specified.
item=opts_in.print_to_con;
if isempty(item)
    opts.boolean.con=print_to_con_default;
else
    [passed,opts.boolean.con]=test_if_scalar_logical(item);
    if ~passed
        ME.message=['Invalid print_to_con parameter:',char(10),...
            'should be a scalar logical.']; %#ok<CHARTEN>
        ME.identifier='HJW:print_to:ValidationFailed';
        opts=[];return
    end
end
end
function warning_(options,varargin)
%Print a warning to the command window, a file and/or the String property of an object
% The lastwarn state will be set if the warning isn't thrown with warning().
% The printed call trace omits this function, but the warning() call does not.
%
% You can also provide a struct (scalar or array) with two fields: 'h' with a function handle, and
% 'data' with arbitrary data passed as third input. These functions will be run with 'warning' as
% first input. The second input is a struct with identifier, message, and stack as fields. This
% function will be run with feval (meaning the function handles can be replaced with inline
% functions or anonymous functions).
%
% The intention is to allow replacement of most warning(___) call with warning_(options,___). This
% does not apply to calls that query or set the warning state.
%
%options.boolean.con: if true print warning to command window with warning()
%options.fid:         file identifier for fprintf (array input will be indexed)
%options.boolean.fid: if true print warning to file (options.fid)
%options.obj:         handle to object with String property (array input will be indexed)
%options.boolean.obj: if true print warning to object (options.obj)
%options.fcn          struct (array input will be indexed)
%options.fcn.h:       handle of function to be run
%options.fcn.data:    data passed as third input to function to be run (optional)
%options.boolean.fnc: if true the function(s) will be run
%
% Syntax:
%   warning_(options,msg)
%   warning_(options,msg,A1,...,An)
%   warning_(options,id,msg)
%   warning_(options,id,msg,A1,...,An)
%   warning_(options,ME)               %rethrow error as warning
%
%examples options struct:
%  % Write to a log file:
%  opts=struct;opts.fid=fopen('log.txt','wt');
%  % Display to a status window and bypass the command window:
%  opts=struct;opts.boolean.con=false;opts.obj=uicontrol_object_handle;
%  % Write to 2 log files:
%  opts=struct;opts.fid=[fopen('log2.txt','wt') fopen('log.txt','wt')];

persistent this_fun
if isempty(this_fun),this_fun=func2str(@warning_);end

[id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs( varargin{:});
if no_op,return,end % Don't waste time parsing the options in case of a no-op.
ME=struct('identifier',id,'message',msg,'stack',stack);

% Parse options struct, allowing an empty input to revert to default.
if isempty(options),options=validate_print_to__options(struct);end
options                   =parse_warning_error_redirect_options(  options  );

% Check if the warning is turned off and exit the function if this is the case.
w=warning;if any(ismember({w(ismember({w.identifier},{id,'all'})).state},'off')),return,end

if options.boolean.con
    if ~isempty(id),warning(id,'%s',msg),else,warning(msg), end
else
    if ~isempty(id),lastwarn(msg,id);    else,lastwarn(msg),end
end

if options.boolean.obj
    msg_=msg;while msg_(end)==10,msg_(end)=[];end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend warning.
        msg_=char2cellstr(['Warning: ' msg_]);
    else              % Only prepend warning.
        msg_=['Warning: ' msg_];
    end
    set(options.obj,'String',msg_)
    for OBJ=options.obj(:).'
        try set(OBJ,'String',msg_);catch,end
    end
end

if options.boolean.fid
    T=datestr(now,31); % Print the time of the warning to the log as well.
    for FID=options.fid(:).'
        try fprintf(FID,'[%s] Warning: %s\n%s',T,msg,trace);catch,end
    end
end

if options.boolean.fcn
    if ismember(this_fun,{stack.name})
        % To prevent an infinite loop, trigger an error.
        error('prevent recursion')
    end
    for FCN=options.fcn(:).'
        if isfield(FCN,'data')
            try feval(FCN.h,'warning',ME,FCN.data);catch,end
        else
            try feval(FCN.h,'warning',ME);catch,end
        end
    end
end
end