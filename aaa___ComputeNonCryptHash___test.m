function aaa___ComputeNonCryptHash___test(varargin)
%test if the syntax checks work and if the output is stable

% suppress v6.5 warning
v=version;v(strfind(v,'.'):end)='';v=str2double(v);
if v<=7
    warning off MATLAB:m_warning_end_without_block;clc
end

warning('create test for memory error in cast_to_uint16_vector')

%syntax tests
try test_syntax;catch,error('syntax test failed');end,clc


[input,hash]=get_test_cases;
for k=1:size(hash{1},2)
    HashLength=numel(hash{1}{1,k})*4;
    for n=1:numel(input)
        for v=1:2
            H=ComputeNonCryptHash(input{n},HashLength,sprintf('-v%d',v));
            if ~strcmp(hash{v}{n,k},H)
                error('hash did not match: k=%d,n=%d (v=%d)\n',k,n,v)
            end
        end
    end
end
clc,disp('test finished successfully')
end
function [cases,hash]=get_test_cases
%clc,for n=1:10,fprintf('    ''%s'',...\n',ComputeNonCryptHash(rand,256,'-v1'));end
str=['6A9A0BB37950EFC198B354CBFB2191E7095C3309A4D233163E62801CDF5A2976',...
    '1124707FC61D818D0DEC78B6784F37AF57144039C2AC5007DBFBA7E6A389D54F',...
    '23BAACA122C3CA3C35D4156FE0DB8BEE8BAE40F2A1B4C1AE37B17CA1CCD08A13',...
    '88FBB55BDA4CE3D6B270EAE270637BA522D5E0F2330BF47EE72A3C406E16BF08',...
    '9B2C731EFB8F2F2D2FEA27E1DD29E52C72803011C73C9285ED076112A4A2B042',...
    '866C8A112B6454F9B0FAF5C953E46744E149772931881F304BBE1EB61A972EA5',...
    'D63997E8BC1596FFF592935A37157DD0AA1D4A0BD0AC7FD8DD9C21A7A48C02BC',...
    'B95E432987FA74501570F808D010DFC9E86B885A1AED6CB65F0D9A146BB44874',...
    '3C670244B256AFE6AAC40178261B71422D5BC43A56F5F285F2D00B8EF0F6D2BC',...
    'D73C5364EE0D722DFE1A66FC777068D0346B8F8CF1A13CA8DA013657AA2C7E27'];
cases{1}=hex2im(str);
cases{2}=uint8(1);
cases{3}=int8(-50);%repeat this value to confirm unique hashes caused by data type difference
cases{4}=uint16(2000);
cases{5}=int16(-50);
cases{6}=uint32(54321);
cases{7}=int32(-50);
cases{8}=uint64(inf);
cases{9}=int64(-50);
cases{10}=[str;str;str(end:-1:1)];
cases{11}=-50;
cases{12}=single(-12345678);
cases{13}=single( 12345678);
cases{14}=(-12345678);
cases{15}=( 12345678);
cases{16}={struct('a',1,'b',str),{''}};
cases{17}={struct('a',1,'b',str),{[]}};

% %compute the 16 bit hash and 256 bit hash and check for collisions
% clc
% hash=cell(numel(cases),2);
% for n=1:numel(cases)
%     hash{n,1}=ComputeNonCryptHash(cases{n},16);
%     hash{n,2}=ComputeNonCryptHash(cases{n},256);
%     fprintf('    ''%s'',''%s'';\n',hash{n,1},hash{n,2})
% end
% for k=1:size(hash,2)
%     if numel(unique(hash(:,k)))~=numel(cases)
%         error('collision in test set')
%     end
% end
hash_v1={...
    '342E','C7AB19F74F03A4B65141A30699AA30C60D2D0B958B737E4F69DFBC09C9CA3293';
    'FC94','274FBFC0CF0014586EF8DA79F52FAB2AF691B933985B167ACFDE8211743003AC';
    '149D','089314ACFE0F41F05A96100746EB57DFC42F08CEB447A5C1190CC54FF98D2B41';
    'C485','B20635EA67FECC19E3B1B1C858FDCF3A1AE360F3F4BB6D319E8504642EF0CE81';
    'E362','EBDB73016028CAF6B31BED64D32FF18390A5DA6627F1C1AFC3055E4D94261A07';
    '7FB2','8BEBA9BDC403634606EC2720F9E49B7B55857CDDF91469451802802DA8644836';
    'E145','1F242D0088DB24E0266DDE69D8B685D9513CDA23D5CCBA1B88F33324667F3273';
    'F0BA','6D6ECD144DFA00CE41F55D6FE16667092FEB6B91D6211BC3A946E8943F726415';
    'FA1E','F19B1B653844968B6AE0BECF36738B3E07D7F57133F52D5DFB380AA0858FEBCF';
    'A31D','5BE85BAE14C0E30A1D58AB612D90644529CACF75B3C8E41CF80997CD9F3154EC';
    '8AF3','1F3C099580D4CA3FA755764E9D8A42D75DF0FE3439A4AF63F71A77464065C095';
    '83A5','C9136FE84F5325C99ACBD883FB140906C01D670C91A9ADAF28829F872F8F9460';
    '6D38','4F78DF4277286136D10EDC8A70ABD8E0AE3BEBAEDC86ED7CCBEB6D971C5CFB3E';
    '9D74','FC630EAE663B4F68D8C00B911D3BD1F7F11267ADFE5549288E4531239A9D08ED';
    '29DF','C2400EAE72CCA7C6A8FC546DEF0DE612A165BF92B776707357A0FD0FC051BF1A';
    '2B5E','61F63070C29CBDF584437FBD531C676D638C7D1A2CD88EA6F347205825541329';
    'B437','E2EFDB2BE3E5548125F98EABE3733E09973EBC81F81E852C8B0487224AFBD634'};
hash_v2={...
    'F466','CDC5D51DC2D7973792646DD13AE5C5A8EDB498176EEB490FC50C03C39B487023';
    '8648','FD172B80031A7B94DA8EF106C3288ED55838AD27229C2959C4263D9B11099215';
    '1896','89496CDEEE67208B8947362288811421167AA82709A9B52F89A4D40E498C5934';
    '97D2','679DDB88F642FD31B8F24694EA5F251AEE6A70F929DD8AC118A4C3039017C036';
    '8E07','42296DDBEF4351D7966B79D5D88D1EE88E92AAA7CF1B068185DBD40EF4D93040';
    'FA97','0000B11AFBB435374EF9143A4A1FC600089AAFD956A8E03952633111FB606EE9';
    'F54A','7E53A713FA36FB424DF3530B2429DD119F6D8DEACA6B8AC13462C3035F0DB839';
    '0072','00000000000000008F24A12F9054A2931D70F2354D4AC9B7ED67CEF0D83C5DE5';
    '3E7E','A0C7000000000000405C85EDAAAEDE96D554E9ACEEB0AE4D59DAE4DDD83C9E23';
    '2256','4F6C5970DCD461A86CC5E25B79952C772C9744A16F893FBA03E8169478B5A667';
    '4403','5637000000000000B664B14404228D9278BD269627D879CEAE598C3C329E3CE5';
    'B5D7','338901AA8A470000A68193CD74E87A07777547393F2D6326AA9009A4A97EE5D4';
    '7697','4EC3955211C40000A68193CD6728F1D17775473917446326AA9009A44F1B3143';
    '624B','6BB6141FDB370000E80EB144E8C86FD9AC7784ED27D8227DAA9009A4A97E346D';
    '4389','76F0A7C762B40000E80EB144A5AF8D92AC7784EDFFF0227DAA9009A44F1B7FDD';
    '6F9B','ED220A20B06AD5302D9E9DFBDC616C02A2E6E3F80D526C094EE18DF80643998F';
    '5C54','22164F2E69A215BFCBA3877670FD0325386D0C43414E584C95A9C740BFCB182C'};
hash={hash_v1,hash_v2};
end
function test_syntax
clc
ThrowError=false;
for n_test=1:10^3%basically inf
    try
        success=eval(sprintf('test%02d',n_test));
        %#ok<*DEFNU,*NASGU,*LERR>
    catch
        break
    end
    
    if success
        fprintf('test %d succeded\n',n_test)
    else
        ME=lasterror;id=ME.identifier;
        fprintf('test %d failed\n(id: %s)\n',n_test,id)
        ThrowError=true;
    end
end
if ThrowError
    error('syntax test failed')
end

end
function success=test01
success=true;
data=@sum;%unsupported data type
try
    ComputeNonCryptHash(data);
    success=false;
catch
    ME=lasterror;
    if ~strcmp(ME.identifier,'HJW:ComputeNonCryptHash:UnwindFailed')
        success=false;
    end
end
end
function success=test02
success=true;
try
    ComputeNonCryptHash;
    success=false;
catch
    ME=lasterror;
    if ~strcmp(ME.identifier,'HJW:ComputeNonCryptHash:InputIncorrect')
        success=false;
    end
end
end
function success=test03
success=true;
try
    data='1';
    ComputeNonCryptHash(data,8,'-v1');
    success=false;
catch
    ME=lasterror;
    if ~strcmp(ME.identifier,'HJW:ComputeNonCryptHash:InputIncorrect')
        success=false;
    end
end
end
function success=test04
success=true;
try
    data='1';
    ComputeNonCryptHash(data,[16 16]);
    success=false;
catch
    ME=lasterror;
    if ~strcmp(ME.identifier,'HJW:ComputeNonCryptHash:InputIncorrect')
        success=false;
    end
end
end
function success=test05
success=true;
try
    data='1';
    ComputeNonCryptHash(data,0);
    success=false;
catch
    ME=lasterror;
    if ~strcmp(ME.identifier,'HJW:ComputeNonCryptHash:InputIncorrect')
        success=false;
    end
end
end
function out=bsxfun_plus(in1,in2)
%implicit expansion for plus(), but without any input validation
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
function data=cast_to_uint16_vector(data,re_encode_char,string_to_cellstr)
%Linearize the input data and convert it to a uint16 vector
if isa(data,'uint16')
    %Append the array size and type to make it influence the hash.
    c='uint16';sz=size(data);
    data=reshape(data,[],1);%linearize
    data=[data;uint16(c.');uint16(mod(sz',2^16))];
    return
end
data=cast_to_uint16_vector__cell({data},re_encode_char,string_to_cellstr);
data([end-1 end])=[];%Remove the [1 1] that is added because of the wrapping in a cell
end
function data=cast_to_uint16_vector__cell(data,re_encode_char,string_to_cellstr)
sz=size(data);data=data(:);
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
            data{n}=cast_to_uint16_vector__integer(data{n});
        case 'char'
            data{n}=cast_to_uint16_vector__char(data{n},re_encode_char);
        case 'string'
            data{n}=cast_to_uint16_vector__string(data{n},re_encode_char,string_to_cellstr);
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n},re_encode_char,string_to_cellstr);
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n},re_encode_char,string_to_cellstr);
        otherwise
            error('HJW:cast_to_uint16_vector:nosupport',...
                'Unsupported data type in nested variable')
    end
end
data=cell2mat(data);%Merge all cell contents
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__floats(data)
sz=size(data);c=class(data);%the rest of the function treats singles as double

%Convert to a uint64, separate it into 4 words and everything merge into a vector
[bit64,bit4]=typecast_double_uint64(double(data));
bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
bit1      =mod(bit64,2^16);
data=[bit1';bit2';bit3';bit4'];
data=uint16(data(:));

%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__logical(data)
sz=size(data);data=data(:);

if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';

data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__integer(data)
%Large values (>2^52) will not have integer precision due to a conversion to double.
%This conversion is done, because of limited availability of operations on ML6.5.
sz=size(data);data=data(:);

c=class(data);
if c(1)~='u'
    %Shift int* values to the uint* data range
    data=double(data)-double(eval([c '(-inf)']));
else
    data=double(data);
end
switch c(end)
    case '8'
        %append a 0 for odd length and merge pairs
        if mod(numel(data),2),data(end+1)=0;end
        data=reshape(data,[],2);
        data=data(:,1)*255+data(:,2);
        data=uint16(data);
    case '6'
        data=uint16(data);
    case '2'
        %split to 2 words
        bit1=floor(data/2^16);
        bit2=mod(data,2^16);
        data=[bit1';bit2'];
        data=uint16(data(:));
    case '4'
        %split to 4 words, bit4 contains a rounding error for data >2^52
        bit64=data;
        bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
        bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
        bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
        bit1      =mod(bit64,2^16);
        
        data=[bit1';bit2';bit3';bit4_round'];
        data=uint16(data(:));
end
%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__char(data,re_encode_char)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if isOctave && re_encode_char
    %Only measure the size after conversion, so it can actually match between Octave and Matlab.
    isRowVector= size(data,1)==1 ;%Store the orientation.
    data=UTF8_to_unicode(data(:).');
    data=unicode_to_char(data,true);%Encode with UTF-16 (output will be uint16, not char).
    if ~isRowVector,data=data.';end
end
sz=size(data);data=data(:);
data=uint16(data);%Chars are 16 bit in Matlab, as they are encoded with UTF-16.
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__string(data,re_encode_char,string_to_cellstr)
if string_to_cellstr
    data2=cell(size(data));for n=1:numel(data2),data2{n}=char(data(n));end
    data=cast_to_uint16_vector__cell(data2,re_encode_char,string_to_cellstr);
else
    data=char(data);%Cast to char instead of a cell array of chars.
    data=cast_to_uint16_vector__char(data,re_encode_char);
end
end
function data=cast_to_uint16_vector__struct(data,re_encode_char,string_to_cellstr)
sz=size(data);data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output,re_encode_char,string_to_cellstr);
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function hash=ComputeNonCryptHash(data,varargin)
%Compute a non-cryptographic hash
%
% This function is intended to be fast, but without requiring a Java or mex implementation to do
% the actual hashing. It was *not* checked for any security flaws and is therefore probably
% vulnerable to most attacks.
% Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like
% storing passwords.
%
%syntax:
%  hash=ComputeNonCryptHash(data)
%  hash=ComputeNonCryptHash(data,HashLength)
%  hash=ComputeNonCryptHash(___,VersionFlag)
%
%data        The data to be hashed. Most common data types are allowed: uint*, int*, char, cell,
%            struct, double, or single (string is cast to a cell array of chars). The contents of
%            the nested data types (i.e. cell and struct) must also be one of the mentioned data
%            types.
%HashLength  The length of the hash (the number of bits). This value must be a multiple of 16. The
%            default is 256 bits. Depending on your input 64 bits might have some collisions, but
%            64 bits and higher should be safe.
%VersionFlag Either '-v1', '-v2'. This is provided for backwards compatibility. Version 1 of this
%            function has many hash collisions for scalar doubles and attempts to cast strings to
%            chars, instead of casting to a cell array of chars. Version 2 also decodes the UTF-8
%            chars from Octave and re-encodes them with UTF-16. That way the output is stable for
%            the Unicode code points.
%
%hash        The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020b     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 2.0.0
% Date:    2020-11-10
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
end
if nargin<2
    HashLength=256;
    UseVersion=2;
elseif nargin<3
    HashLength=varargin{1};
    UseVersion=2;
else
    HashLength=varargin{1};
    UseVersion=varargin{2};
    try
        if isa(UseVersion,'string'),UseVersion=char(UseVersion);end
    	if ~isa(UseVersion,'char'), error('trigger'); end
        UseVersion=str2double(UseVersion(3:end));
        if isnan(UseVersion) || round(UseVersion)~=UseVersion || UseVersion>2
            error('trigger');
        end
    catch
            error('HJW:ComputeNonCryptHash:InputIncorrect',...
                'Version input incorrect. Must be ''-v1'', ''-v2''.')
    end
end
if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
    error('HJW:ComputeNonCryptHash:InputIncorrect',...
        'Second input (hash length) must be a multiple of 16.')
end
    
try
    %Convert the input to an uint16 array (Nx1).
    re_encode_char_on_Octave=UseVersion>=2;%Chars will be decoded and encoded with UTF-16.
    string_to_cellstr=UseVersion>=2;%Cast strings to a cell of chars, instead of a char.
    data=cast_to_uint16_vector(data,re_encode_char_on_Octave,string_to_cellstr);
catch
    ME=lasterror; %#ok<LERR>
    if strcmp(ME.identifier,'MATLAB:nomem')
        %rethrow memory error
        rethrow(ME)
    else
        error('HJW:ComputeNonCryptHash:UnwindFailed',...
            'The nested input contains an unsupported data type.')
    end
end

%Extend to a multiple of HashLength bits. Padding with zeros is generally not advised, and the
%performance penalty for this extension (compared to padding with zeros) should be negligible.
if mod(numel(data),HashLength/16)
    extra=uint16(1:HashLength/16).'; extra(1:mod(numel(data),HashLength/16))=[];
    data=[data;extra];
end

%Add perturbation to the data and convert to 16xN logical. Then further perturb the intermediate
%result by circshifting every column (by col_index-1 positions) and doing an XOR in blocks with
%itself (by reshaping and transposing).
if UseVersion==1
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

%Reshape to HashLength cols and collapse the key size down to the hash length by counting the
%number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if nargin>3
    hash=data;%Return uint16 for the salting.
    return
end

%Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
%few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data,UseVersion);

%Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=circshift_by_col(data)
%Circshift every column by col_index-1 positions.
persistent LUT
sz=size(data);
if isempty(LUT) || any(size(LUT)<sz) || isempty(LUT{sz(1),sz(2)})
    %keep a record of ind, which speeds up similar sizes
    [x,y]=ndgrid(1:size(data,2),1:size(data,1));
    y=y.';x=x.';
    z=mod(x+y-2,size(data,1))+1;
    ind=sub2ind(size(data),z,x);
    if prod(sz)<=1000 %to prevent a memory-hog, only keep ind for small sizes
        LUT{sz(1),sz(2)}=ind;
    end
else
    ind=LUT{sz(1),sz(2)};
end
data=data(ind);
end
function data=ComputeNonCryptHash_add_salt(data,UseVersion)
%Apply a few transformations to the hash to increase the spread.
%If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
%A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
%Avoid an infinite recursion by using a fourth input:
salt=ComputeNonCryptHash(data,saltHashLength,'-v1',[]);
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
if UseVersion>1
    salt=salt(end:-1:1);
end
data=mod(double(data).*double(salt),1+2^16);
data=uint16(data);
end
function hash=ComputeNonCryptHash_dec2hex(data)
%Look up the precomputed dec2hex for faster conversion.
persistent LUT
if isempty(LUT)
    LUT=upper(dec2hex(0:(-1+2^16),4));%even though the default should already upper case
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
hash=LUT(data,:);
end
function data=ComputeNonCryptHash_logical_to_uint16(data)
if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';
end
function data=ComputeNonCryptHash_shuffle_uint16(data)
%input should be uint16
base=65537;%base=(1+(2^16));
key=479001600;%key=1*2*3*4*5*6*7*8*9*10*11*12;
data = uint16(mod(double(data) * key , base));
end
function data=ComputeNonCryptHash_shuffle_uint16_inv(data)
base=65537;%base=(1+(2^16));
%key=1*2*3*4*5*6*7*8*9*10*11*12;
% %solution suggested by John D'Errico, https://www.mathworks.com/matlabcentral/answers/81859
% [G,C]=gcd(key,base);invKey=mod(C,base);
invKey=1919;
data=uint16(mod(double(data) * invKey,base));
end
function data=ComputeNonCryptHash_uint16_to_logical(data)
%uint16 Nx1 vector in, logical 16xN array out
persistent LUT
if isempty(LUT)
    LUT=dec2bin(0:(-1+2^16))=='1';
    LUT=LUT.';
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
data=LUT(:,data);
end
function error_(options,varargin)
%Print an error to the command window, a file and/or the String property of an object.
%The error will first be written to the file and object before being actually thrown.
%
%The intention is to allow replacement of every error(___) call with error_(options,___).
%
% NB: the error trace that is written to a file or object may differ from the trace displayed by
% calling the builtin error function. This was only observed when evaluating code sections. 
%
%options.fid.boolean: if true print error to file (options.fid.fid)
%options.obj.boolean: if true print error to object (options.obj.obj)
%
%syntax:
%  error_(options,msg)
%  error_(options,msg,A1,...,An)
%  error_(options,id,msg)
%  error_(options,id,msg,A1,...,An)
%  error_(options,ME)               %equivalent to rethrow(ME)

%parse input to find id, msg, stack and the trace str
if isempty(options),options=struct;end%allow empty input to revert to default
if ~isfield(options,'fid'),options.fid.boolean=false;end
if ~isfield(options,'obj'),options.obj.boolean=false;end
if nargin==2
    %  error_(options,msg)
    %  error_(options,ME)
    if isa(varargin{1},'struct') || isa(varargin{1},'MException')
        ME=varargin{1};
        try
            stack=ME.stack;%use original call stack if possible
            trace=get_trace(0,stack);
        catch
            [trace,stack]=get_trace(2);
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
        [trace,stack]=get_trace(2);
        [id,msg]=deal('',varargin{1});
    end
else
    [trace,stack]=get_trace(2);
    if ~isempty(strfind(varargin{1},'%')) % id can't contain a percent symbol
        %  error_(options,msg,A1,...,An)
        id='';
        A1_An=varargin(2:end);
        msg=sprintf(varargin{1},A1_An{:});
    else
        %  error_(options,id,msg)
        %  error_(options,id,msg,A1,...,An)
        id=varargin{1};
        msg=varargin{2};
        if nargin>3
            A1_An=varargin(3:end);
            msg=sprintf(msg,A1_An{:});
        end
    end
end
ME=struct('identifier',id,'message',msg,'stack',stack);

%print to object
if options.obj.boolean
    msg_=msg;while msg_(end)==10,msg_(end)='';end%crop trailing newline
    if any(msg_==10)  % parse to cellstr and prepend error
        msg_=regexp_outkeys(['Error: ' msg_],char(10),'split'); %#ok<CHARTEN>
    else              % only prepend error
        msg_=['Error: ' msg_];
    end
    set(options.obj.obj,'String',msg_)
end

%print to file
if options.fid.boolean
    fprintf(options.fid.fid,'Error: %s\n%s',msg,trace);
end

%Actually throw the error.
rethrow(ME)
end
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers=1;end
if nargin<2, stack=dbstack;end
stack(1:skip_layers)=[];

%parse ML6.5 style of dbstack (name field includes full file location)
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp=stack(n).name;
        if strcmp(tmp(end),')')
            %internal function
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

%parse Octave style of dbstack (file field includes full file location)
persistent IsOctave,if isempty(IsOctave),IsOctave=exist('OCTAVE_VERSION', 'builtin');end
if IsOctave
    for n=1:numel(stack)
        [ignore,stack(n).file]=fileparts(stack(n).file); %#ok<ASGLU>
    end
end

%create actual char array with a (potentially) modified stack
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
function im=hex2im(str)
%only supports 0-9,A-F and space
glyphs=hex2im_glyphs;
glyphs{32}=ones(size(glyphs{48}));
im=logical(cell2mat(glyphs(double(str))));
end
function TxtIm=hex2im_glyphs
%https://www.mathworks.com/matlabcentral/fileexchange/19896-convert-text-to-an-image
TxtIm=cell(1,70);
TxtIm{48}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{49}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{50}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{51}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{52}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{53}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{54}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{55}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{56}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{57}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{65}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
    1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{66}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{67}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{68}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{69}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
TxtIm{70}=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
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
        %k contains the indices of the beginning of each match
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
        k=k(:);%enforce direction
        remove=(n+1):numel(pattern);
        idx=bsxfun_plus(k,remove-1);
        out(idx(:))=[];
    end
end
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
%Input is implicitly converted to a row-vector.

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if nargin==1
    encode_as_UTF16=~isOctave;
end
if encode_as_UTF16
    if all(unicode<65536)
        str=uint16(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %encode as UTF-16
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
        str=char(str);% Conversion to char could trigger a conversion range error in Octave.
    end
else
    if all(unicode<128)
        str=char(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %encode as UTF-8
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
U=double(unicode)-65536;%convert to double for ML6.5
U=dec2bin(U,20);
str=bin2dec(['110110' U(1:10);'110111' U(11:20)]).';
end
function str=unicode_to_UTF8(unicode)
%Convert a single character to UTF-8 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-8
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
bytes=find(pers.limits.lower<unicode & unicode<pers.limits.upper);
str=pers.scheme{bytes};
scheme_pos=pers.scheme_pos{bytes};
b=dec2bin(unicode,pers.bits(bytes));
str(scheme_pos)=b;
str=bin2dec(str.').';
end
function [unicode,isUTF8,assumed_UTF8]=UTF8_to_unicode(UTF8,print_to)
%Convert UTF-8 to the code points stored as uint32
%Plane 16 goes up to 10FFFF, so anything larger than uint16 will be able to hold every code point.
%
%If there a second output argument, this function will not return an error if there are encoding
%error. The second output will contain the attempted conversion, while the first output will
%contain the original input converted to uint32.
%
%The second input can be used to also print the error to a GUI element or to a text file.
if nargin<2,print_to=[];end
return_on_error= nargout==1 ;

UTF8=uint32(UTF8);
[assumed_UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error);
if strcmp(flag,'success')
    isUTF8=true;
    unicode=assumed_UTF8;
elseif strcmp(flag,'error')
    isUTF8=false;
    if return_on_error
        error_(print_to,ME)
    end
    unicode=UTF8;%return input unchanged (apart from casting to uint32)
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
    multibyte=UTF8>=val & UTF8<256;%Exclude the already converted chars
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
            %To view the binary data, you can use this: bin=bin(:)';
            %Remove binary header (3 byte example):
            %1110xxxx10xxxxxx10xxxxxx
            %    xxxx  xxxxxx  xxxxxx
            if ~strcmp(header_bits,bin(header_locs))
                %Check if the byte headers match the UTF8 standard
                flag='error';
                if return_on_error,return,end
                continue %leave unencoded
            end
            bin(header_locs)='';
            if ~isOctave
                S3=uint32(bin2dec(bin  ));
            else
                S3=uint32(bin2dec(bin.'));%Octave needs an extra transpose
            end
            %Perform actual replacement
            UTF8=PatternReplace(UTF8,S2{n},S3);
        end
    end
end
end