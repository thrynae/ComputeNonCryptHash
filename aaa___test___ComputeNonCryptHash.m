function aaa___test___ComputeNonCryptHash(varargin)
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
for k=1:size(hash,2)
    HashLength=numel(hash{1,k})*4;
    for n=1:numel(input)
        H=ComputeNonCryptHash(input{n},HashLength);
        if ~strcmp(hash{n,k},H)
            error('hash did not match: k=%d,n=%d\n',k,n)
        end
    end
end
clc,disp('test finished successfully')
end
function [cases,hash]=get_test_cases
%clc,for n=1:10,fprintf('    ''%s'',...\n',ComputeNonCryptHash(rand,256));end
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
hash={'342E','C7AB19F74F03A4B65141A30699AA30C60D2D0B958B737E4F69DFBC09C9CA3293';
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
    ComputeNonCryptHash(data,8);
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
function data=cast_to_uint16_vector(data)
%Linearize the input data and convert it to a uint16 vector
data=cast_to_uint16_vector__cell({data});
data([end-1 end])=[];%Remove the [1 1] that is added because of the wrapping in a cell
end
function data=cast_to_uint16_vector__cell(data)
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
            data{n}=cast_to_uint16_vector__char(data{n});
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n});
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n});
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
function data=cast_to_uint16_vector__char(data)
sz=size(data);data=data(:);
data=uint16(data);%Normal chars are already a uint16 internally (in Matlab).
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__struct(data)
sz=size(data);data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output);
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
%
%data       The data to be hashed. Most common data types are allowed: uint*, int*, char, cell,
%           struct, double, or single (string is cast to char). The contents of the nested data
%           types (i.e. cell and struct) must also be one of the mentioned data types.
%HashLength The length of the hash (the number of bits). This value must be a multiple of 16. The
%           default is 256 bits. Depending on your input 64 bits might have some collisions, but
%           64 bits and higher should be safe.
%
%hash       The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.1
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( http://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
else
    try
        if isa(data,'string'),data=char(data);end
    catch
        error('HJW:ComputeNonCryptHash:ConvertToChar',...
            'The required conversion from string to char failed.')
    end
end
if nargin<2
    HashLength=256;
else
    HashLength=varargin{1};
    if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
        error('HJW:ComputeNonCryptHash:InputIncorrect',...
            'Second input (hash length) must be a multiple of 16.')
    end
end

try
    %Convert the input to an uint16 array (Nx1).
    data=cast_to_uint16_vector(data);
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
%result by doing an XOR in blocks with itself (by reshaping and transposing).
data=ComputeNonCryptHash_shuffle_uint16(data);
data=ComputeNonCryptHash_uint16_to_logical(data);
data=xor(data,reshape(data,[],16).');

%Reshape to HashLength cols and collapse the key size down to the hash length by counting the
%number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if nargin>2
    hash=data;%Return uint16 for the salting.
    return
end

%Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
%few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data);

%Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=ComputeNonCryptHash_add_salt(data)
%Apply a few transformations to the hash to increase the spread.
%If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
%A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
salt=ComputeNonCryptHash(data,saltHashLength,[]);%avoid infinite recursion by using a third input
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
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