[![View ComputeNonCryptHash on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/76832-computenoncrypthash)

This function is intended to be fast, but without requiring a Java or mex implementation to do the actual hashing. It was **not** checked for any security flaws and is therefore probably vulnerable to most attacks.

Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like storing passwords.

This function will transform most common data types to a uint16 vector to apply the hash in an array operation. Changing the data type should change the hash. The allowed data types are uint*, int*, char, cell, struct, double, single, and string (which is cast to cell array of chars). The contents of the nested data types (i.e. cell and struct) must also be one of the mentioned data types.

Version 1.x of this algorithm attempts to cast string to char, instead of a cell array of chars. Version 1.x also has many hash collisions for scalar doubles. Version 2 will transcode the UTF-8 chars on Octave to UTF-16 (the Matlab standard), which ensures that the same Unicode code points as input will return the same hash.

Performance was tested with a 216553 items long English word list in both upper and lower case ([this list](http://web.archive.org/web/20061231134037/http://www.sitopreferito.it/html/all_english_words.html) with the two duplicates removed) and the numbers 0-1e6 as char and double. An additional test was performed with the images from the [Stanford Dog Dataset](http://web.archive.org/web/20191106220700/http://vision.stanford.edu/aditya86/ImageNetDogs/) containing 20580 images (the 89 duplicates were removed from [this](http://web.archive.org/web/20191106220631/http://vision.stanford.edu/aditya86/ImageNetDogs/images.tar) tar file before running this test). Timings below were determined on R2020b on Windows 10. For a comparison with other hash functions, see [this SE thread](https://softwareengineering.stackexchange.com/a/145633). Note that these tests are different from the relative performance comparison.

|Hash length|English words|Numbers (in char)|Numbers (in double)|Images|
|--|--|--|--|--|
|16 bits|50 &mu;s/hash<br>&thinsp;391&thinsp;630 collisions|50 &mu;s/hash<br>&thinsp;958&thinsp;488 collisions|56 &mu;s/hash<br>&thinsp;958&thinsp;487 collisions|106&thinsp;129 &mu;s/hash<br>4&thinsp;976 collisions|
|32 bits|58 &mu;s/hash<br>305 collisions|55 &mu;s/hash<br>31&thinsp;056 collisions|61 &mu;s/hash<br>120 collisions|105&thinsp;571 &mu;s/hash<br>0 collisions|
|48 bits|64 &mu;s/hash<br>144 collisions|63 &mu;s/hash<br>0 collisions|67 &mu;s/hash<br>16&thinsp;033 collisions|105&thinsp;288 &mu;s/hash<br>0 collisions|
|64 bits|60 &mu;s/hash<br>1 collisions|56 &mu;s/hash<br>0 collisions|62 &mu;s/hash<br>0 collisions|105&thinsp;366 &mu;s/hash<br>0 collisions|
|128 bits|65 &mu;s/hash<br>0 collisions|58 &mu;s/hash<br>0 collisions|69 &mu;s/hash<br>0 collisions|106&thinsp;289 &mu;s/hash<br>0 collisions|
|192 bits|70 &mu;s/hash<br>0 collisions|68 &mu;s/hash<br>0 collisions|72 &mu;s/hash<br>0 collisions|105&thinsp;932 &mu;s/hash<br>0 collisions|
|256 bits|83 &mu;s/hash<br>0 collisions|85 &mu;s/hash<br>0 collisions|92 &mu;s/hash<br>0 collisions|106&thinsp;187 &mu;s/hash<br>0 collisions|

Licence: CC by-nc-sa 4.0