This function is intended to be fast, but without requiring a Java or mex implementation to do the actual hashing. It was **not** checked for any security flaws and is therefore probably vulnerable to most attacks.

Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like storing passwords.

This function will transform most common data types to a uint16 vector to apply the hash in an array operation. Changing the data type should change the hash. The allowed data types are uint*, int*, char, cell, struct, double, single, and string (which is cast to char). The contents of the nested data types (i.e. cell and struct) must also be one of the mentioned data types.

Performance was tested with a 433110 items long English word list ([this list](http://web.archive.org/web/20061231134037/http://www.sitopreferito.it/html/all_english_words.html) with the two duplicates removed) and the numbers 1-216553 as char. Timings below were determined on R2020a on Windows 10. Time per hash may increase by a factor of 5 for older releases and 20 for Octave, and by a factor of 200 for data types other than uint16 and char. For a comparison with other hash functions, see [this SE thread](https://softwareengineering.stackexchange.com/a/145633).

|Hash length|English words|Numbers|
|--|--|--|
|16 bits|42 &mu;s/hash<br>407198 collisions|42 &mu;s/hash<br>192790 collisions|
|32 bits|50 &mu;s/hash<br>442 collisions|50 &mu;s/hash<br>15155 collisions|
|48 bits|56 &mu;s/hash<br>14 collisions|57 &mu;s/hash<br>0 collisions|
|64 bits|54 &mu;s/hash<br>0 collisions|52 &mu;s/hash<br>0 collisions|
|128 bits|57 &mu;s/hash<br>0 collisions|54 &mu;s/hash<br>0 collisions|
|256 bits|118 &mu;s/hash<br>0 collisions|115 &mu;s/hash<br>0 collisions|

Licence: CC by-nc-sa 4.0