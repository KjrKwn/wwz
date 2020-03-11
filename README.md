# weighted wavelet Z-transform(WWZ)

reference: Foster (1996)

developed by Kojiro Kawana 190108

Main codes are written in C (wwz_c.c) for fast computing.

Note:  
There is another implementation of WWZ https://github.com/eaydin/WWZ.  
My code can compute faster because of the C-implemententaion.

## How to build
choose icc or gcc in `Makefile`. Then, type

```
make
```

Then library file `libwwz.so` is created. The library is loaded by
wwz_c.py

## Tutorial for python wrapper
You can use wwz by
```
import wwz
wwz.wwz()
```
