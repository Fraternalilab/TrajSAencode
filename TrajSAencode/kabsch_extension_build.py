from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef("double wrmsd_kabsch(unsigned int size,  float (*Xarray)[3], float (*Yarray)[3]);")

ffibuilder.set_source("_kabsch", """ #include "kabsch.h" """, sources=["TrajSAencode/kabsch.c"],
                      include_dirs=["./TrajSAencode"], libraries=["m", "gsl"])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
