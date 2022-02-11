from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef("void encode_frame(unsigned int n_windows, unsigned int n_fragments, unsigned int f_size, float (*MDframe)[3], double (*Fragments)[3], int *Encoding);")

ffibuilder.set_source("_encodeframe", """ #include "encodeframe.h" """, sources=["TrajSAencode/kabsch.c", "TrajSAencode/encodeframe.c"],
                      include_dirs=["./TrajSAencode"], libraries=["m", "gsl"])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)