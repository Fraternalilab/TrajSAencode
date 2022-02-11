import setuptools

setuptools.setup(
    name="TrajSAencode",
    version="0.0.1",
    packages=setuptools.find_packages(exclude=["docs","tests", ".gitignore", "README.rst","DESCRIPTION.rst"]),
    python_requires=">=3.6",
    cffi_modules=["TrajSAencode/kabsch_extension_build.py:ffibuilder",
                  "TrajSAencode/encodeframe_extension_build.py:ffibuilder"],
    install_requires=['cffi', "mdtraj"]
)