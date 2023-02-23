from glob import glob
from pybind11.setup_helpers import Pybind11Extension, build_ext

def build(setup_kwargs):
    ext_modules = [
        Pybind11Extension(
            # NOTE changes in this section MUST be reflect in CMakeLists.txt
            "libshorah", 
            sources=sorted(glob("lib/src/*.cpp")),
            include_dirs=["lib/include"],
            libraries=["hts"],
            undef_macros=["HAVE_POPCNT"],
            extra_compile_args = ["-std=c++11"]
        ),
    ]

    setup_kwargs.update({
        "ext_modules": ext_modules,
        "cmdclass": {"build_ext": build_ext},
        "zip_safe": False,
    })
