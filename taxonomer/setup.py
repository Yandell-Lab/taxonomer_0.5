from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(
    [
    Extension("scripts/cython/binner_PY",sources=["scripts/cython/binner_PY.pyx","src/binner.c","src/load_dbs.c","src/query_db.c","src/kmer_utils.c"],extra_link_args=["-lz","-lm","-fopenmp"]),
    Extension("scripts/cython/build_DBpy", sources=["scripts/cython/build_DBpy.pyx","src/build_db.c","src/kmer_counts.c","src/kmer_utils.c","src/query_db.c","src/load_dbs.c"],extra_link_args=["-lz","-lm","-fopenmp"]),
    Extension("scripts/cython/classify_READpy",sources=["scripts/cython/classify_READpy.pyx","src/classify.c","src/kmer_counts.c","src/kmer_utils.c","src/query_db.c","src/load_dbs.c","src/process_results.c"],extra_link_args=["-lz","-lm", "-fopenmp"]),
    
    #Extension("scripts/cython/binner_PY",sources=["scripts/cython/binner_PY.pyx","src/binner.c","src/load_dbs.c","src/query_db.c","src/kmer_utils.c"],extra_link_args=["-lz","-lm"]),
    #Extension("scripts/cython/build_DBpy", sources=["scripts/cython/build_DBpy.pyx","src/build_db.c","src/kmer_counts.c","src/kmer_utils.c","src/query_db.c","src/load_dbs.c"],extra_link_args=["-lz","-lm"]),
    #Extension("scripts/cython/classify_READpy",sources=["scripts/cython/classify_READpy.pyx","src/classify.c","src/kmer_counts.c","src/kmer_utils.c","src/query_db.c","src/load_dbs.c","src/process_results.c"],extra_link_args=["-lz","-lm"]),

    Extension("scripts/cython/binner_splitpy",sources=["scripts/cython/binner_splitpy.pyx"]),
    Extension("scripts/cython/binner_mergepy",sources=["scripts/cython/binner_mergepy.pyx"])
    ]
    )
)
