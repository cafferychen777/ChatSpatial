"""Restore the sparse-matrix ``.A`` alias that SpaGCN 1.2.7 relies on.

SciPy 1.14 removed it; ChatSpatial <= 1.2.x restored it before importing the
original SpaGCN. ChatSpatial runs SpaGCN in a child process, so the alias is
installed at interpreter start-up via this module on PYTHONPATH.
"""

try:
    import scipy.sparse as _sp

    for _cls in (_sp.csr_matrix, _sp.csc_matrix, _sp.coo_matrix):
        if not hasattr(_cls, "A"):
            _cls.A = property(lambda self: self.toarray())
except ImportError:  # pragma: no cover
    pass
