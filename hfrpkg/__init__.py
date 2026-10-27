"""hfrpkg: homodesmotic family reaction enthalpies.

Kept import-light on purpose: the version must be readable at build time
without RDKit / AaronTools installed. Legacy names are resolved lazily.
"""

__version__ = "0.2.0"

_LAZY = {
    "hfr_main": ("hfrpkg.cli", "main"),
    "compute_folder_main": ("hfrpkg.compute_single", "main"),
    "multi_main": ("hfrpkg.multi.hfr_driver", "main_cli"),
    "multi_run_main": ("hfrpkg.multi.run_multiple", "main_cli"),
    "multi_compute_main": ("hfrpkg.multi.compute_multiple", "main_cli"),
}

__all__ = ["__version__", *_LAZY]


def __getattr__(name):
    if name in _LAZY:
        import importlib

        mod, attr = _LAZY[name]
        return getattr(importlib.import_module(mod), attr)
    raise AttributeError(name)
