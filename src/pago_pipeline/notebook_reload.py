from __future__ import annotations

import importlib
import sys
from typing import List

PIPELINE_PACKAGE_NAME = "src.pago_pipeline"


def purge_pipeline_module_cache(
    *,
    package_name: str = PIPELINE_PACKAGE_NAME,
) -> List[str]:
    """
    Remove a package and all of its submodules from the module cache.

    `importlib.reload` re-executes only the module it is given, so modules that
    are imported only transitively keep the version loaded earlier in a
    running Jupyter kernel, and the order of reload calls decides which
    version each module binds. Removing every cached module of the package
    instead makes the next `import` statements re-read all of the package's
    code from disk, with Python's import system resolving the dependency
    order.

    Objects created before the purge (for example class instances, `Enum`
    members, or values stored in notebook variables) keep referring to the
    old code. Changes to class definitions, dataclass fields, or `Enum`
    members therefore still require a kernel restart.

    Returns the names of the removed modules, in sorted order.
    """
    package_prefix = f"{package_name}."
    cached_module_names = sorted(
        module_name
        for module_name in sys.modules
        if module_name == package_name or module_name.startswith(package_prefix)
    )

    for module_name in cached_module_names:
        del sys.modules[module_name]

    # Drop finder caches so files added since the last import are found.
    importlib.invalidate_caches()

    return cached_module_names
