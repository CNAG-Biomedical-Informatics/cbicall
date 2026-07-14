"""Adapter for the resource-bundle installer shipped with CBIcall."""

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import List, Optional

from .paths import resource_installer_path


def _load_installer(path: Optional[Path] = None) -> ModuleType:
    installer = Path(path) if path is not None else resource_installer_path()
    if not installer.is_file():
        raise RuntimeError(f"CBIcall resource installer was not found: {installer}")

    spec = importlib.util.spec_from_file_location("_cbicall_resource_installer", installer)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load the CBIcall resource installer: {installer}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_resource_installer(argv: Optional[List[str]] = None) -> int:
    module = _load_installer()
    return int(module.main(argv))
