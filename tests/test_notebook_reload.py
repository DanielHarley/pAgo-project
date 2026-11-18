from __future__ import annotations

import importlib
import sys
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.notebook_reload import purge_pipeline_module_cache

TEST_PACKAGE_NAME = "pago_notebook_reload_test_package"
SIBLING_PACKAGE_NAME = f"{TEST_PACKAGE_NAME}_sibling"


def _write_helpers_module(*, package_directory: Path, value: str) -> None:
    (package_directory / "helpers.py").write_text(
        f"def get_value():\n    return {value!r}\n",
        encoding="utf-8",
    )


class PurgePipelineModuleCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.root_directory = Path(temporary_directory.name)

        # Without bytecode files, an edit made within the same second as the
        # first import cannot be masked by a stale .pyc file.
        original_dont_write_bytecode = sys.dont_write_bytecode
        sys.dont_write_bytecode = True
        self.addCleanup(
            setattr, sys, "dont_write_bytecode", original_dont_write_bytecode
        )

        self.package_directory = self.root_directory / TEST_PACKAGE_NAME
        self.package_directory.mkdir()
        (self.package_directory / "__init__.py").write_text("", encoding="utf-8")
        _write_helpers_module(
            package_directory=self.package_directory,
            value="old",
        )
        (self.package_directory / "consumer.py").write_text(
            f"from {TEST_PACKAGE_NAME}.helpers import get_value\n"
            "\n"
            "\n"
            "def read_value():\n"
            "    return get_value()\n",
            encoding="utf-8",
        )

        sibling_package_directory = self.root_directory / SIBLING_PACKAGE_NAME
        sibling_package_directory.mkdir()
        (sibling_package_directory / "__init__.py").write_text(
            "", encoding="utf-8"
        )

        sys.path.insert(0, str(self.root_directory))
        self.addCleanup(sys.path.remove, str(self.root_directory))
        self.addCleanup(self._remove_test_modules_from_cache)
        importlib.invalidate_caches()

    @staticmethod
    def _remove_test_modules_from_cache() -> None:
        for module_name in list(sys.modules):
            if module_name.startswith(TEST_PACKAGE_NAME):
                del sys.modules[module_name]

    def test_reload_alone_keeps_transitively_imported_module_stale(
        self,
    ) -> None:
        consumer_module = importlib.import_module(f"{TEST_PACKAGE_NAME}.consumer")
        self.assertEqual(consumer_module.read_value(), "old")

        _write_helpers_module(
            package_directory=self.package_directory,
            value="new-value",
        )
        consumer_module = importlib.reload(consumer_module)

        self.assertEqual(consumer_module.read_value(), "old")

    def test_purge_makes_imports_reread_transitively_imported_modules(
        self,
    ) -> None:
        consumer_module = importlib.import_module(f"{TEST_PACKAGE_NAME}.consumer")
        self.assertEqual(consumer_module.read_value(), "old")

        _write_helpers_module(
            package_directory=self.package_directory,
            value="new-value",
        )
        purge_pipeline_module_cache(package_name=TEST_PACKAGE_NAME)
        consumer_module = importlib.import_module(f"{TEST_PACKAGE_NAME}.consumer")
        helpers_module = importlib.import_module(f"{TEST_PACKAGE_NAME}.helpers")

        self.assertEqual(consumer_module.read_value(), "new-value")
        self.assertIs(consumer_module.get_value, helpers_module.get_value)

    def test_purge_removes_only_the_package_and_its_submodules(self) -> None:
        importlib.import_module(f"{TEST_PACKAGE_NAME}.consumer")
        importlib.import_module(SIBLING_PACKAGE_NAME)

        purged_module_names = purge_pipeline_module_cache(
            package_name=TEST_PACKAGE_NAME,
        )

        self.assertEqual(
            purged_module_names,
            [
                TEST_PACKAGE_NAME,
                f"{TEST_PACKAGE_NAME}.consumer",
                f"{TEST_PACKAGE_NAME}.helpers",
            ],
        )
        self.assertNotIn(TEST_PACKAGE_NAME, sys.modules)
        self.assertIn(SIBLING_PACKAGE_NAME, sys.modules)


if __name__ == "__main__":
    unittest.main()
