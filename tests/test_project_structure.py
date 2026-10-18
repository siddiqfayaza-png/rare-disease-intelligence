"""Initial structure and import tests for rare_disease_ai package."""

import importlib
import unittest


class TestProjectStructure(unittest.TestCase):
    """Verify package scaffolding, exports, and basic metadata."""

    def test_package_import(self) -> None:
        """Verify rare_disease_ai package can be imported successfully."""
        pkg = importlib.import_module("rare_disease_ai")
        self.assertIsNotNone(pkg)

    def test_version_defined(self) -> None:
        """Verify __version__ string is defined and non-empty."""
        import rare_disease_ai

        self.assertTrue(hasattr(rare_disease_ai, "__version__"))
        self.assertIsInstance(rare_disease_ai.__version__, str)
        self.assertTrue(len(rare_disease_ai.__version__) > 0)

    def test_package_exports(self) -> None:
        """Verify __all__ contains declared exports."""
        import rare_disease_ai

        self.assertTrue(hasattr(rare_disease_ai, "__all__"))
        self.assertIn("__version__", rare_disease_ai.__all__)


if __name__ == "__main__":
    unittest.main()
