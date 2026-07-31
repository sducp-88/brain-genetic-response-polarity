from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
import platform
import sys


REQUIRED_PACKAGES = {
    "anndata": "anndata",
    "scanpy": "scanpy",
    "pyarrow": "pyarrow",
    "synapseclient": "synapseclient",
}


def main() -> None:
    failures: list[str] = []

    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")

    for distribution, module in REQUIRED_PACKAGES.items():
        try:
            import_module(module)
            installed_version = version(distribution)
            print(f"{distribution}: {installed_version}")
        except (ImportError, PackageNotFoundError) as error:
            failures.append(f"{distribution}: {error}")

    if failures:
        raise SystemExit(
            "Python environment verification failed:\n- "
            + "\n- ".join(failures)
        )

    print("Python environment verification passed.")


if __name__ == "__main__":
    main()
