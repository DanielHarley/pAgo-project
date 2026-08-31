from __future__ import annotations

import platform
import sys

import Bio
import numpy
import pandas
import plotly
import pyhmmer
import scipy
import sklearn


def main() -> None:
    print("Python:", sys.version)
    print("Executable:", sys.executable)
    print("Platform:", platform.platform())
    print("numpy:", numpy.__version__)
    print("pandas:", pandas.__version__)
    print("scikit-learn:", sklearn.__version__)
    print("scipy:", scipy.__version__)
    print("biopython:", Bio.__version__)
    print("plotly:", plotly.__version__)
    print("pyhmmer:", pyhmmer.__version__)


if __name__ == "__main__":
    main()
