from . import baf
from . import rdr
from . import tools
from . import utils
from .config import VERSION


# what `__all__` does:
# https://stackoverflow.com/questions/44834/what-does-all-mean-in-python
# https://stackoverflow.com/questions/42950256/how-to-import-private-functions-with-import
__all__ = ["__version__"]
__version__ = VERSION
