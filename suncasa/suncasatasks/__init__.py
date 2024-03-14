# __init__.py of suncasa.suncasatasks
import warnings

# Always available modules
from .ptclean6 import ptclean6
from .subvs import subvs
from .concateovsa import concateovsa

__all__ = ['ptclean6', 'subvs', 'concateovsa']

# Conditional import for modules that depend on 'aipy-eovsa'
try:
    from .importeovsa import importeovsa
    __all__.append('importeovsa')
except ImportError:
    warnings.warn(
        "Failed to import 'importeovsa'. This module requires 'eovsapy'. "
        "Please install 'eovsapy' to use 'importeovsa'."
    )

try:
    from .calibeovsa import calibeovsa
    __all__.append('calibeovsa')
except ImportError:
    warnings.warn(
        "Failed to import 'calibeovsa'. This module requires 'eovsapy'. "
        "Please install 'eovsapy' to use 'calibeovsa'."
    )
