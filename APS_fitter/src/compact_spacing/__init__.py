"""Mean period-spacing inference for pulsating compact stars."""
from .core import FitResult, Splitting, correct_to_m0, fit_sequence, gls_line

__all__ = ["FitResult", "Splitting", "correct_to_m0", "fit_sequence", "gls_line"]
__version__ = "0.1.0"
