"""
Utilities for managing output and logging.
"""

import io
import sys
import warnings
from contextlib import contextmanager, redirect_stdout, redirect_stderr


@contextmanager
def suppress_output():
    """Context manager to suppress stdout, stderr, and warnings during analysis."""
    # Suppress warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        # Capture stdout and stderr
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        
        try:
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                yield
        finally:
            # Optionally log captured output for debugging
            # This can be enabled if needed for development
            pass


class ProcessingError(Exception):
    """Custom exception for processing errors in spatial analysis."""
    pass 