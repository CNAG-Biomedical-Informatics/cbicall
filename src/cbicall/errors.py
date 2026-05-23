class CBIcallError(Exception):
    """Base class for CBIcall-specific failures."""


class ParameterValidationError(ValueError, CBIcallError):
    """Raised when user parameters are invalid or semantically incompatible."""


class WorkflowResolutionError(RuntimeError, CBIcallError):
    """Raised when a workflow cannot be resolved from the validated registry/config."""


class WorkflowExecutionError(RuntimeError, CBIcallError):
    """Raised when a resolved workflow cannot be launched or exits unsuccessfully."""
