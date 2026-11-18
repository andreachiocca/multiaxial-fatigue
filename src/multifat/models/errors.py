from __future__ import annotations


class MissingParamsError(ValueError):
    """A model cannot run on a material because material constants are missing."""

    def __init__(self, model: str, missing, hint: str = ""):
        self.model, self.missing = model, sorted(missing)
        msg = f"Missing required material params for {model}: {self.missing}"
        super().__init__(msg + (f". {hint}" if hint else ""))


class NativeLifeError(ValueError):
    """A model's native life law cannot be fitted on this material (its damage values still can)."""


class LoadingScopeError(ValueError):
    """A case lies outside the loading scope of a criterion (the model is skipped, not failed)."""


class MaterialScopeError(ValueError):
    """The material constants lie outside the validity range of a criterion (model skipped).

    Raised when the source restricts a criterion to a range of material properties, for example
    RFB to 0.577 < t_-1/s_-1 < 0.866. It is reported like a loading-scope skip (eligibility
    ``out_of_scope``), but it is not a ``LoadingScopeError``: it does not depend on the load, so
    the load-scaling searches must not treat it as an out-of-scope trial load.
    """


class DesignSetError(ValueError):
    """A material lacks the design tests a model's master curve is fitted on (model skipped)."""
