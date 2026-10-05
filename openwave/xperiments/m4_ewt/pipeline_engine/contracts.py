"""
Domain-agnostic processor contracts.

SourceTermInterface     -- a per-step additive contributor to one
                           buffer (psi_new). Domain-agnostic: does not
                           reference any feature.
ScatteringOperatorInterface -- a per-step magnitude-preserving
                           exchanger between two buffers. Domain-agnostic.

Both are IProcessor specializations. Neither is enforced at runtime;
the contract is verified by the tests that use them.
"""

from .pipeline import BaseProcessor, Stage


class SourceTermInterface(BaseProcessor):
    """
    A per-step additive source.

    Subclass contract:
      - Set name, order, requires.
      - Override contribute(ctx). The override adds via +=; it must not
        assign via = and must not write psi or psi_prev.
      - Do not override process(). The base class calls contribute().

    Additivity is verified by the tests, not enforced here.
    """

    name = "SourceTermInterface"
    stage = Stage.UPDATE
    order = 7

    def contribute(self, ctx) -> None:
        raise NotImplementedError(f"{type(self).__name__} must implement contribute(ctx)")

    def process(self, ctx) -> None:
        self.contribute(ctx)


class ScatteringOperatorInterface(BaseProcessor):
    """
    A per-step magnitude-preserving exchanger between two buffers.

    Subclass contract:
      - Set name, order, requires.
      - Override scatter(ctx). The override reads one field, writes
        another, and preserves the total magnitude per voxel.
      - Do not override process(). The base class calls scatter().

    Unitarity is verified by the tests, not enforced here.
    """

    name = "ScatteringOperatorInterface"
    stage = Stage.UPDATE
    order = 7

    def scatter(self, ctx) -> None:
        raise NotImplementedError(f"{type(self).__name__} must implement scatter(ctx)")

    def process(self, ctx) -> None:
        self.scatter(ctx)
