# pipeline_engine

Domain-agnostic staged pipeline engine for M4. Runs a list of stateless
processors against a shared Context. Knows nothing about physics, waves,
or Taichi.

## What it is

The engine is a composition root for a simulation: you build a Pipeline
from processors, hand it to a Runner, and the Runner executes it step by
step. Each processor has one responsibility and is stateless; runtime
state lives in a typed FeatureBag inside the Context.

## Why

M4's wave_engine.py is a single monolithic solver where the update,
boundary condition, wave-center drive, trackers, and buffer swap all live
in one kernel. Adding a mode means editing several files; comparing two
modes means two branches of code. The engine offers composition instead:
a pipeline is a list of processors, and swapping one processor is the only
change needed to try a different physics model. wave_engine.py and the
launcher are left as they are; the engine sits alongside them.

## Quick start

    from openwave.xperiments.m4_ewt.pipeline_engine import (
        Pipeline, Runner, Stage, BaseProcessor,
        JsonSessionSink, LogProcessor,
    )

    class Increment(BaseProcessor):
        name = "Increment"
        stage = Stage.UPDATE
        order = 10
        # requires/provides default to empty

        def process(self, ctx):
            # runtime state goes in ctx.data, not self
            pass

    class MyPipeline(Pipeline):
        def __init__(self):
            super().__init__()
            self.add(Increment())

    runner = Runner({"session": JsonSessionSink("out/session.json")})
    ctx = runner.run(MyPipeline(), name="demo", max_steps=100)
    print(ctx.diag.errors)

## Design rules

Z1. Zero module-level state. Everything lives in Context or in instances.
Z2. Every mutation goes through a processor. No side effects from the runner.
Z3. Explicit requires/provides graph, validated at Pipeline.build().
Z4. setup/teardown are lifecycle hooks; process() is per-step.
Z5. The runner is an interface (IRunner). Runner is the default implementation.
Z6. Sinks are an interface (ILogSink). Several implementations included.
Z7. Measurement and logging are separate processors.
Z8. Parameters are typed and grouped separately from run identity.

## Stages

    PRE_UPDATE   -- runs before evolution. Example: seed, boundary pin, source.
    UPDATE       -- evolution. Example: Laplacian, nonlinearity, integrator.
    POST_UPDATE  -- runs after evolution. Example: boundary enforcement.
    MEASURE      -- read-only. Example: trackers, logging, visualization.

Processors with stage=None are lifecycle only: they run setup() once before
the first step and teardown() once after the last (allocators, subprocess
spawns, resource acquisition). There is no SETUP member in Stage.

## Context

Every processor sees one object, Context, with five sub-objects:

    run     -- immutable run identity (name, output_dir, seed)
    params  -- typed read-only parameter access (get_int, require_float, ...)
    sim     -- mutable simulation clock (dt, t, step, should_stop)
    data    -- FeatureBag: typed container, key is the class
    log     -- sinks by name (LogContext)
    diag    -- errors, warnings, timings (Diagnostics)

A processor reaches state through ctx.data.require(SomeFeature). Adding a
new feature never changes any signature.

## Statelessness

Processors are stateless by default. The engine raises if process()
rebinds an instance attribute (attributes are compared by object identity
before and after the call); mutating an attribute in place, such as
appending to a list held on self, is not detected. Runtime state belongs
in ctx.data. When memoization is genuinely required, set
`stateless = False` on the class.

The check runs on every process() call of a processor with
`stateless = True`. The `check_stateless` argument of `Runner` is not read.

## Error handling

ErrorPolicy selects what happens when a processor raises:

    FAIL_FAST  -- re-raise, stop immediately
    SOFT_STOP  -- log, set should_stop, run teardown, exit cleanly
    CONTINUE   -- log, skip to the next processor

Default is SOFT_STOP. The error is recorded in ctx.diag.errors either way.

## Folder layout

    context.py    -- FeatureBag, Params, Context, Diagnostics, sub-contexts
    pipeline.py   -- Stage, ErrorPolicy, IProcessor, BaseProcessor, Pipeline
    sinks.py      -- ILogSink, JsonSessionSink, LiveJsonSink, InMemorySink
    loggers.py    -- LogProcessor (writes a ctx-derived payload to a sink)
    runner.py     -- IRunner, Runner
    processors/   -- domain-agnostic processors (checkpoint, monitor spawn, ...)
    physics/      -- physics processors and feature types (separate layer)

The engine never imports physics. The physics layer imports the engine.

## Testing

Smoke test (no Taichi required):

    python -m openwave.xperiments.m4_ewt.pipeline_engine._smoke_test

Runs four scenarios: a normal pipeline, a pipeline with an injected failure
(SOFT_STOP in action), a pipeline with file sinks, and a pipeline reading a
caller-supplied feature (external_provides).

Physics demo (requires Taichi):

    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo --k 10 --geometry tetrahedron_10_locked --no-window --steps 3000

## Adding a new processor

1. Subclass BaseProcessor.
2. Set name, stage, order, requires, provides.
3. Put configuration in __init__ (immutable).
4. Put runtime state in ctx.data (typed feature).
5. Implement setup/process/teardown as needed.

If the processor needs a feature that is not provided by any other
processor, either declare it via Pipeline(external_provides=...) and pass
an instance through Runner.run(..., initial_features=[...]), or add a
lifecycle processor that provides it in setup().
