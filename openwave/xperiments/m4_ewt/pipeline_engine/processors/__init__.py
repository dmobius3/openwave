"""
Generic, reusable processors. Domain-agnostic.

Examples of what belongs here (future):
    CheckpointProcessor    -- save/restore ctx.data
    SpawnMonitorProcessor  -- launch an external monitor subprocess
    MetricsProcessor       -- write timings from ctx.diag to a sink
    ThrottleProcessor      -- decorator to run a processor every N steps

Physics-domain processors live in ../physics/ instead.
"""