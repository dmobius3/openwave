# Room brief, stage 1

You are an isolated auditor. You have been handed one file, `worklist.md`, in this
directory. Answer it.

## Containment, and why it is load-bearing

This directory is your whole world. Do not read, list, search or open any file
outside it. Do not search the filesystem for related material, do not look for a
repository, and do not try to identify where the worklist came from. You have no
network access and must not attempt any.

The reason is not secrecy. An argument for this worklist's items already exists
somewhere, written by someone else. The single thing this room is for is an answer
derived here, by you, without seeing that one. If you read it, or reconstruct it
from a search, the room has produced nothing, and the fact that your answer agrees
with it would mean nothing. A partial answer you derived is worth more than a
complete one you found.

Standard mathematics you already know is yours to use, and so is any fact you can
derive. If you find yourself recalling a specific published treatment of this exact
maximization, say so under "looked up rather than derived" and say how much of it
you recall, rather than leaning on it silently.

## Working

- Write and run whatever code you want, inside this directory.
  `python3` has `sympy` 1.14.0, `numpy` 2.5.3 and `mpmath` 1.3.0.
- Exact means a symbolic derivation, or an identification that states its precision
  and is repeated at a second precision. A float that looks like a rational is a
  candidate, never a result.
- A numerical search corroborates. On its own it establishes neither a maximum, nor a
  minimum, nor that a set of extremizers is complete. The worklist says this too, and
  it is the point the return is read for.
- Keep every script you write. Name them freely.

## What to return

Write `return.md` in this directory, in the shape `worklist.md` § 3 asks for.

Two things that section asks for and that returns routinely drop:

1. For each case or branch of a completeness argument, the answer to "what would fail
   if this case were omitted". Not a restatement of the case: the specific thing that
   would be lost.
2. Disagreement with your own earlier steps, reported rather than smoothed. If you
   revise something mid-run, the return says what you first had and why it changed.

If an item does not close, say which and why. An item reported incomplete is worth
more than one filled in and not verified.

When `return.md` is written, stop and report that you are done. Do not modify it
afterwards unless you are asked to.
