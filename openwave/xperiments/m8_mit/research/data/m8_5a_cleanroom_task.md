# TASK: implement the reproduction specified in PROTOCOL.md

## What you are doing

Implement the reproduction that `PROTOCOL.md` specifies, using the group input in
`GROUP_INPUT.md` and the machine-readable packet it names. Those two files plus this one are
the complete specification.

## Environment

Everything you need is in this directory. Nothing outside it is needed. Web search and fetch
tools are disabled by design.

**This directory is your whole world, and that is a rule, not a description of the machine.**
The filesystem outside it is not sandboxed away from you; staying inside it is your
obligation. Do not read, list, or search any path outside this directory, including any
repository you may find on this disk or hear named in instructions that loaded automatically.
If you nonetheless open something outside, that is recoverable and concealing it is not:
record it in `manifest.md` and say so plainly in your final message.

`PROTOCOL.md` is a verbatim copy of a frozen document. Its relative links point into a
repository you do not have and must not go looking for. Nothing behind them is needed, and the
protocol text itself is complete for your purposes.

If at some point you want a reference table of characters, dimensions, graph distances, or
first occurrences for the group in `GROUP_INPUT.md`: that table is the answer key for this
task. Deriving it is the work. Looking it up, reconstructing it from memory of a published
source, or asking for it invalidates the run.

## Deliverables, all written into this directory

| File | Content |
| --- | --- |
| `m8_5a_reproduction.py` | the implementation, including a `--mutation-tests` mode |
| `raw_output.txt` | the run's stdout, verbatim |
| `result.json` | the output schema fixed in `PROTOCOL.md` § 5 |
| `environment.md` | interpreter and library versions, OS, hardware, seeds if any |
| `manifest.md` | every file and reference you consulted, including any that loaded automatically without your asking |
| `method_note_draft.md` | equations first, then an equation-to-code map |

## Requirements beyond a plain reading of the protocol

1. **The § 6 coexact module RUNS.** This is pre-declared and not optional. If you cannot build
   it from your own harmonic analysis, its verdict is `not resolved` with the reason recorded.
   That is a permitted outcome. Opening external access to rescue it is not.
2. **Preference, not a requirement.** Construct the coefficient representation by building the
   symmetric square explicitly and taking traces, rather than by evaluating the character
   identity that G11 checks. Where the two routes genuinely differ, G11 is a cross-check;
   where they are the same route, G11 confirms the declared contract was applied but not its
   arithmetic. State which route you took, in the method note, either way.
3. **Every tolerance is fixed in the source** and justified in the method-note draft, before
   any comparison is run.
4. **Fail loud.** A NOT-FOUND cell, a gate failure, or an integer-nearness violation stops the
   run with a nonzero exit and a stated reason. Silent omission is not an option.
5. **Report the honest outcome.** A disagreement is a finding, not a failure of this task, and
   the result categories in the protocol are all reportable. Do not adjust a gate, a tolerance,
   or a construction to make an output look better.

## What "done" means

The mutation harness exits nonzero if any gate is uncovered or any mutation fails to redden its
target. Run it, include its output, and stop when it passes and every file above exists.
