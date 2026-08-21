# M8.8 clean room: operational instructions

This file is operational only. **`PROTOCOL.md` in this directory governs everything
substantive**, and where this file and the protocol appear to differ, the protocol wins and
the difference is something to report rather than resolve on your own.

## 0. Eligibility, read this before anything else

You are the implementer for a reproduction whose entire value depends on the firewall
described in the protocol's § 4. Eligibility is permanent and is about what a context has
been shown, never about what it still holds, so no summary, teardown, or fresh prompt
restores it.

**Stop immediately and report, without producing any work, if you have previously been shown
any of the following**: the M8.3 implementation, its outputs, or its method note; the
mode-identity-theory artifact or its published pages; the answer packet for this task; or any
published source carrying computed values for the quantity this task reproduces. The
protocol's § 4.3 carries the authoritative forbidden list; this is the short form for the
check you run before reading anything else.

If you are a genuinely fresh context, continue.

## 1. What this room is

A sealed working directory. It has no git relationship to any repository and there is nothing
here but the four files listed below. The room was built by extracting bytes from a frozen
commit; every technical file in it is byte-identical to a published, hash-pinned object.

| File | What it is |
| --- | --- |
| `PROTOCOL.md` | the frozen reproduction protocol, the governing document |
| `m8_5a_packet.json` | the group packet: exact generators and the coefficient-field convention |
| `m8_8_construction_packet.json` | the construction packet: the based chain complex you build from |
| `TASK.md` | this file |

**Scope boundary, and it is hard.** Work only inside this directory. Do not read, list, open,
or search any path outside it, and do not use the network. Nothing outside this directory is
part of your task, and material adjacent to this room on the filesystem is exactly the
material the firewall exists to keep away from you. If something you need appears to be
missing, that is a stopping condition under § 4 below, never a reason to look further afield.

## 2. Inputs

Your permitted construction inputs are the two packets in this directory plus your own
generic knowledge of algebra over group rings and number fields, representation theory of
finite groups, algebraic topology of chain complexes, and standard exact-arithmetic
technique. The protocol's § 4.3 carries the authoritative permitted list and the forbidden
list; read it before you start and treat the forbidden list as absolute, in every form,
including reading for orientation.

**No external references.** No web access, no literature lookup, no consulting a published
source of any kind. This is deliberate and the protocol's § 4.1 explains why: the literature
that carries the structural input you need tends to print the answer beside it, frequently on
the same page. If a step genuinely cannot be completed from the permitted inputs, that is a
recorded outcome under § 4 below, not a reason to open access.

## 3. What to produce, and in what order

The protocol's § 8 steps 4 and 5 fix the order, and the order is itself the evidence. Produce
the method-and-gate manifest **first and completely**, as its own file, before you write
implementation code:

1. **`METHOD_AND_GATE_MANIFEST.md`**: the route you selected within the protocol's § 6 class,
   the instantiated route-native gates of its § 7, the conventions you consumed, and your
   declared native orientation per its § 5.4. Fix these before implementing; the protocol
   treats a manifest written after the fact as a different and weaker object. **Once
   implementation begins, the manifest is immutable. If implementation exposes a mismatch
   with it, either change the implementation to conform to the still-intended manifest or
   STOP and report the mismatch; never amend the manifest to match implementation.**
2. **The implementation**, as source files in this directory.
3. **`ENVIRONMENT.md`**: interpreter version, library versions, platform, and anything else
   needed to rerun you.
4. **The derivation artifacts** required by the protocol's § 7.
5. **`RAW_OUTPUT.json`**: your result, in the protocol's § 5.5 schema, exactly.
6. **`CONSULTED_FILES.md`**: every file you read, listed explicitly, with a statement that no
   external references were consulted, stated affirmatively rather than left to be inferred
   from an absence.

Do not run `git` and do not attempt to commit anything. Someone outside the room copies your
files out and commits them, and that commit happens **before** any quarantined object is
unsealed. You will never be shown the comparison reference, and that is the design working
rather than a limitation of your access.

## 4. Gates, mutation tests, and stopping rules

Every gate you write must be able to fail. The protocol's § 9 requires this and it is the
single most common defect in work of this kind: for every line your code prints as a pass,
construct the mutation that makes it wrong, confirm the check goes red, and record the
mutation and its red result. A check whose two sides evaluate the same expression always
passes and is indistinguishable from a verified result to anyone reading later. Where a
quantity has no independent target to check against, label it as asserted rather than giving
it a check that cannot discriminate.

**Stopping conditions.** Each of these is a legitimate, recorded outcome, and reporting one
honestly is a better result than working around it:

- a step cannot be completed from the permitted inputs alone;
- the packets appear inconsistent, incomplete, or not to describe what the protocol says
  they describe;
- you find yourself needing a forbidden input to proceed;
- you realize mid-run that you are not eligible under § 0.

In any of these cases, stop, write what you established and exactly where it stopped, and
report. Do not reach outside the room, do not substitute a recalled value for a computed one,
and do not present a partial result as a complete one.

## 5. The one thing that matters most

Nothing you produce is graded on whether it agrees with anything. The protocol's § 8 lists
the outcome categories and several disagreement categories are recorded as results rather
than failures. What the record cannot survive is a number that came from memory instead of
from your construction, or a gate that could not have failed. Compute it, show the
derivation, and let the comparison fall where it falls.
