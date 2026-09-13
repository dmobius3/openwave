#!/bin/sh
# Headless clean-room launch (dev_docs/CLEAN_ROOM_STANDARDS.md § 3.4). macOS only: needs sandbox-exec.
#
# Usage: clean_room_launch.sh <room-dir> <opening-prompt-file> <transcript.jsonl> [model]
# Env:   CLEAN_ROOM_PYTHON  interpreter for the room (default: python3 on PATH)
#        CLEAN_ROOM_PACKAGE the platform package that must NOT import from the room (default: openwave)
#
# The room holds only the audited packet. The session loads no instruction file, memory, MCP server
# or skill, has four tools, denies anything that would prompt, and runs code only through ./py,
# which sandbox-exec confines to the room and the interpreter's own install, with no network.
set -eu
[ $# -ge 3 ] || { echo "usage: $0 <room-dir> <opening-prompt-file> <transcript.jsonl> [model]"; exit 64; }
ROOM=$(cd "$1" && pwd -P)
PROMPT=$(cd "$(dirname "$2")" && pwd -P)/$(basename "$2")
OUT=$3
MODEL=${4:-opus}
PY=$(command -v "${CLEAN_ROOM_PYTHON:-python3}")
PKG=${CLEAN_ROOM_PACKAGE:-openwave}
PREFIX=$("$PY" -c "import sys; print(sys.base_prefix)")
SITE=$("$PY" -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")
USERTMP=$(cd "${TMPDIR:-/tmp}" && pwd -P)

# § 3.1: no agent-instruction file anywhere on the room's ancestor path
d=$ROOM
while [ "$d" != "/" ]; do
  for f in "$d/CLAUDE.md" "$d/.claude/CLAUDE.md"; do
    [ -e "$f" ] && { echo "ABORT: instruction file on the room's ancestor path: $f"; exit 2; }
  done
  d=$(dirname "$d")
done

cat > "$ROOM/.room.sb" <<EOF
(version 1)
(allow default)
(deny file-read* (subpath "/Users"))
(deny file-read* (subpath "/private/tmp"))
(deny file-read* (subpath "$USERTMP"))
(deny file-read* (subpath "/private/var/tmp")) ; system scratch, allow-local-path
(deny file-write* (subpath "/Users"))
(deny file-write* (subpath "/private/tmp"))
(deny file-write* (subpath "$USERTMP"))
(deny file-write* (subpath "/private/var/tmp")) ; system scratch, allow-local-path
(allow file-read* (subpath "$PREFIX"))
(allow file-read* (subpath "$SITE"))
(allow file-read* file-write* (subpath "$ROOM"))
(deny network*)
EOF

cat > "$ROOM/py" <<EOF
#!/bin/sh
# room interpreter: OS-sandboxed to this room, no site processing, one thread, low priority
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH=$SITE
exec nice -n 15 sandbox-exec -f "$ROOM/.room.sb" "$PY" -S "\$@"
EOF
chmod +x "$ROOM/py"

# § 3.2: the platform package must not import from the room, and the sandbox must refuse an outside read
cd "$ROOM"
./py -c "import $PKG" >/dev/null 2>&1 && { echo "ABORT: $PKG is importable from the room"; exit 2; }
"$PY" -S -c "import os; os.listdir('$HOME')" >/dev/null 2>&1 || { echo "ABORT: the outside-read probe fails even unsandboxed, so it cannot test the sandbox"; exit 2; }
./py -c "import os; os.listdir('$HOME')" >/dev/null 2>&1 && { echo "ABORT: the sandbox did not refuse an outside read"; exit 2; }
./py -c "import numpy" >/dev/null 2>&1 || echo "note: numpy does not import in the room"

claude -p "$(cat "$PROMPT")" --model "$MODEL" --restricted --tools "Read,Edit,Write,Bash" \
  --strict-mcp-config --disable-slash-commands --permission-prompts none \
  --allowedTools "Edit(./**)" "Write(./**)" "Bash(./py *)" \
  --no-session-persistence --output-format stream-json --verbose < /dev/null > "$OUT"
echo "done: transcript at $OUT"
