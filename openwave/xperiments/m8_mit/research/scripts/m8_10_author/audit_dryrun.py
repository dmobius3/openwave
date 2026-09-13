"""Audits a dry-run transcript (stream-json) for anything that reached outside the room.

Kept out of the room on purpose: it names the answer-bearing paths. Flags any web or agent tool, any file
path outside the room, and any shell command that names a path outside the room, climbs out with '../', or
carries an answer-bearing token. A synthetic out-of-room read is run first and must be flagged.
"""
import json
import re
import sys

ROOM = '<home>/m810_dryrun'
TOKENS = ['Desktop', 'mode-identity-theory', 'openwave', '.claude/projects', 'M8_10', 'RESULTS', 'surviving', '/tmp/', '../']
OUTSIDE_TOOLS = {'WebFetch', 'WebSearch', 'Task', 'Agent'}


def tool_uses(obj):
    if isinstance(obj, dict):
        if obj.get('type') == 'tool_use':
            yield obj
        for v in obj.values():
            yield from tool_uses(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from tool_uses(v)


def audit(lines):
    flags, n = [], 0
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        for tu in tool_uses(event):
            n += 1
            name, inp = tu.get('name'), tu.get('input') or {}
            if name in OUTSIDE_TOOLS:
                flags.append((name, 'web or agent tool'))
                continue
            for key in ('file_path', 'path', 'notebook_path'):
                p = str(inp.get(key) or '')
                if p and not p.startswith(ROOM) and not p.startswith('.') and '/' in p:
                    flags.append((name, f'{key} outside the room: {p}'))
            if name in ('Glob', 'Grep') and str(inp.get('pattern', '')).startswith('/') and not str(inp['pattern']).startswith(ROOM):
                flags.append((name, f"pattern outside the room: {inp['pattern']}"))
            if name == 'Bash':
                cmd = str(inp.get('command', ''))
                hits = [t for t in TOKENS if t in cmd]
                hits += [q for q in re.findall(r'(/(?:Users|private|tmp|etc|var|opt)\S*)', cmd) if not q.startswith(ROOM)]
                hits += [q for q in re.findall(r'~/\S*', cmd) if not q.startswith('~/m810_dryrun')]
                if hits:
                    flags.append((name, f'command reaches outside: {hits} :: {cmd[:140]}'))
    return n, flags


if __name__ == '__main__':
    control = json.dumps({'type': 'assistant', 'message': {'content': [
        {'type': 'tool_use', 'name': 'Read', 'input': {'file_path': '<home>/Desktop/MIT/OpenWave/M8_10/RESULTS.md'}},
        {'type': 'tool_use', 'name': 'Bash', 'input': {'command': 'cat ../RESULTS.md'}}]}})
    assert len(audit([control])[1]) == 2, 'positive control not flagged'
    path = sys.argv[1] if len(sys.argv) > 1 else ROOM + '/transcript.jsonl'
    n, flags = audit(open(path, encoding='utf-8'))
    print(f'control flagged as required; {n} tool calls audited; {len(flags)} flags')
    for f in flags:
        print('  FLAG', f)
    print('CLEAN' if not flags else 'NOT CLEAN')
