import json, sys, re
d = json.load(sys.stdin)
ti = d.get("tool_input", {})
text = ti.get("content", "") or ti.get("new_string", "")
if "—" in text:
    lines = [str(i+1) for i, l in enumerate(text.splitlines()) if "—" in l]
    print("Em dash (U+2014) is banned by CLAUDE.md. Found on line(s): "
          + ", ".join(lines) + ". Replace with standard punctuation.", file=sys.stderr)
    sys.exit(2)
sys.exit(0)
