"""Where the TurboGAP source tree is.

TURBOGAP_ROOT wins; otherwise the tree this package sits in (tools/mcp/turbogap_mcp).
"""

import os
from pathlib import Path


def root():
    env = os.environ.get("TURBOGAP_ROOT")
    path = Path(env) if env else Path(__file__).resolve().parents[3]
    if not (path / "tools" / "keyword_docs.py").is_file():
        raise RuntimeError(
            f"{path} is not a TurboGAP tree (no tools/keyword_docs.py); set TURBOGAP_ROOT")
    return path
