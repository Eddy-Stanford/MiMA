import os
import sys

if __package__ in (None, ""):
    # run as `python3 tools/mimadoc`: make the package importable
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.dont_write_bytecode = True

from mimadoc.cli import main  # noqa: E402

sys.exit(main())
