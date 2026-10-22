"""Generated regions in the Markdown pages.  A region is

    <!-- mimadoc:KIND [ARG] -->

    ...generated text...

    <!-- mimadoc:end -->

and everything outside the regions is written by hand."""

import re

REGION_RE = re.compile(
    r"(<!-- mimadoc:(?P<kind>[a-z]+)(?: (?P<arg>\S+))? -->\n)(?P<body>.*?)(<!-- mimadoc:end -->)",
    re.S)


def regions(text):
    """[(kind, arg)] of the regions in a page, in order."""
    return [(m.group("kind"), m.group("arg")) for m in REGION_RE.finditer(text)]


def fill(text, generate):
    """Text with every region replaced by generate(kind, arg)."""
    def sub(m):
        body = generate(m.group("kind"), m.group("arg"))
        return "%s\n%s\n\n%s" % (m.group(1), body.strip("\n"), m.group(5))
    return REGION_RE.sub(sub, text)
