from __future__ import annotations

import nox
from nox.sessions import Session


@nox.session
def docs(session: Session) -> None:
    """Build the documentation. Simply execute `nox -rs docs -- serve` to locally build and serve the docs."""
    session.install(".[docs]")
    session.chdir("docs")
    session.run("sphinx-build", "-M", "html", "source", "_build")

    if session.posargs:
        if "serve" in session.posargs:
            print("Launching docs at http://localhost:8000/ - use Ctrl-C to quit")
            session.run("python", "-m", "http.server", "8000", "-d", "_build/html")
        else:
            print("Unsupported argument to docs")
