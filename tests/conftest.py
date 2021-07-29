"""Define app here for pytest-flask."""

import pytest

# Path required for VSCode test debugger to work (see GitHub rdkit issue #1276)
# Change the below <...> and uncomment if you want to use VSCode test debugger
# import os
# os.environ['PATH'] += r";<...>\Anaconda3\envs\MINE\Library\bin"
from api.run import create_app


@pytest.fixture
def app():
    """Create app. This fixture is required for pytest-flask plugin."""
    application = create_app()
    return application
