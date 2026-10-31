import pytest

from turbogap_mcp.keywords import KeywordTable
from turbogap_mcp.tree import root


@pytest.fixture(scope="session")
def tree():
    return root()


@pytest.fixture(scope="session")
def table(tree):
    return KeywordTable(tree)
