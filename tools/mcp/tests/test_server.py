"""The server as a client sees it: the installed entry point, over stdio."""

import asyncio
import json
import shutil
import sys
from pathlib import Path

import pytest
from mcp import Client, StdioServerParameters

ENTRY = Path(sys.executable).parent / "turbogap-mcp"


def run(session):
    async def go():
        params = StdioServerParameters(command=str(ENTRY), args=[])
        async with Client(params) as client:
            return await session(client)

    return asyncio.run(go())


@pytest.fixture(scope="module", autouse=True)
def entry_point():
    if not ENTRY.exists() and shutil.which("turbogap-mcp") is None:
        pytest.fail(f"{ENTRY} is missing: install the package (uv sync) before testing it")


def test_lists_the_phase_one_surface():
    async def session(client):
        tools = await client.list_tools()
        templates = await client.list_resource_templates()
        resources = await client.list_resources()
        return tools, templates, resources

    tools, templates, resources = run(session)
    assert {t.name for t in tools.tools} == {"keywords", "validate_input"}
    assert {str(r.uri) for r in resources.resources} == {"turbogap://keywords.md", "turbogap://examples"}
    assert [t.uri_template for t in templates.resource_templates] == ["turbogap://examples/{case}"]


def test_keywords_and_validate_over_the_wire():
    async def session(client):
        described = await client.call_tool("keywords", {"name": "exp_energy_scales_beg"})
        checked = await client.call_tool("validate_input",
                                         {"input_text": "md_step=1\n", "mode": "md"})
        return described, checked

    described, checked = run(session)
    assert described.structured_content["found"][0]["name"] == "exp_energy_scales"
    result = checked.structured_content
    assert not result["ok"]
    assert {i["keyword"] for i in result["issues"]} == {"md_step=1", "atoms_file", "pot_file"}


def test_example_resource_is_a_runnable_recipe():
    async def session(client):
        return await client.read_resource("turbogap://examples/co_predict")

    recipe = json.loads(run(session).contents[0].text)
    assert recipe["mode"] == "predict"
    assert {"source": "atoms_7176.xyz", "target": "atoms.xyz"} in recipe["links"]
    assert recipe["input"].startswith('atoms_file = "atoms.xyz"')


def test_mad_prompt_points_at_examples():
    async def session(client):
        return await client.get_prompt("mad_setup", {"observable": "xrd"})

    text = run(session).messages[0].content.text
    assert "turbogap://examples/xrd_mad" in text
    assert "validate_input" in text
