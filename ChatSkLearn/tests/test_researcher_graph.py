import pytest
import os
import sys

# Add root path
current_file_path = os.path.abspath(__file__)
project_root = os.path.abspath(os.path.join(current_file_path, "../.."))
if project_root not in sys.path:
    sys.path.append(project_root)

from app.graphs.researcher import create_researcher_graph
from app.graphs.states import ResearcherState

@pytest.mark.asyncio
async def test_researcher_graph_basic():
    graph = create_researcher_graph()
    input_state = ResearcherState(question="How does cross-validation work in scikit-learn?")
    output = await graph.ainvoke(input_state)
    
    assert isinstance(output, dict)
    assert "documents" in output
    assert isinstance(output["documents"], list)
