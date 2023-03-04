import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--answer_store",
        action="store",
        default=False,
        help="Generate new answers, but don't test.",
    )


@pytest.fixture()
def answer_store(request):
    return request.config.getoption("--answer_store")


@pytest.fixture
def info_with_some_attr(info):
    if not hasattr(info, "some_attr"):
        pytest.skip("Missing attribute 'some_attr'")
    yield info
