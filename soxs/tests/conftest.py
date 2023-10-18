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
