import pytest


def pytest_addoption(parser):
    parser.addoption("--answer_store", action="store_true",
                     help="Generate new answers, but don't test. "
                          "Argument is the directory to store the answers to. "
                          "Default: None, which performs the test.")


@pytest.fixture()
def answer_store(request):
    return request.config.getoption('--answer_store')