import os

import pytest


def pytest_addoption(parser):
    parser.addoption("--answer_dir", help="Directory where answers are stored.")
    parser.addoption(
        "--answer_store",
        action="store_true",
        help="Generate new answers, but don't test.",
    )


@pytest.fixture()
def answer_store(request):
    return request.config.getoption("--answer_store")


@pytest.fixture(scope="session")
def answer_dir(request):
    ad = request.config.getoption("--answer_dir")
    if ad is None:
        pytest.skip("No answer directory")
    ad = os.path.abspath(ad)
    if not os.path.exists(ad):
        os.makedirs(ad)
    return ad


@pytest.fixture
def info_with_some_attr(info):
    if not hasattr(info, "some_attr"):
        pytest.skip("Missing attribute 'some_attr'")
    yield info
