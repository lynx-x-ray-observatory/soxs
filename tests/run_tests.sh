set -x

micromamba activate test-env
python -m pytest -vv soxs/tests
