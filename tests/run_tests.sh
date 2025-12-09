set -x

micromamba activate test-env
if [[ ${npver} == "1" && ${pyver} == "3.11" ]]; then
  export SPEX90=$PWD/SPEX-3.08.02-Linux
  source $SPEX90/spexdist.sh
fi
python -m pytest -vv soxs/tests
