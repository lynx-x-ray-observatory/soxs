set -x   # Show which command is being run

if [[ ${mode} == "testing" ]]; then

  # Download answers

  curl -OL http://hea-www.cfa.harvard.edu/~jzuhone/soxs_test_data_${ANSWER_VER}.tar.gz
  tar -zxvf soxs_test_data_${ANSWER_VER}.tar.gz

  # Set location of SOXS data
  mkdir -p $HOME/.config/soxs
  echo "[soxs]" > $HOME/.config/soxs/soxs.cfg
  echo "soxs_data_dir = ${GITHUB_WORKSPACE}/soxs_data" >> $HOME/.config/soxs/soxs.cfg
  echo "soxs_answer_dir = ${GITHUB_WORKSPACE}/soxs_test_data" >> $HOME/.config/soxs/soxs.cfg
  cat $HOME/.config/soxs/soxs.cfg

  if ! [[ ${whichos} == "windows-latest" ]]; then
    # Set location of ATOMDB data
    mkdir -p $HOME/atomdb
    touch $HOME/atomdb/userdata # to avoid
  fi

fi

# Install dependencies using mamba and pip

eval "$(micromamba shell hook --shell bash)"
micromamba shell init --shell bash --root-prefix=~/micromamba
micromamba activate test-env
micromamba install --yes -c conda-forge numpy pytest pip astropy scipy cython h5py tqdm pyyaml appdirs pandas regions

if [[ ${mode} == "wheels" ]]; then
  micromamba install --yes wheel setuptools
fi
# special case for SPEX
if [[ ${mode} == "testing" && ${npver} == "1" && ${pyver} == "3.11" ]]; then
  curl -OL https://zenodo.org/records/17313851/files/spex-3.08.02-Linux-Intel.tar.gz
  tar xvfz spex-3.08.02-Linux-Intel.tar.gz
fi

# Install soxs
if [[ ${mode} == "testing" ]]; then
  python -m pip install -e .
  if ! [[ ${whichos} == "windows-latest" || ${pyver} == "3.11" ]]; then
    python -m pip install pyatomdb
  fi
fi
