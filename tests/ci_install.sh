set -x   # Show which command is being run

if [[ ${mode} == "testing" ]]; then

  # Download test data
  curl -OL http://yt-project.org/data/GasSloshingLowRes.tar.gz
  tar -zxf GasSloshingLowRes.tar.gz

  # Download answers

  curl -O http://hea-www.cfa.harvard.edu/~jzuhone/soxs_test_data_${ANSWER_VER}.tar.gz
  tar -zxvf soxs_test_data_${ANSWER_VER}.tar.gz

  # Set location of yt test data
  mkdir -p $HOME/.config/yt
  echo "[yt]" > $HOME/.config/yt/yt.toml
  echo "test_data_dir = \"${PWD}\"" >> $HOME/.config/yt/yt.toml
  cat $HOME/.config/yt/yt.toml

  # Set location of soxs data

  mkdir -p $HOME/.config/soxs
  echo "[soxs]" > $HOME/.config/soxs/soxs.cfg
  echo "soxs_data_dir = ${GITHUB_WORKSPACE}/soxs_data" >> $HOME/.config/soxs/soxs.cfg
  echo "soxs_answer_dir = ${GITHUB_WORKSPACE}/soxs_test_data" >> $HOME/.config/soxs/soxs.cfg
  cat $HOME/.config/soxs/soxs.cfg

fi

# Install dependencies using conda and pip

conda install --yes numpy pytest pip astropy scipy cython h5py tqdm pyyaml appdirs regions

if [[ ${mode} == "wheels" ]]; then
  conda install --yes wheel setuptools
fi

# Install soxs
if [[ ${mode} == "testing" ]]; then
  python -m pip install -e .
fi

# Install pyxsim
if [[ ${mode} == "testing" ]]; then
  git clone https://github.com/jzuhone/pyxsim
  cd pyxsim
  git checkout soxs_read
  pip install .
  cd ..
fi
