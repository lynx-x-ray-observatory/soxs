set -x   # Show which command is being run

# Download answers

wget -q http://hea-www.cfa.harvard.edu/~jzuhone/soxs_test_data_${ANSWER_VER}.tar.gz
tar -zxf soxs_test_data_${ANSWER_VER}.tar.gz

# Set location of soxs data

mkdir -p $HOME/.config/soxs
echo "[soxs]" > $HOME/.config/soxs/soxs.cfg
echo "soxs_data_dir = \"${PWD}/soxs_data\"" >> $HOME/.config/soxs/soxs.cfg
cat $HOME/.config/soxs/soxs.cfg

# Install dependencies using conda

conda install --yes numpy pytest pip astropy scipy cython h5py tqdm pyyaml regions

# Install soxs

python -m pip install -e .
