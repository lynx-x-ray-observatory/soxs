#!/bin/sh

DEST_DIR=${HOME}/xrs_tools

PY_VERSION=2 # 2 installs Python 2.7, 3 installs Python 3.5.

INST_SIMX=0 # 1 means we install SIMX, 0 means we don't.
INST_PYXSIM=0 # 1 means we install pyXSIM, 0 means we don't.
INST_GIT=1 # 1 means we install git, 0 means we don't.

###############################################################
# You should not need to modify anything below this line.
###############################################################

echo "Beginning installation of xrs_tools."

MYARCH=`uname -m`       # A guess at the OS
MYOS=`uname -s`

MINICONDA_URLBASE="http://repo.continuum.io/miniconda"
MINICONDA_VERSION="latest"

if [ ! -e ${DEST_DIR} ]
then
    mkdir $DEST_DIR
fi

ORIG_DIR=`pwd`

cd $DEST_DIR

if type -P curl &>/dev/null
then
    echo "Using curl to download files."
    export GETFILE="curl -OL"
else
    echo "Using wget to download files."
    export GETFILE="wget"
fi

if [ $MYOS = "Darwin" ]
then
    echo "Looks like you're running Mac OS."
    MINICONDA_OS="MacOSX"
    MINICONDA_ARCH="x86_64"
elif [ $MYOS = "Linux" ]
then
    echo "Looks like you're running Linux."
    MINICONDA_OS="Linux"
    if [ $MYARCH = "i386" ]
    then
        MINICONDA_ARCH="x86"
    elif [ $MYARCH = "i686"  ]
    then
        MINICONDA_ARCH="x86"
    elif [ $MYARCH = "x86_64"  ]
    then
        MINICONDA_ARCH="x86_64"
    else
        echo "Not sure which architecture you are running."
        echo "Going with x86_64 architecture."
        MINICONDA_OS="Linux-x86_64"
    fi
fi

MINICONDA_PKG="Miniconda${PY_VERSION}-${MINICONDA_VERSION}-${MINICONDA_OS}-${MINICONDA_ARCH}.sh"

if [ ! -e $DEST_DIR/anaconda ]
then
    if [ -f ${MINICONDA_PKG} ]
    then
        rm ${MINICONDA_PKG}
    fi
    echo "Downloading ${MINICONDA_URLBASE}/${MINICONDA_PKG}"
    ${GETFILE} ${MINICONDA_URLBASE}/${MINICONDA_PKG}
fi

if [ ! -f $DEST_DIR/anaconda_done ]
then
    echo "Installing the Miniconda python environment."
    bash ./${MINICONDA_PKG} -b -p $DEST_DIR -f
    touch $DEST_DIR/anaconda_done
fi

export PATH=${DEST_DIR}/bin:$PATH

if [ $INST_GIT -eq 1 ]
then
    conda install --yes git
fi

if [ $INST_PYXSIM -eq 1 ]
then
    conda install --yes -c jzuhone pyxsim
fi

# Install SIMX

SIMX="simx-2.4.1"

if [ $INST_SIMX -eq 1 ]
then
    if [ ! -e $SIMX ]
    then
        echo "Downloading SIMX..."
        ${GETFILE} http://hea-www.harvard.edu/simx/$SIMX.tar.gz
        tar xfz $SIMX.tar.gz
    fi
    if [ ! -e $SIMX/done ]
    then
        echo "Installing SIMX..."
        cd $SIMX
        ( ./configure --prefix=`pwd` 2>&1 ) 
        ( make install 2>&1 )
        touch done
        cd ..
    fi
fi

ACTIVATE_SH=${DEST_DIR}/bin/activate.sh

if [ -f $ACTIVATE_SH ]
then
    rm $ACTIVATE_SH
fi
echo '#!/bin/sh' >> $ACTIVATE_SH
echo 'PATH='$DEST_DIR'/bin:$PATH' >> $ACTIVATE_SH
if [ $INST_SIMX -eq 1 ]
then
    echo '. '$DEST_DIR/$SIMX'/bin/simx-setup.sh' >> $ACTIVATE_SH
fi

ACTIVATE_CSH=${DEST_DIR}/bin/activate.csh

if [ -f $ACTIVATE_CSH ]
then
    rm $ACTIVATE_CSH
fi
echo '#!/bin/csh' >> $ACTIVATE_CSH
echo 'set path = ('$DEST_DIR'/bin $path)' >> $ACTIVATE_CSH
if [ $INST_SIMX -eq 1 ]
then
    echo 'source '$DEST_DIR/$SIMX'/bin/simx-setup.csh' >> $ACTIVATE_CSH
fi

cd $ORIG_DIR

echo "All done!"
