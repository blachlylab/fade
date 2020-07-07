# Installation

Install [htslib](http://www.htslib.org/download/) and [parasail](https://github.com/jeffdaily/parasail#compiling-and-installing) and download a binary fade release.
<details>
<summary> Prerequisites </summary>

#### Install htslib

Install [htslib](http://www.htslib.org/download/) prerequisites.
```
Debian / Ubuntu
---------------

sudo apt-get update  # Ensure the package list is up to date
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

Note: libcurl4-openssl-dev can be used as an alternative to libcurl4-gnutls-dev.

RedHat / CentOS
---------------

sudo yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel

Alpine Linux
------------

sudo apk update  # Ensure the package list is up to date
sudo apk add autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev

OpenSUSE
--------

sudo zypper install autoconf automake make gcc perl zlib-devel libbz2-devel xz-devel libcurl-devel libopenssl-devel
```
In addition, please make sure ```wget``` and ```bzip2``` are installed in order to follow the rest of the instructions.
Now download and make htslib. This specifically must be version [1.9](https://github.com/samtools/htslib/releases/tag/1.9).
```
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -xjf htslib-1.9.tar.bz2
cd htslib-1.9
./configure
make 
sudo make install
```
#### Install parasail

Download and install parasail.
```
wget https://github.com/jeffdaily/parasail/releases/download/v2.4.2/parasail-2.4.2-manylinux1_x86_64.tar.gz
tar -xzf parasail-2.4.2-manylinux1_x86_64.tar.gz
cd parasail-2.4.2-manylinux1_x86_64
cd lib/
sudo cp * /usr/local/lib/
```
#### Other
Make sure ```/usr/local/lib``` is on your ```LD_LIBRARY_PATH```.
```
# add to your bashrc 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```
</details>
<details>
<summary> Install FADE </summary>

Go to the [latest release](https://github.com/blachlylab/fade/releases/latest).
```
wget https://github.com/blachlylab/fade/releases/download/v0.1.1/fade
sudo cp fade /usr/local/bin
```
</details>
<details>
<summary> Building FADE from source </summary>

Build from source using dub and a D compiler. We recommend ldc2 as the compiler, however dmd should work as well. For more information on D compilers, visit [here](https://dlang.org/download.html).
#### Install dub and ldc2 (Preferred)
```
curl -fsS https://dlang.org/install.sh | bash -s ldc
```
#### Install dub and dmd
```
curl -fsS https://dlang.org/install.sh | bash -s dmd
```
#### Build FADE
```
source ~/dlang/*compiler*/activate
git clone https://github.com/blachlylab/fade.git 
cd fade
LIBRARY_PATH=/usr/local/lib/ dub build -b release

# deactivate dlang environment
deactivate
```
</details>
<details>
<summary> Non-root installation </summary>

Non-root installs can be tricky. It basically sums up to making sure all necessary shared libraries can be found by the fade executable.
An easy way of install htslib in a non-root capacity is via [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Install htslib.
```
conda install -c bioconda htslib
```
Download and install parasail.
```
wget https://github.com/jeffdaily/parasail/releases/download/v2.4.2/parasail-2.4.2-manylinux1_x86_64.tar.gz
tar -xzf parasail-2.4.2-manylinux1_x86_64.tar.gz
cd parasail-2.4.2-manylinux1_x86_64
cd lib/
cp * minconda3/lib/
```
Make sure ```miniconda3/lib``` is on your ```LD_LIBRARY_PATH```.
```
# add to your bashrc 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```

Now follow previous instructions to install FADE. If compiling change this line.
```
LIBRARY_PATH=~/miniconda/lib/ dub build -b release
```
</details>

<details>
<summary> Fade installation with multiple htslib versions</summary>

FADE requires htslib version 1.9, though the latest version is 1.10.2 (as of writing).
Your htslib version may be more up to date than the one our instructions would have you install.

In the case of using a provided binary, this should have no effect. FADE's binary will be able to 
find the correct shared library for htslib.

However, if building fade from source you should ensure that your htslib.so symbolic link 
under /usr/local/lib points to htslib.so.2 .

</details>
