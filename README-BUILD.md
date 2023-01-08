# VNA Library Build Instructions

The following sections give detailed build instructions on various platforms.

## Fedora, RedHat, Suse, etc. Linux (RPM-based distros)

### Preparing the Source From a .tar.gz file (tarball)

```
sudo yum install -y gcc libyaml-devel make man m4
tar xzf libvna-X.Y.Z.tar.gz
cd libvna-X.Y.Z
./configure
```

### Alternate Method: Preparing the Source From Git

```
sudo yum install -y autoconf automake gcc git libtool libyaml-devel make man m4
git clone https://github.com/scott-guthridge/libvna.git
cd libvna
./bootstrap
./configure
```

### Local Install

Local install installs typically in /usr/local/\*

```
make -j12 check
sudo make install
make distclean
```

### Alternative Method: RPM Install

RPM install has the advantage that binaries are installed in the standard
system directories: /usr/bin, /usr/lib\*, /usr/share.  RPM also keeps track
of the install version and locations of all installed files.

First, build the RPMs:

```
sudo yum install -y rpm-build
make rpm
```

Next, install:

```
sudo yum install -y rpm/*/*.rpm
```

Finally, clean up:

```
make distclean
```

## Ubuntu, Debian, Arch, etc. Linux (Debian-based distros)

### Preparing the Source From a .tar.gz file (tarball)

```
sudo apt-get update
sudo apt-get install -y gcc libyaml-dev make man m4
tar xzf libvna-X.Y.Z.tar.gz
cd libvna-X.Y.Z
./configure
```

### Alternate Method: Preparing the Source From Git

```
sudo apt-get update
sudo apt-get install -y autoconf automake gcc git libtool libyaml-dev \
	make man m4
git clone https://github.com/scott-guthridge/libvna.git
cd libvna
./bootstrap
./configure
```

### Local Install

Local install installs typically in /usr/local/\*

```
make -j12 check
sudo make install
make distclean
```

### Alternative Method: APT Install

APT install has the advantage that binaries are installed in the standard
system directories: /usr/bin, /usr/lib\*, /usr/share.  APT also keeps track
of the install version and locations of all installed files.

First, build the packages:

```
sudo apt-get update
sudo apt-get install -y build-essential devscripts debhelper-compat
make deb
```

Next, install:

```
( cd deb; sudo apt install ./\*.deb )
```

Finally, clean up:

```
make distclean
```

## MacOS Brew

### Preparing the Source From a .tar.gz file (tarball)

```
brew install autoconf automake libtool libyaml
tar xzf libvna-X.Y.Z.tar.gz
cd libvna-X.Y.Z
./configure
```

### Alternate Method: Preparing the Source From Git

```
brew install autoconf automake libyaml
git clone https://github.com/scott-guthridge/libvna.git
cd libvna
./bootstrap
./configure
```

### Local Install

Local install installs in /usr/local/\*

```
make -j12 check
sudo make install
make distclean
```

### Alternative Method: Brew Install

TODO: fill in details

---

# Creating PDF Versions of the Man Pages

```
make pdfman
```

Generated files appear in src/\*.pdf
