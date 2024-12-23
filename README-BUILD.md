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
sudo yum install -y package/rpm/*/*.rpm
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
( cd package/deb; sudo apt install ./\*.deb )
```

Finally, clean up:

```
make distclean
```

## MacOS Brew

If you haven't already installed "brew", use the following steps to
install it (see https://phoenixnap.com/kb/install-homebrew-on-mac):

- `xcode-select --install`
- Select Install and accept the license
- `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
- Enter admin password when prompted
- Press return to install Homebrew


### Preparing the Source From a .tar.gz file (tarball)

```
brew install libyaml
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

To install using brew, first update your brew environment:

```
brew update
```

Next, create/update the Formula for libvna.  This will put you into
your editor (default is vi).

```
brew create https://github.org/scott-guthridge/libvna/archive/refs/tags/libvna-X.Y.Z.tar.gz
```

Make the file contents look like the following (with version numbers
and SHA256 sum updated apporopriately):

```
class Libvna < Formula
  desc "Vector Network Analyzer Library"
  homepage "https://github.com/scott-guthridge/libvna"
  url "https://github.com/scott-guthridge/libvna/archive/refs/tags/libvna-0.3.6.tar.gz"
  sha256 "d88c9fc74612f460755fc62eafd0165a1016b5dda5518f3b3fc9a31dcf4bbbee"
  license "GPL-3"
  depends_on "libyaml"

  def install
    system "./configure", *std_configure_args
    system "make", "-j12", "check"
    system "make", "install"
  end
end
```

Install:

```
brew install libvna
```

---

# Creating PDF Versions of the Man Pages

```
make pdfman
```

Generated files appear in src/\*.pdf
