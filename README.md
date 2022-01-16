# a.multiminer

Fork Preface
------------

a.multiminer is a fork of the multiminer created by bogdanadnan 
with added support for the algorithm Argon2d16000 used by Alterdot 
and without the dev fee as no further development has been done by the initial developer for over 2 years.

The rest of this README is left mostly untouched, as it was left by its creator.

Description
------------

multiminer is a fork of cpuminer-opt by Jay D Dee. (https://github.com/JayDDee/cpuminer-opt)

I've added support for GPU mining a couple of Argon2D coins:
Zumy/Credits (argon2d250), Dynamic (argon2d500) and Argentum/Unitus (argon2d4096).
Support for UraniumX (argon2ad) was also added and is in beta testing, with a lot of
hashrate optimizations for other algos. (**please checkout beta branch to use it - 
keep in mind that it might not be stable tough**)

Miner programs are often flagged as malware by antivirus programs. This is
a false positive, they are flagged simply because they are cryptocurrency 
miners. The source code is open for anyone to inspect. If you don't trust 
the software, don't use it.

Requirements
------------

1. A x86_64 architecture CPU with a minimum of SSE2 support. This includes
Intel Core2 and newer and AMD equivalents. In order to take advantage of AES_NI
optimizations a CPU with AES_NI is required. This includes Intel Westbridge
and newer and AMD equivalents. Further optimizations are available on some
algoritms for CPUs with AVX and AVX2, Sandybridge and Haswell respectively.

Older CPUs are supported by cpuminer-multi by TPruvot but at reduced
performance.

ARM CPUs are not supported.

2. 64 bit Linux OS. Ubuntu and Fedora based distributions, including Mint and
Centos, are known to work and have all dependencies in their repositories.
Others may work but may require more effort. Older versions such as Centos 6
don't work due to missing features. 
64 bit Windows OS is supported with mingw_w64 and msys or pre-built binaries.

MacOS, OSx and Android are not supported.

3. Stratum pool. Some algos may work wallet mining using getwork or GBT. YMMV.

Windows Building Process
---------------------
The following instructions have been supplied by frozen80 from the Alterdot community and they can be used to compile a.multiminer for Windows with OpenCL only. It will work for both Nvidia and AMD cards. CUDA support can be enabled on Linux for now.

1. Download MSYS2 from here https://www.msys2.org/ and install the necessary tools. 
Follow the tutorial on the page to update the system packages.
    
    First step:
    ```sh
    pacman -Syu
    ```
    Second step(after restart of msys):
    ```sh
    pacman -Su
    ```
    Third step: Install some tools for compilation:
    ```sh
    pacman -S --needed base-devel mingw-w64-x86_64-toolchain
    pacman -S mingw-w64-x86_64-cmake
    pacman -S mingw-w64-x86_64-opencl-headers
    pacman -S mingw-w64-x86_64-opencl-icd
    ```

2. Open "MSYS2 MinGW x64" shortcut and check if you have all the tools properly installed
    ```sh
    gcc -v
    g++ -v
    cmake --version
    ```
   Now we need to create a symlink to mingw32-make.exe:
    ```sh
    cd /mingw64/bin/ && ln mingw32-make.exe make.exe
    make -v
    ```

3. Go to Windows environment variables and the msys64/mingw64/bin folder (usually C:/msys64/mingw64/bin) to your path.

4. Download the a.multiminer source code and unzip it. Create a folder named build inside. Open a command prompt and navigate to the build folder. Before we can compile we need to make a few changes to the source code:
    
    Edit file CMakeLists.txt.

    Uncomment these lines after line 4 (delete the #):
    ```sh
    set(DEPRECATION_FLAGS "-Xlinker --allow-multiple-definition -Wno-error=deprecated")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEPRECATION_FLAGS}")
    set(C_DEP_FLAGS "-Wno-pointer-sign -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEP_FLAGS}")
    set(NO_CUDA TRUE)
    ```

    Go to line 313 and change:
    ```sh
    target_compile_options(a.multiminer PRIVATE -O3 -mtune=native -march=native -I.)
    ```
    to 
    ```sh
    target_compile_options(a.multiminer PRIVATE -O1 -mtune=native -march=native -I.)
    ```
    Now go back to command prompt and make sure you are in the build folder created earlier.
    ```sh
    cmake .. -G "MinGW Makefiles"
    make -j 4
    ```

4. Your a.multiminer.exe file should be available in the build folder.

Linux Building Process
---------------------

These build instructions are for Ubuntu 20.04, 18.04 and 16.04. For any other distribution you might
need to adapt them accordingly (especially CUDA installation). Also keep in mind CUDA
toolkit installed in the following manner might change your drivers. The miner should
compile ok without CUDA, case in which it will use only OpenCL. Just skip any CUDA mention
from following instructions if you don't need it.

1. Install required dependencies:
```sh
sudo apt-get install git cmake gcc g++ libjansson-dev libcurl4-openssl-dev libssl-dev libgmp-dev ocl-icd-opencl-dev  
```
2. Install CUDA toolkit (vers. 9.x or newer). This process differes on different Ubuntu versions.

For Ubuntu 16.04, default nvidia toolkit is too old so please follow the instructions from 
Nvidia site to get a newer version: 
https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

For Ubuntu 18.04, you can use Nvidia instructions or install it from ubuntu repository using this:
```sh
sudo apt-get install gcc-6 g++-6 nvidia-cuda-toolkit  
```
If installed from Ubuntu repository, CUDA version in Ubuntu 18.04 is 9.1. CUDA version 9.x works only 
with gcc/g++ 6.x, while default compiler version is 7. You can check the version you are running 
using this command:
```sh
gcc --version
```

For Ubuntu 20.04, you can use CUDA 10.x.

If compiler version is different than 6 than do this (next commands suppose you have 7 as default version, 
change accordingly if this is not the case):
```sh
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 10
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 20

sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-6 10
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 20

sudo update-alternatives --config gcc
sudo update-alternatives --config g++
```
3. Clone the repository:
```sh
git clone http://github.com/Alterdot/a.multiminer
```
4. Build the source code:
```sh
cd a.multiminer
mkdir build
cd build
cmake ..
make
```
5. You should now have a binary called a.multiminer in current folder.

Usage
-----
This miner has the same general options as cpuminer-opt by Jay D Dee.
Please check regular usage on its GitHub repository:
https://github.com/JayDDee/cpuminer-opt

In order to enable gpu mining there are additional options:

```sh
          --use-gpu=CUDA|OPENCL Use GPU for algorithms supporting it (Argon2d)
          --gpu-id=N        Use GPU device with specific index in detected devices - default 1
          --gpu-batchsize=N Specify batch size - default 1
```

The regular threads option (-t) when used in combination with --use-gpu, will control 
number of threads started on each card you have. You will have to play with threads and
batchsize numbers to get best hashrate. Please note that batchsize is required to get a
normal speed. It should be given in powers of 2 (though that is not mandatory, but is better).
The batchsize represents the number of hashes calculated by the card in a call. Threads multiplied
by batchsize multiplied by argon2d memory requirements (16000 KiB for Alterdot, 250KiB for Zumy, 500KiB for Dynamic
and 4096KiB for Argentum/Unitus) should not exceed available card memory.
Sample command line:

```sh
./a.multiminer -a argon2d16000 -o stratum+tcp://pooladdress:port -u walletaddress -p c=ADOT,workername --use-gpu CUDA -t 2 --gpu-batchsize 64
```

When using OpenCL keep in mind that the program will take a LOT of time to start
hashing (from tens of seconds to minutes, depending on how many threads you want to start).
Please be patient :). 
Also, if you have Nvidia cards, always use CUDA, it gives best performance.

Supported Algorithms
--------------------

                          allium        Garlicoin
                          anime         Animecoin
                          argon2        Argon2 coin (AR2)
                          argon2d250    argon2d-crds, Credits (CRDS) & Zumy (ZMY)
                          argon2d500    argon2d-dyn,  Dynamic (DYN)
                          argon2d4096   argon2d-uis, Unitus, (UIS)
                          argon2d16000  argon2d-adot, Alterdot (ADOT)\n\
                          argon2ad      UraniumX, (URX) - in beta
                          axiom         Shabal-256 MemoHash
                          bastion
                          blake         Blake-256 (SFR)
                          blakecoin     blake256r8
                          blake2s       Blake-2 S
                          bmw           BMW 256
                          c11           Chaincoin
                          decred
                          deep          Deepcoin (DCN)
                          dmd-gr        Diamond-Groestl
                          drop          Dropcoin
                          fresh         Fresh
                          groestl       Groestl coin
                          heavy         Heavy
                          hmq1725       Espers
                          hodl          Hodlcoin
                          jha           Jackpotcoin
                          keccak        Maxcoin
                          keccakc       Creative coin
                          lbry          LBC, LBRY Credits
                          luffa         Luffa
                          lyra2h        Hppcoin
                          lyra2re       lyra2
                          lyra2rev2     lyra2v2, Vertcoin
                          lyra2z        Zcoin (XZC)
                          lyra2z330     Lyra2 330 rows, Zoin (ZOI)
                          m7m           Magi (XMG)
                          myr-gr        Myriad-Groestl
                          neoscrypt     NeoScrypt(128, 2, 1)
                          nist5         Nist5
                          pentablake    Pentablake
                          phi1612       phi, LUX coin
                          pluck         Pluck:128 (Supcoin)
                          polytimos     Ninja
                          quark         Quark
                          qubit         Qubit
                          scrypt        scrypt(1024, 1, 1) (default)
                          scrypt:N      scrypt(N, 1, 1)
                          scryptjane:nf
                          sha256d       Double SHA-256
                          sha256t       Triple SHA-256, Onecoin (OC)
                          shavite3      Shavite3
                          skein         Skein+Sha (Skeincoin)
                          skein2        Double Skein (Woodcoin)
                          skunk         Signatum (SIGT)
                          timetravel    Machinecoin (MAC)
                          timetravel10  Bitcore
                          tribus        Denarius (DNR)
                          vanilla       blake256r8vnl (VCash)
                          veltor        (VLT)
                          whirlpool
                          whirlpoolx
                          x11           Dash
                          x11evo        Revolvercoin
                          x11gost       sib (SibCoin)
                          x12           Galaxie Cash (GCH)
                          x13           X13
                          x13sm3        hsr (Hshare)
                          x14           X14
                          x15           X15
                          x16r          Ravencoin (RVN)
                          x16s          pigeoncoin (PGN)
                          x17
                          xevan         Bitsend (BSD)
                          yescrypt      Globalboost-Y (BSTY)
                          yescryptr8    BitZeny (ZNY)
                          yescryptr16   Yenten (YTN)
                          yescryptr32   WAVI
                          zr5           Ziftr

Errata
------

Neoscrypt crashes on Windows, use legacy version.

AMD CPUs older than Piledriver, including Athlon x2 and Phenom II x4, are not
supported by cpuminer-opt due to an incompatible implementation of SSE2 on
these CPUs. Some algos may crash the miner with an invalid instruction.
Users are recommended to use an unoptimized miner such as cpuminer-multi.

cpuminer-opt does not work mining Decred algo at Nicehash and produces
only "invalid extranonce2 size" rejects.

Benchmark testing does not work for x11evo.

Bugs
---

Please create an issue on GitHub.
All problem reports must be accompanied by a proper definition.
This should include how the problem occurred, the command line and
output from the miner showing the startup and any errors.

Happy mining!
