# Powershell script

# Install conda environment
conda config --set always_yes yes
conda create --name desyre c-compiler cxx-compiler cmake fmt boost boost-cpp symengine python=3.10
conda activate desyre

# Define environment variables for clang ...
# ... and make them persistent 
cmd.exe /c "call `"C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat`" && set > %temp%\vcvars.txt"

Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
  if ($_ -match "^(.*?)=(.*)$") {
    Set-Content "env:\$($matches[1])" $matches[2]
  }
}

#build and run the desyre tests
mkdir build
cd build
cmake `
    -G "Ninja" `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\desyre `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\desyre `
    -DDSYRE_BUILD_TESTS=yes ..

cmake --build . --target install --config Release
ctest -j4 -V -C Release

