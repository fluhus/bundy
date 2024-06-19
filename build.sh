# Compiles bundy executables for all platforms.
set -e

rm -fr build
mkdir build

# Linux
go build -o build ./bundy ./bundyx
zip -j build/bundy_linux_amd64.zip build/bundy build/bundyx

# Mac
GOOS=darwin GOARCH=arm64 go build -o build ./bundy ./bundyx
zip -j build/bundy_macos_arm64.zip build/bundy build/bundyx

# Windows
GOOS=windows go build -o build ./bundy ./bundyx
zip -j build/bundy_win_amd64.zip build/bundy.exe build/bundyx.exe

rm build/bundy build/bundyx build/bundy.exe build/bundyx.exe
