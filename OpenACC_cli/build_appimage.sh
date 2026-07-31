#!/bin/bash
set -e

if [ ! -f main ] || [ ! -f libblackhole.so ]; then
    echo "main and/or libblackhole.so not found - run 'make' first." >&2
    exit 1
fi

echo "Cleaning up old build..."
rm -rf AppDir
mkdir -p AppDir/usr/lib AppDir/src AppDir/usr/python

echo "Copying application files..."
cp main libblackhole.so *.py kep.dat icon_black_hole.bmp AppDir/src/
cp -r static AppDir/src/
# if there is a web_images dir, copy it too
if [ -d "web_images" ]; then
    cp -r web_images AppDir/src/
fi
if [ -d "cache" ]; then
    cp -r cache AppDir/src/
fi

echo "Extracting libraries from nvcc container..."
distrobox-enter nvcc -e bash -c 'cp /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libacchost.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libaccdevaux.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libaccdevice.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libnvhpcatm.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libnvomp.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libnvcpumath.so /opt/nvidia/hpc_sdk/Linux_x86_64/26.3/compilers/lib/libnvc.so AppDir/usr/lib/'

echo "Downloading portable Python 3.11..."
wget -q -O python.tar.gz "https://github.com/astral-sh/python-build-standalone/releases/download/20260508/cpython-3.11.15%2B20260508-x86_64-unknown-linux-gnu-install_only.tar.gz"
tar -xzf python.tar.gz -C AppDir/usr/python --strip-components=1
rm python.tar.gz

echo "Installing Python dependencies..."
AppDir/usr/python/bin/python -m pip install fastapi uvicorn pillow streamlit numpy scikit-image filelock opencv-python-headless streamlit-image-coordinates

echo "Creating AppRun..."
cat << 'EOF' > AppDir/AppRun
#!/bin/bash
HERE="$(dirname "$(readlink -f "${0}")")"
export LD_LIBRARY_PATH="$HERE/usr/lib:$LD_LIBRARY_PATH"
export PATH="$HERE/usr/python/bin:$PATH"

# Create a writable working directory in user's cache
WORKDIR="${XDG_CACHE_HOME:-$HOME/.cache}/OpenACC_cli"
mkdir -p "$WORKDIR/web_images"
mkdir -p "$WORKDIR/cache/images"

cd "$WORKDIR"

# Symlink all read-only files and dirs from the AppImage into the writable
# WORKDIR (e.g. static/). web_images/ and cache/ already exist as real,
# writable dirs from the mkdir -p above, so the existence check leaves them
# alone instead of shadowing them with a read-only symlink.
for item in "$HERE/src"/*; do
    base=$(basename "$item")
    if [ ! -e "$base" ]; then
        ln -sf "$item" "./$base"
    fi
done

PORT=8000
( sleep 1; xdg-open "http://127.0.0.1:$PORT" >/dev/null 2>&1 || true ) &
exec "$HERE/usr/python/bin/uvicorn" web_server:app --host 127.0.0.1 --port "$PORT" "$@"
EOF
chmod +x AppDir/AppRun

echo "Creating Desktop file..."
cat << 'EOF' > AppDir/OpenACC_cli.desktop
[Desktop Entry]
Name=OpenACC_cli
Exec=AppRun
Icon=icon_black_hole
Type=Application
Categories=Science;
EOF

echo "Copying icon..."
# appimagetool requires a .png or .svg in the root of AppDir matching the Icon name
# Let's convert bmp to png using imagemagick
distrobox-enter nvcc -e bash -c 'apt-get update && apt-get install -y imagemagick' || true
distrobox-enter nvcc -e convert icon_black_hole.bmp AppDir/icon_black_hole.png || cp icon_black_hole.bmp AppDir/icon_black_hole.png # fallback

echo "Downloading appimagetool..."
wget -q -O appimagetool-x86_64.AppImage "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
chmod +x appimagetool-x86_64.AppImage

echo "Building AppImage..."
# Extract appimagetool to avoid FUSE issues in some environments
./appimagetool-x86_64.AppImage --appimage-extract
./squashfs-root/AppRun AppDir OpenACC_cli-x86_64.AppImage

echo "Done! You can run ./OpenACC_cli-x86_64.AppImage"
