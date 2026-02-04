from PIL import Image
import streamlit as st
from streamlit_image_coordinates import streamlit_image_coordinates

import subprocess
import threading
import tempfile
import hashlib
import json
import sqlite3
import shutil
from pathlib import Path

import cv2
from skimage.measure import shannon_entropy

import numpy as np

# =========================
# Thread-safe database access
# =========================

TEMP_IMAGE = Path("./web_images/blackhole_cli.png")
CACHE_DIR = Path("./cache")
CACHE_IMG_DIR = CACHE_DIR / "images"
CACHE_DB = CACHE_DIR / "cache.db"

CACHE_IMG_DIR.mkdir(parents=True, exist_ok=True)

# Global lock for database access
_db_lock = threading.Lock()

def init_cache_db():
    with _db_lock:
        with sqlite3.connect(CACHE_DB, timeout=30.0) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS renders (
                    hash TEXT PRIMARY KEY,
                    image_path TEXT NOT NULL
                )
            """)
            
def cache_lookup(h):
    with _db_lock:
        with sqlite3.connect(CACHE_DB, timeout=30.0, check_same_thread=False) as conn:
            row = conn.execute(
                "SELECT image_path FROM renders WHERE hash = ?",
                (h,)
            ).fetchone()
        return row[0] if row else None

def cache_store(h, image_path):
    with _db_lock:
        with sqlite3.connect(CACHE_DB, timeout=30.0, check_same_thread=False) as conn:
            conn.execute(
                "INSERT OR REPLACE INTO renders (hash, image_path) VALUES (?, ?)",
                (h, image_path)
            )
            conn.commit()

init_cache_db()

def render_hash(params: dict) -> str:
    key = json.dumps(params, sort_keys=True)
    return hashlib.sha256(key.encode()).hexdigest()

# =========================
# Thread-safe rendering
# =========================

def run_render(params, output_path):
    """Run the rendering process with isolated output file"""
    prec_str = "--double" if params["precision"] == "double" else "--float"
    
    # Run main executable
    subprocess.run([
        "./main", 
        "--a", str(params["a"]),
        "--rs", str(params["rs"]), 
        "--Q", str(params["Q"]), 
        "--de0", str(params["de0"]), 
        "--errormax", str(params["errormax"]),
        "--SZELES", str(params["SZELES"]), 
        "--MAGAS", str(params["MAGAS"]), 
        "--kepernyoSZELES", str(params["kepernyoSZELES"]), 
        "--kepernyoMAGAS", str(params["kepernyoMAGAS"]), 
        "--ikezd", str(params["ikezd"]), 
        "--jkezd", str(params["jkezd"]), 
        "--iveg", str(params["iveg"]), 
        prec_str
    ], check=True)
    
    # Run image maker
    subprocess.run(["python", "cli_imagemaker.py"], check=True)
    
    # Copy from temp location to output
    temp_image = Path("./web_images/blackhole_cli.png")
    if temp_image.exists():
        shutil.copy(temp_image, output_path)
    else:
        raise FileNotFoundError(f"Render failed: {temp_image} not found")

def get_or_render_image(params):
    """Thread-safe function to get cached or render new image"""
    h = render_hash(params)
    
    # Check cache first (fast path, no lock needed for read)
    cached_image = cache_lookup(h)
    if cached_image and Path(cached_image).exists():
        return str(cached_image), True
    
    # Not in cache, need to render
    # Use hash as part of temp filename to avoid collisions
    temp_output = CACHE_IMG_DIR / f"{h}_temp.png"
    final_output = CACHE_IMG_DIR / f"{h}.png"
    
    # Check again if another thread just created it
    if final_output.exists():
        cache_store(h, str(final_output))
        return str(final_output), True
    
    try:
        # Render to temp file
        run_render(params, temp_output)
        
        # Move to final location atomically
        shutil.move(str(temp_output), str(final_output))
        
        # Store in cache
        cache_store(h, str(final_output))
        
        return str(final_output), False
        
    except Exception as e:
        # Clean up temp file on error
        if temp_output.exists():
            temp_output.unlink()
        raise e

# =========================
# Streamlit App
# =========================

st.title("Black Hole Ray Tracer: Click Tracker")

if "image_version" not in st.session_state:
    st.session_state.image_version = 0

de0_def = 0.01
errormax_def = 0.001
kepernyoSZELES_def = 9223372036854775808
SZELES = 640
MAGAS = 320

# Initialize session state
if "prec_prev" not in st.session_state:
    st.session_state.prec_prev = False
if "prec_double" not in st.session_state:
    st.session_state.prec_double = False
if "click_count" not in st.session_state:
    st.session_state.click_count = 0
if "kepernyoSZELES" not in st.session_state:
    st.session_state.kepernyoSZELES = kepernyoSZELES_def
if "kepernyoMAGAS" not in st.session_state:
    st.session_state.kepernyoMAGAS = kepernyoSZELES_def//2
if "SZELES" not in st.session_state:
    st.session_state.SZELES = 640
if "MAGAS" not in st.session_state:
    st.session_state.MAGAS = 320
if "ikezd" not in st.session_state:
    st.session_state.ikezd = 0
if "jkezd" not in st.session_state:
    st.session_state.jkezd = 0
if "iveg" not in st.session_state:
    st.session_state.iveg = kepernyoSZELES_def
if "subkepernyoSZELES" not in st.session_state:
    st.session_state.subkepernyoSZELES = kepernyoSZELES_def
if "errormax" not in st.session_state:
    st.session_state.errormax = errormax_def
if "de0" not in st.session_state:
    st.session_state.de0 = de0_def
if "fast" not in st.session_state:
    st.session_state.fast = True 
if "rs" not in st.session_state:
    st.session_state.rs = 0.05
if "a" not in st.session_state:
    st.session_state.a = 0.0
if "fast_spining" not in st.session_state:
    st.session_state.fast_spining = False
if "Q" not in st.session_state:
    st.session_state.Q = 0.0

# Update parameters based on settings
if st.session_state.fast == True:
    st.session_state.de0 = de0_def
    st.session_state.errormax = errormax_def
else:
    st.session_state.de0 = 0.0001
    st.session_state.errormax = 0.000001

if st.session_state.fast_spining == False:
    st.session_state.a = 0.0
else:
    st.session_state.a = st.session_state.rs/2 - 0.001

# Build render parameters
render_params = {
    "SZELES": st.session_state.SZELES,
    "MAGAS": st.session_state.MAGAS,
    "kepernyoSZELES": st.session_state.kepernyoSZELES,
    "kepernyoMAGAS": st.session_state.kepernyoMAGAS,
    "ikezd": st.session_state.ikezd,
    "jkezd": st.session_state.jkezd,
    "iveg": st.session_state.iveg,
    "precision": "double" if st.session_state.prec_double else "float",
    "de0": st.session_state.de0,
    "errormax": st.session_state.errormax,
    "rs": st.session_state.rs,
    "a": st.session_state.a,
    "Q": st.session_state.Q
}

# Get or render image (thread-safe)
try:
    IMAGE_PATH, was_cached = get_or_render_image(render_params)
    
    if was_cached:
        st.info("📦 Cache hit – image loaded from disk")
    else:
        st.success("✨ Newly rendered image")
        
    # Display image
    img = Image.open(IMAGE_PATH)
    st.success(str(IMAGE_PATH))
    value = streamlit_image_coordinates(img, key=f"image_{st.session_state.image_version}")

    # Calculate image metrics
    img_cv = cv2.imread(IMAGE_PATH)
    gray = cv2.cvtColor(img_cv, cv2.COLOR_RGB2GRAY)
    entropy_value = shannon_entropy(gray)
    variance = float(np.var(gray)) * 0.2 / 255.0
    edges = cv2.Canny(gray, 50, 150)
    edge_density = edges.mean()
    
    st.success(f"Shannon entropy: {entropy_value:.3f}, Variance: {variance:.3f}, Edge density: {edge_density:.3f}")
    st.success(f"ikezd={st.session_state.ikezd}, jkezd={st.session_state.jkezd}, iveg={st.session_state.iveg}")

except Exception as e:
    st.error(f"Error rendering image: {e}")
    st.stop()

# Controls
st.checkbox("🔬 Use double precision", key="prec_double")
st.checkbox("⚡ Fast mode — reduces runtime by using larger steps ⚡", key="fast")
st.checkbox("🌀 Critically spinning black hole", key="fast_spining")

if st.button("🔄 Reset view"):
    st.session_state.iveg = kepernyoSZELES_def
    st.session_state.ikezd = 0
    st.session_state.jkezd = 0
    st.session_state.subkepernyoSZELES = kepernyoSZELES_def
    st.session_state.image_version += 1
    st.success("View reset")
    st.rerun()

if st.button("🧮 Double vs Float view"):
    st.session_state.iveg = 451518464
    st.session_state.ikezd = 451517952
    st.session_state.jkezd = 206130944
    st.success("View changed")
    st.rerun()

# Handle clicks
if value:
    click_x = value["x"]
    click_y = value["y"]
    
    if click_x < SZELES//3:
        st.session_state.iveg = st.session_state.ikezd + st.session_state.subkepernyoSZELES//2
    elif click_x < 2*SZELES//3:
        st.session_state.iveg = st.session_state.ikezd + 3*st.session_state.subkepernyoSZELES//4
        st.session_state.ikezd = st.session_state.ikezd + st.session_state.subkepernyoSZELES//4
    else:
        st.session_state.iveg = st.session_state.ikezd + st.session_state.subkepernyoSZELES
        st.session_state.ikezd = st.session_state.ikezd + st.session_state.subkepernyoSZELES//2
        
    if click_y < MAGAS//3:
        pass
    elif click_y < 2*MAGAS//3:
        st.session_state.jkezd = st.session_state.jkezd + st.session_state.subkepernyoSZELES//8
    else:
        st.session_state.jkezd = st.session_state.jkezd + st.session_state.subkepernyoSZELES//4
        
    st.session_state.subkepernyoSZELES = st.session_state.subkepernyoSZELES//2
    st.success(f"User clicked at: X={click_x}, Y={click_y}")
    st.session_state.image_version += 1
    st.rerun()
else:
    st.info("Click anywhere on the image to track coordinates.")
