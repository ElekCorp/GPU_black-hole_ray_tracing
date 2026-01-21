from PIL import Image

import subprocess

import hashlib
import json
import sqlite3
import shutil
from pathlib import Path

from itertools import product
import cv2
from skimage.measure import shannon_entropy

import numpy as np

def is_black(path):#if it's not interesting or totaly black we block further zooming
    if path and Path(path).exists():
        img_cv = cv2.imread(path)
        gray = cv2.cvtColor(img_cv, cv2.COLOR_RGB2GRAY)
        #entropy_value = shannon_entropy(gray)
        #variance = float(np.var(gray))
        edges = cv2.Canny(gray, 50, 150)
        edge_density = edges.mean()
        print(edge_density)
        if edge_density <4.13:
            return True
        #if edge_density>2.5:
        #    return True
    else:
        return False

# =========================
# Paths
# =========================

TEMP_IMAGE = Path("./web_images/blackhole_cli.png")
CACHE_DIR = Path("./cache")
CACHE_IMG_DIR = CACHE_DIR / "images"
CACHE_DB = CACHE_DIR / "cache.db"

CACHE_IMG_DIR.mkdir(parents=True, exist_ok=True)

def init_cache_db():
    with sqlite3.connect(CACHE_DB) as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS renders (
                hash TEXT PRIMARY KEY,
                image_path TEXT NOT NULL
            )
        """)
def cache_lookup(h):
    with sqlite3.connect(CACHE_DB) as conn:
        row = conn.execute(
            "SELECT image_path FROM renders WHERE hash = ?",
            (h,)
        ).fetchone()
    return row[0] if row else None

def cache_store(h, image_path):
    with sqlite3.connect(CACHE_DB) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO renders (hash, image_path) VALUES (?, ?)",
            (h, image_path)
        )

init_cache_db()

def render_hash(params: dict) -> str:
    key = json.dumps(params, sort_keys=True)
    return hashlib.sha256(key.encode()).hexdigest()




de0_def=0.01
errormax_def=0.001
kepernyoSZELES_def = 2147483648*8
SZELES=640
MAGAS=320

def tile_number_to_click(x_in): #return click_x,click_y
#0 1 2
#3 4 5
#6 7 8
    x=x_in%9
    if x==0:
        return 1,1
    if x==1:
        return SZELES//2,1
    if x==2:
        return SZELES-1,1
    if x==3:
        return 1,MAGAS//2
    if x==4:
        return SZELES//2,MAGAS//2
    if x==5:
        return SZELES-1,MAGAS//2
    if x==6:
        return 1,MAGAS-1
    if x==7:
        return SZELES//2,MAGAS-1
    if x==8:
        return SZELES-1,MAGAS-1


def check_path_is_black(partial_path, fast_spining):
    """Check if a partial path leads to a black image"""
    prec_prev = False
    prec_double = False
    click_count = 0
    kepernyoSZELES = kepernyoSZELES_def
    kepernyoMAGAS  = kepernyoSZELES_def//2
    SZELES = 640
    MAGAS  = 320
    ikezd=0
    jkezd=0
    iveg=kepernyoSZELES_def
    subkepernyoSZELES=kepernyoSZELES_def
    errormax = errormax_def
    de0 = de0_def
    fast = True 
    rs=0.05
    a = 0.0
    Q=0.0
    
    if prec_double == True:
        prec_str = "--double"
    else:
        prec_str = "--float"
    if fast == True:
        de0=de0_def
        errormax=errormax_def
    else:
        de0=0.0001
        errormax=0.000001
    if fast_spining==False:
        a=0.0
    else:
        a=rs/2-0.001
    
    for k in partial_path:
        click_x, click_y = tile_number_to_click(k)
        if click_x < SZELES//3:
            iveg=ikezd+subkepernyoSZELES//2
        elif click_x < 2*SZELES//3:
            iveg=ikezd+3*subkepernyoSZELES//4
            ikezd=ikezd+subkepernyoSZELES//4
        else:
            iveg=ikezd+subkepernyoSZELES
            ikezd=ikezd+subkepernyoSZELES//2
        if click_y < MAGAS//3:
            None
        elif click_y < 2*MAGAS//3:
            jkezd=jkezd+subkepernyoSZELES//8
        else:
            jkezd=jkezd+subkepernyoSZELES//4
        subkepernyoSZELES=subkepernyoSZELES//2
    
    render_params_ck = {
        "SZELES": SZELES,
        "MAGAS": MAGAS,
        "kepernyoSZELES": kepernyoSZELES,
        "kepernyoMAGAS": kepernyoMAGAS,
        "ikezd": ikezd,
        "jkezd": jkezd,
        "iveg": iveg,
        "precision": "double" if prec_double else "float",
        "de0" : de0,
        "errormax" : errormax,
        "rs" : rs,
        "a" : a,
        "Q" : Q
    }
    h = render_hash(render_params_ck)
    cached_image = cache_lookup(h)
    return is_black(cached_image)


def filtered_product_with_blackcheck(ranges, repeat, fast_spining):
    """Generate product combinations, skipping branches where check_path_is_black returns True."""
    def build(prefix):
        if len(prefix) == repeat:
            yield prefix
            return
        
        # Check if current prefix leads to black image
        if len(prefix) > 0 and check_path_is_black(prefix, fast_spining):
            return  # Skip this entire branch
        
        for val in ranges:
            candidate = prefix + (val,)
            yield from build(candidate)
    
    yield from build(())


n_depth=16
for i in range(1,n_depth):
    for (fast_spining,) in product([True,False],repeat=1):
        for p in filtered_product_with_blackcheck(range(9), i, fast_spining):
            prec_prev = False
            prec_double = False
            click_count = 0
            kepernyoSZELES = kepernyoSZELES_def
            kepernyoMAGAS  = kepernyoSZELES_def//2
            SZELES = 640
            MAGAS  = 320
            ikezd=0
            jkezd=0
            iveg=kepernyoSZELES_def
            subkepernyoSZELES=kepernyoSZELES_def
            errormax = errormax_def
            de0 = de0_def
            fast = True 
            rs=0.05
            a = 0.0
            Q=0.0
            
            if prec_double == True:
                prec_str = "--double"
            else:
                prec_str = "--float"
            if fast == True:
                de0=de0_def
                errormax=errormax_def
            else:
                de0=0.0001
                errormax=0.000001
            if fast_spining==False:
                a=0.0
            else:
                a=rs/2-0.001
            kepernyoSZELES = kepernyoSZELES_def
            kepernyoMAGAS  = kepernyoSZELES_def//2
            SZELES = 640
            MAGAS  = 320
            iveg=kepernyoSZELES_def
            ikezd=0
            jkezd=0
            subkepernyoSZELES=kepernyoSZELES_def
            for k in p:
                click_x, click_y = tile_number_to_click(k)
                if click_x < SZELES//3:
                    iveg=ikezd+subkepernyoSZELES//2
                elif click_x < 2*SZELES//3:
                    iveg=ikezd+3*subkepernyoSZELES//4
                    ikezd=ikezd+subkepernyoSZELES//4
                else:
                    iveg=ikezd+subkepernyoSZELES
                    ikezd=ikezd+subkepernyoSZELES//2
                if click_y < MAGAS//3:
                    None
                elif click_y < 2*MAGAS//3:
                    jkezd=jkezd+subkepernyoSZELES//8
                else:
                    jkezd=jkezd+subkepernyoSZELES//4
                subkepernyoSZELES=subkepernyoSZELES//2
            
            if prec_double == True:
                prec_str = "--double"
            else:
                prec_str = "--float"
            render_params = {
                    "SZELES": SZELES,
                    "MAGAS": MAGAS,
                    "kepernyoSZELES": kepernyoSZELES,
                    "kepernyoMAGAS": kepernyoMAGAS,
                    "ikezd": ikezd,
                    "jkezd": jkezd,
                    "iveg": iveg,
                    "precision": "double" if prec_double else "float",
                    "de0" : de0,
                    "errormax" : errormax,
                    "rs" : rs,
                    "a" : a,
                    "Q" : Q
                    }
            h = render_hash(render_params)
            cached_image = cache_lookup(h)
            print(p)#,render_params) 
            IMAGE_PATH = f"./web_images/blackhole_cli.png"
            if cached_image:# and Path(cached_image).exists():
                IMAGE_PATH = cached_image
                print("done by cache")
            else:
                subprocess.run(["./main", "--a", str(a),"--rs",str(rs), "--Q", str(Q), "--de0", str(de0), "--errormax", str(errormax),"--SZELES", str(SZELES), "--MAGAS", str(MAGAS), "--kepernyoSZELES", str(kepernyoSZELES), "--kepernyoMAGAS", str(kepernyoMAGAS), "--ikezd", str(ikezd), "--jkezd", str(jkezd), "--iveg", str(iveg), prec_str ])
                subprocess.run(["python", "cli_imagemaker.py"])
                IMAGE_PATH = f"./web_images/blackhole_cli.png"
                cached_path = CACHE_IMG_DIR / f"{h}.png"
                shutil.copy(IMAGE_PATH, cached_path)
                cache_store(h, str(cached_path))
