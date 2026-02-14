#!/bin/bash
# Download and prepare multiple EOS tables from stellarcollapse.org
# for robustness testing (Week 3-4)

set -e  # Exit on error

# Base URL for EOS tables
BASE_URL="https://stellarcollapse.org"

# Output directory
SCRIPT_DIR=`dirname -- $0`
EOS_DIR="${SCRIPT_DIR}/../data/eos_tables"
mkdir -p "$EOS_DIR"

echo "======================================================================="
echo "DOWNLOADING EOS TABLES FROM STELLARCOLLAPSE.ORG"
echo "======================================================================="
echo ""
echo "This script downloads EOS tables for robustness testing."
echo "Files are ~200-350 MB compressed, ~2-4 GB uncompressed."
echo ""

# Function to download and extract EOS table
download_eos() {
    local name=$1
    local file=$2
    local url=$3
    
    echo "-----------------------------------------------------------------------"
    echo "Processing: $name"
    echo "-----------------------------------------------------------------------"
    
    # Check if already exists
    if [ -f "$EOS_DIR/${file%.bz2}" ]; then
        echo "✓ Already exists: $EOS_DIR/${file%.bz2}"
        echo "  Skipping download."
        return
    fi
    
    # Download
    echo "→ Downloading: $file"
    if [ -f "$EOS_DIR/$file" ]; then
        echo "  (compressed file exists, skipping download)"
    else
        wget -q --show-progress -O "$EOS_DIR/$file" "$url" || {
            echo "✗ Download failed: $url"
            return 1
        }
    fi
    
    # Extract
    echo "→ Extracting..."
    bunzip2 -k "$EOS_DIR/$file"
    
    # Get file size
    size=$(du -h "$EOS_DIR/${file%.bz2}" | cut -f1)
    echo "✓ Complete: ${file%.bz2} ($size)"
    
    # Optional: Remove compressed file to save space
    # rm "$EOS_DIR/$file"
}

# ===== TIER 1: MUST TEST =====
echo "===== TIER 1: CRITICAL FOR ROBUSTNESS ====="
echo ""

# LS180 - SOFT
download_eos \
    "LS180 (K=180 MeV, SOFT)" \
    "LS180_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2" \
    "$BASE_URL/EOS/LS180_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2"

# LS220 - MEDIUM (already have, but document for completeness)
if [ -f "$EOS_DIR/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5" ]; then
    echo "-----------------------------------------------------------------------"
    echo "LS220 (K=220 MeV, MEDIUM)"
    echo "-----------------------------------------------------------------------"
    echo "✓ Already validated in Week 1-2"
    echo "  M_max = 2.44 M☉ (measured)"
    echo ""
fi

# SFHo - STIFF
download_eos \
    "SFHo (Steiner-Fischer-Hempel-Oertel, STIFF)" \
    "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5.bz2" \
    "$BASE_URL/~evanoc/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5.bz2"

# DD2 - STIFF
download_eos \
    "DD2 (Density-Dependent RMF, STIFF)" \
    "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2" \
    "$BASE_URL/~evanoc/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5.bz2"

# ===== TIER 2: NICE TO HAVE =====
echo ""
echo "===== TIER 2: EXTENDED TESTING (OPTIONAL) ====="
echo ""

# LS375 - VERY STIFF
download_eos \
    "LS375 (K=375 MeV, VERY STIFF)" \
    "LS375_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2" \
    "$BASE_URL/EOS/LS375_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2"

# HShen - DIFFERENT FRAMEWORK
download_eos \
    "HShen (H. Shen et al. 2011)" \
    "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5.bz2" \
    "$BASE_URL/~evanoc/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5.bz2"

echo ""
echo "======================================================================="
echo "DOWNLOAD COMPLETE"
echo "======================================================================="
echo ""
echo "Next steps:"
echo "  1. Run: python scripts/extract_all_eos.py"
echo "     (extracts cold beta-equilibrium slices from all tables)"
echo ""
echo "  2. Run: python tests/test_multi_eos_robustness.py"
echo "     (compares TOV results across all EOS)"
echo ""
echo "Location: $EOS_DIR/"
ls -lh "$EOS_DIR"/*.h5 2>/dev/null || echo "  (no .h5 files yet)"
echo ""
