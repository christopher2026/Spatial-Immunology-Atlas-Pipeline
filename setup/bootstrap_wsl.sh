#!/usr/bin/env bash
# Bootstrap the Linux toolchain for the spatial immune atlas pipeline.
# Idempotent: safe to re-run. See setup/WSL_SETUP.md for context.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_BIN="$HOME/.local/bin"
NXF_WORK_DIR="$HOME/nxf-work/stpipe"
STPIPE_DATA_DIR="$HOME/nxf-data/stpipe"
ENV_NAME="stpipe"

log() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }

mkdir -p "$LOCAL_BIN" "$NXF_WORK_DIR" "$STPIPE_DATA_DIR"
export PATH="$LOCAL_BIN:$PATH"

# --- system packages ---------------------------------------------------------
log "Installing system packages (needs sudo)"
sudo apt-get update -qq
sudo apt-get install -y -qq \
  build-essential curl wget unzip git ca-certificates \
  graphviz default-jre-headless bzip2

# --- nextflow ----------------------------------------------------------------
if command -v nextflow >/dev/null 2>&1; then
  log "Nextflow already installed: $(nextflow -version 2>&1 | grep -oP 'version \K[0-9.]+' || echo present)"
else
  log "Installing Nextflow into $LOCAL_BIN"
  curl -fsSL https://get.nextflow.io | bash
  mv nextflow "$LOCAL_BIN/nextflow"
  chmod +x "$LOCAL_BIN/nextflow"
fi

# --- micromamba --------------------------------------------------------------
# micromamba is a single static binary that speaks the conda ecosystem. Preferred
# over a full Miniconda install: no base environment, no shell hijacking, ~10x faster solves.
if command -v micromamba >/dev/null 2>&1; then
  log "micromamba already installed: $(micromamba --version)"
else
  log "Installing micromamba into $LOCAL_BIN"
  curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$HOME/.local" bin/micromamba
  chmod +x "$LOCAL_BIN/micromamba"
fi

export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$LOCAL_BIN/micromamba" shell hook --shell bash)"

# --- analysis environment ----------------------------------------------------
if micromamba env list | grep -qE "^\s+${ENV_NAME}\s"; then
  log "Updating existing '${ENV_NAME}' environment"
  micromamba env update -y -n "$ENV_NAME" -f "$REPO_DIR/env/environment.yml"
else
  log "Creating '${ENV_NAME}' environment (this takes a few minutes)"
  micromamba create -y -n "$ENV_NAME" -f "$REPO_DIR/env/environment.yml"
fi

# --- shell configuration -----------------------------------------------------
BASHRC="$HOME/.bashrc"
MARKER="# >>> stpipe bootstrap >>>"
if grep -qF "$MARKER" "$BASHRC" 2>/dev/null; then
  log "~/.bashrc already configured"
else
  log "Appending environment exports to ~/.bashrc"
  cat >> "$BASHRC" <<EOF

$MARKER
export PATH="\$HOME/.local/bin:\$PATH"
export MAMBA_ROOT_PREFIX="\$HOME/micromamba"
eval "\$(micromamba shell hook --shell bash)"
# Keep Nextflow's work dir and downloaded data on native ext4, not the /mnt/c bridge.
export NXF_WORK="$NXF_WORK_DIR"
export STPIPE_DATA="$STPIPE_DATA_DIR"
alias stpipe='cd $REPO_DIR && micromamba activate $ENV_NAME'
# <<< stpipe bootstrap <<<
EOF
fi

# --- verification ------------------------------------------------------------
log "Verification"
printf 'java      : %s\n' "$(java -version 2>&1 | head -n1)"
printf 'nextflow  : %s\n' "$(nextflow -version 2>&1 | grep -i version | head -n1 | tr -s ' ')"
printf 'micromamba: %s\n' "$(micromamba --version)"
printf 'NXF_WORK  : %s\n' "$NXF_WORK_DIR"
printf 'STPIPE_DATA: %s\n' "$STPIPE_DATA_DIR"

if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
  printf 'docker    : %s\n' "$(docker --version)"
else
  warn "Docker is not available inside WSL."
  warn "Install Docker Desktop, then enable Settings -> Resources -> WSL Integration -> Ubuntu."
fi

log "Done. Run 'source ~/.bashrc' then 'stpipe' to start working."
