# Phase 0 setup: WSL2, Docker, Nextflow

Nextflow is a JVM application that assumes a POSIX shell, POSIX file permissions, and symlinks. It has
no native Windows support — not as an oversight, but because its execution model (staging task inputs as
symlinks into per-task working directories) has no clean Windows equivalent. So the pipeline runs inside
WSL2, which is a real Linux kernel, not an emulation layer.

This is also what you want for the interview: every bioinformatics pipeline you will touch professionally
runs on Linux.

Estimated time: 30–45 minutes, including one reboot.

---

## Where files live, and why it matters

| What | Where | Why |
|---|---|---|
| Pipeline source (this repo) | `C:\Users\Chris\Documents\stpipe` → `/mnt/c/Users/Chris/Documents/stpipe` in WSL | So your editor and the agent can edit it directly |
| Nextflow `work/` directory | `~/nxf-work/stpipe` (native Linux ext4) | Nextflow symlinks task inputs; symlinks on the `/mnt/c` bridge are slow and semantically lossy |
| Downloaded data | `~/nxf-data/stpipe` (native Linux ext4) | Multi-GB files; `/mnt/c` IO is roughly an order of magnitude slower |

`setup/bootstrap_wsl.sh` sets `NXF_WORK` and a `STPIPE_DATA` variable for you so you never have to think
about this again. **Be able to explain this tradeoff** — it is a concrete, credible answer to "what
environment issues did you hit?"

---

## Step 1 — Install WSL2 + Ubuntu

Open PowerShell **as Administrator** (right-click Start → Terminal (Admin)) and run:

```powershell
wsl --install -d Ubuntu
```

Then **reboot**. On first launch Ubuntu asks you to create a UNIX username and password — this is local
to the WSL distro and unrelated to your Windows account. Remember the password; you need it for `sudo`.

Verify afterwards, in a normal PowerShell:

```powershell
wsl --list --verbose
```

You want to see `Ubuntu` with `VERSION 2`. If it says `VERSION 1`, run `wsl --set-version Ubuntu 2`.

**If `wsl --install` fails** it is almost always hardware virtualisation being disabled. Check Task Manager
→ Performance → CPU → "Virtualization: Enabled". If disabled, enable Intel VT-x / AMD-V in your BIOS.

---

## Step 2 — Install Docker Desktop

Download from https://www.docker.com/products/docker-desktop/ and install with the
**"Use WSL 2 instead of Hyper-V"** option checked.

Then, critically: **Settings → Resources → WSL Integration → enable `Ubuntu`**, and Apply & Restart.
Without this the `docker` command does not exist inside Ubuntu.

Also give Docker enough memory for Tangram — **Settings → Resources → Memory: at least 12 GB** (16 GB if
you have 32 GB of RAM). The 73k-cell reference is the memory-hungry part of this project.

Verify inside Ubuntu:

```bash
docker run --rm hello-world
```

---

## Step 3 — Bootstrap the Linux toolchain

Open Ubuntu (Start menu → "Ubuntu", or `wsl` from PowerShell) and run:

```bash
cd /mnt/c/Users/Chris/Documents/stpipe
bash setup/bootstrap_wsl.sh
```

This installs, idempotently:

- build essentials + `curl`, `unzip`, `graphviz` (for rendering the pipeline DAG)
- a JDK (Nextflow needs Java 17+)
- Nextflow, into `~/.local/bin`
- `micromamba`, a fast standalone conda-compatible package manager, into `~/.local/bin`
- the `stpipe` analysis environment from `env/environment.yml`
- `NXF_WORK` / `STPIPE_DATA` exports appended to `~/.bashrc`

Then reload your shell and verify:

```bash
source ~/.bashrc
nextflow -version
micromamba --version
micromamba activate stpipe && python -c "import scanpy, squidpy; print(scanpy.__version__, squidpy.__version__)"
```

---

## Step 4 — Phase 0 done

You are finished with Phase 0 when all three of these succeed inside Ubuntu:

```bash
nextflow -version                      # 24.10 or newer
docker run --rm hello-world            # "Hello from Docker!"
micromamba activate stpipe && python -c "import scanpy, squidpy, torch"
```

Tick the Phase 0 boxes in [`../PROJECT_PLAN.md`](../PROJECT_PLAN.md) and move to Phase 1.

---

## Day-to-day workflow from here

```bash
wsl                                          # from PowerShell, or open Ubuntu directly
cd /mnt/c/Users/Chris/Documents/stpipe
micromamba activate stpipe

python bin/fetch_data.py --all               # once, Phase 1
python bin/sc_qc.py --input ... --outdir ...  # standalone script iteration
nextflow run main.nf -profile docker         # the pipeline
nextflow run main.nf -profile test,docker    # the fast smoke test
```

Edit files in Cursor on the Windows side; run them in WSL. Both see the same files.

---

## Troubleshooting

**`docker: command not found` inside Ubuntu** — WSL Integration is not enabled for Ubuntu in Docker
Desktop settings, or Docker Desktop is not running.

**Nextflow errors about file permissions or `Operation not permitted`** — you are running with `work/` on
`/mnt/c`. Confirm `echo $NXF_WORK` points at a `/home/...` path.

**`Cannot allocate memory` / killed Tangram process** — raise the Docker Desktop memory limit, and confirm
with `free -h` inside WSL.

**WSL eats disk and never gives it back** — the ext4 virtual disk grows but does not auto-shrink. Reclaim
with `wsl --shutdown` then `Optimize-VHD` (Windows Pro) or `diskpart compact vdisk`.

**Slow first `import scanpy`** — normal, it is compiling nothing but importing a large dependency tree.
Subsequent imports are cached.
