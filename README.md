# Origami Analysis Tests

Simulation scripts for origami-inspired structures using **OpenSeesPy**, with a custom OpenSees core provided in the `ZibraMax/OpenSees` repository.

---

## Overview

This project provides:

-   **Bar-and-hinge** and **shell-and-hinge** origami models.
-   Example patterns: **Miura-ori** and a **simple fold**.

**custom-built OpenSees** implementation is required, hosted at [ZibraMax/OpenSees](https://github.com/ZibraMax/OpenSees).

---

## Requirements

-   **Python 3.11**
-   **CMake**, compilers (gcc/g++ or MSVC), and build tools
-   **Dependencies**: `numpy`, `matplotlib`

Install Python packages:

```bash
pip install numpy matplotlib
```

---

## Custom OpenSees Build

### 1. Clone the custom OpenSees repository

```bash
git clone https://github.com/ZibraMax/OpenSees.git
cd OpenSees
```

### 2. Build the OpenSees core and Python module

Follow the build instructions appropriate for your platform (Linux/macOS/Windows), as outlined in the [official OpenSees documentation](https://opensees.github.io/OpenSeesDocumentation/developer/build.html). Typically:

```bash
mkdir build
cd build
cmake ..
cmake --build . --target OpenSees -j4
cmake --build . --target OpenSeesPy -j4
```

### Using WSL:

I used WSL for building opensees. The file `installation_script.sh` of the `ZibraMax/OpenSees` repository includes the step by step commands for compiling the OpenSees code.

### 3. Use the custom module in Python

-   **Option A** — Use `PYTHONPATH` to point Python at your custom build:

    ```bash
    export PYTHONPATH=/path/to/OpenSees/build:$PYTHONPATH
    ```

-   **Option B** — (not recomended) Replace existing `openseespy` installation:

    ```bash
    pip install openseespy
    # Identify its location
    python3 -c "import opensees, inspect; print(inspect.getfile(opensees))"
    # Replace `opensees.so` with your custom-built one
    ```

Note: Ensure the Python interpreter used to build OpenSees matches the one used when running analysis, otherwise segmentation faults may occur.

---

## Running Origami Simulations

Once the custom OpenSees build is set up, you can access to the OriHinge element to create the origami simulations. There are a couple of examples of both <i>bar and hinge</i> and <i>shell and hinge models</i>/

### Miura-ori (bar‑and‑hinge)

```bash
cd bar_and_hinge/miura
python miura.py
```

### Simple Fold (shell‑and‑hinge)

```bash
cd shell_and_hinge/simple_fold
python fold.py
```

## License

MIT License — see the [LICENSE](LICENSE) file for full details.
