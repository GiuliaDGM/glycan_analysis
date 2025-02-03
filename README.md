# glycan_analysis

# Environment Setup with UV

For this project, we will be using **UV** to create and manage our environments.

### Why UV and not Conda?
1. **Speed** – UV is significantly faster than Conda for package installation and dependency resolution.
2. **Lightweight** – Unlike Conda, UV does not require a separate package manager and works efficiently with `pyproject.toml`.
3. **Better Dependency Management** – UV offers improved dependency resolution and reproducibility compared to Conda.

For more details, visit the **[UV Documentation](https://astral.sh/uv/)**.

---

## How to Install UV

### **On macOS and Linux:**
```sh
curl -LSsf https://astral.sh/uv/install.sh | sh
```

### **On Windows:**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### **Or install from PyPI:**
```sh
# Using pip
pip install uv

# Using pipx
pipx install uv
```

After installation, you can update UV with:
```sh
uv self update
```

For more installation details, check the [UV GitHub page](https://github.com/astral-sh/uv).

---

## Clone the Git Repository

```sh
git clone git@github.com:GiuliaDGM/glycan_analysis.git
cd glycan_analysis
```

---

## Sync Dependencies

Once inside the repository, sync the environment with:
```sh
uv sync
```

This will install all necessary dependencies from `pyproject.toml`.

---

## Running the Code

After setting up your environment with `uv`, you can run the project using:

```sh
uv run src/main.py
```

### **Difference Between `python src/main.py` and `uv run src/main.py`**

1. **`python src/main.py`**
   - Runs the script using the globally installed Python interpreter or the one in your virtual environment.
   - Requires you to manually activate the virtual environment before running (`source .venv/bin/activate` or `.\.venv\Scripts\activate` on Windows).

2. **`uv run src/main.py`**
   - Runs the script using the environment managed by `uv`.
   - Ensures that dependencies from `pyproject.toml` are correctly resolved and installed before running the script.
   - Does **not** require activating a virtual environment manually.

### **Configuring the Output File**

By default, the output data file is set in **`src/main.py`** at **line 27**:

```python
output_file = "results/IL3_GLY_R1/glycan_ranges.txt"
```

If you want to change the output location, update this line before running the script.
