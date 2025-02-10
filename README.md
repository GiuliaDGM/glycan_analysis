# glycan_analysis

# 1) **Clone the git repository**

```sh
git clone git@github.com:GiuliaDGM/glycan_analysis.git
cd glycan_analysis
```

# 2) Environment setup with UV

For this project, we will be using **UV** to create and manage our environments.

### Why UV and not Conda?
1. **Speed** – UV is significantly faster than Conda for package installation and dependency resolution.
2. **Lightweight** – Unlike Conda, UV does not require a separate package manager and works efficiently with `pyproject.toml` (file that automatically manages all libraries and dependencies needed).
3. **Better Dependency Management** – UV offers improved dependency resolution and reproducibility compared to Conda.

For more details, visit the **[UV Documentation](https://astral.sh/uv/)**.

---

## How to install UV

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

# 3) **Sync dependencies**

Once inside the repository, sync the environment with:
```sh
uv sync
```

This will install all necessary dependencies from `pyproject.toml`.



# 4) **Running the code**

After setting up your environment with `uv`, you can run the project using:

```sh
uv run src/main.py
```

By default, the output data file is set in **`src/main.py`** at **line 27**:

```python
output_file = "results/IL3_GLY_R1/glycan_ranges.txt"
```

If you want to change the output location, update this line before running the script.
