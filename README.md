# glycan_analysis

Here is the **README** section based on your request:

---

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
git clone <repository_url>
cd <repository_name>
```

---

## Sync Dependencies

Once inside the repository, sync the environment with:
```sh
uv sync
```

This will install all necessary dependencies from `pyproject.toml`.

---

This ensures a **fast**, **lightweight**, and **reproducible** setup for your project. 🚀