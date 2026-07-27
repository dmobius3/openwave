# OpenWave System Architecture

## Modular Design & Development Roadmap

This diagram illustrates the architecture of the OpenWave system, broken down into the following system modules:

![ROADMAP](images/modules.png)

## Scalability & Performance

- Support increasing simulation resolution to handle extreme granularity of Planck-scale interactions
- Efficient handling of large particle counts and ultra-small wavelength resolution
- GPU optimized parallel processing for computational performance

## Tech Stack

- **Primary Language**:
  - Python (>=3.12)
- **Parallel Processing**:
  - Taichi Python Acceleration: GPU optimization for computationally intensive wave simulations
- **Math/Physics Libraries**:
  - NumPy, SciPy
- **Visualization**:
  - Taichi: 3D rendering
  - Matplotlib: numerical analysis plots and cross-sectional graphs
  - Export of 3D images and GIFs for visual inspection
- **Data Output**:
  - Numerical datasets, graphs, and analysis reports in open formats (CSV, JSON, PNG, STL)

---

**Deep readers and AI agents**: the full map of OpenWave's key documents, and the order to read them in, is in [`CLAUDE.md`](CLAUDE.md). The AI-collaboration contract is [`AI_HYGIENE.md`](AI_HYGIENE.md).
