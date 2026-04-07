# **OpenIBIS: Open Implementation of BIS Algorithm**
![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)
![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-orange?logo=mathworks&logoColor=white)

---

## Overview
 
The BIS (Bispectral Index) is one of the most widely used clinical indicators for monitoring depth of anesthesia during surgery. Despite its widespread adoption, the underlying algorithm has historically been proprietary — limiting reproducibility, clinical understanding, and research into closed-loop anesthesia control.

**OpenIBIS** provides a transparent, documented Python and MATLAB implementation of the BIS algorithm, enabling:
 
- Researchers to reproduce and validate BIS scores from arbitrary EEG recordings
- Clinicians and engineers to understand the physiological basis of the index
- Development of improved anesthesia monitoring and closed-loop control systems
- Reduction of postoperative cognitive dysfunction (POCD) risk through better monitoring

This implementation is based on the forensic disassembly of an A-2000 BIS monitor described by Connor (2022), which allowed BIS scores to be calculated from arbitrary EEG data for the first time.

## Key features
 
- **BIS score computation** from raw single-channel EEG time-series
- **Python implementation** (`openibis.py`) — actively developed
- **MATLAB implementation** (`openibis.m`) — fully functional
- Signal preprocessing: bandpass filtering, artifact rejection
- Spectral feature extraction: power spectral density, bispectrum, SynchFastSlow
- BIS index computation: weighted combination of spectral features

 ---
 
## Quickstart
 
### Requirements
 
```
Python >= 3.9
numpy
scipy
matplotlib
```

Install dependencies:
 
```bash
git clone https://github.com/galadriana/openibis.git
cd openibis
pip install -r requirements.txt
```
 
### Basic usage (Python)
 
```python
import numpy as np
from openibis import compute_bis
 
# Load your EEG signal: 1D array, single channel
# fs: sampling frequency in Hz (typically 128 or 256 Hz for BIS monitors)
eeg_signal = np.load("your_eeg_data.npy")
fs = 128  # Hz
 
# Compute BIS score (returns value between 0–100)
bis_score = compute_bis(eeg_signal, fs=fs)
print(f"BIS score: {bis_score:.1f}")
```
 
### Basic usage (MATLAB)
 
```matlab
% Load EEG signal as a column vector
load('your_eeg_data.mat');   % variable: eeg_signal
fs = 128;                    % sampling frequency in Hz
 
% Compute BIS score
bis_score = openibis(eeg_signal, fs);
fprintf('BIS score: %.1f\n', bis_score);
```
 
---
 
## BIS score interpretation
 
| BIS range | Clinical state |
|-----------|---------------|
| 100 | Fully awake |
| 80–100 | Light/no sedation |
| 60–80 | Light–moderate anesthesia |
| 40–60 | General anesthesia (target range) |
| < 40 | Deep anesthesia / burst suppression |
| 0 | Isoelectric EEG (flat line) |
 
---
 
## Algorithm overview
 
The BIS algorithm combines three main spectral features extracted from the EEG:
 
1. **Beta ratio** — log ratio of power in beta band (30–47 Hz) vs alpha/theta bands (11–20 Hz). Captures cortical activation.
2. **SynchFastSlow** — bispectral measure of synchrony between fast (40–47 Hz) and slow (0.5–47 Hz) components.
3. **Burst suppression ratio (BSR)** — fraction of isoelectric (suppressed) epochs in the past 63 seconds.
 
These are combined via a weighted polynomial into a single 0–100 index. See Connor (2022) for the full derivation.
 
---
 
## Repository structure
 
```
openibis/
├── openibis.py          # Python implementation
├── openibis.m           # MATLAB implementation
├── requirements.txt     # Python dependencies
├── tests/
│   └── test_openibis.py # Unit tests (pytest)
├── notebooks/
│   └── demo.ipynb       # Example usage + visualisation
├── .github/
│   └── workflows/
│       └── ci.yml       # GitHub Actions CI
└── README.md
```
 
---
 
## Validation & limitations
 
**What this implementation reproduces:**
- BIS scores within ±3 points of published reference values on standard EEG test signals (see Connor 2022, Table 1).
 
**Known limitations:**
- The Python version is still under active development — some edge cases in artifact rejection differ from the MATLAB reference implementation.
- This implementation does not include the proprietary A-2000 signal quality index (SQI) used in commercial monitors.
- Not validated for clinical use. This is a research tool only.
 
**Not for clinical decision-making.** OpenIBIS is intended for research and educational purposes. It has not been evaluated as a medical device under EU MDR / FDA 510(k) regulations.
 
---
 
## Citation
 
If you use OpenIBIS in your research, please cite both the original algorithm paper and this repository:
 
**Algorithm:**
```bibtex
@article{connor_open_2022,
  title   = {Open Reimplementation of the {BIS} Algorithms for Depth of Anesthesia},
  author  = {Connor, Christopher W.},
  journal = {Anesthesia \& Analgesia},
  volume  = {135},
  number  = {4},
  pages   = {855--864},
  year    = {2022},
  doi     = {10.1213/ANE.0000000000006119}
}
```
 
**This repository:**
```bibtex
@software{vargas_openibis_2024,
  author  = {Vargas, Gabriela},
  title   = {{OpenIBIS}: Open Implementation of the BIS Algorithm},
  url     = {https://github.com/galadriana/openibis},
  year    = {2024}
}
```
 
---
 
## Related work & context
 
This project is part of a broader effort to make anesthesia monitoring more transparent and reproducible. Related topics:
 
- Closed-loop anesthesia control (CLAC)
- EEG-based consciousness monitoring
- Postoperative cognitive dysfunction (POCD) prevention
- Bispectrum analysis in biomedical signal processing
 
---
 
## Contributing
 
Contributions are welcome! Please open an issue first to discuss what you'd like to change.
 
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m 'Add your feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a pull request
 
Please include tests for any new functionality.
 
---
 
## License
 
This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
 
---
 
## Contact
 
**Gabriela Vargas González**
Computational neuroscience · Biomedical data science
📧 gadriana.diez@gmail.com
🔗 [GitHub](https://github.com/galadriana) · [LinkedIn](https://www.linkedin.com/in/gabriela-vargas-gonz%C3%A1lez-1344b2118/)




