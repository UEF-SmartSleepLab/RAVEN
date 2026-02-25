# RAVEN
**Rhythmic Analysis of Variations in EEG Neural Activity**

RAVEN is an automated analysis tool for processing overnight EEG recordings with a focus on sleep microstructures. It detects and characterizes:

- **Sleep Spindles**
- **K-complexes**
- **Delta Wave Sequences**
- **Cyclic Alternating Pattern (CAP)**

RAVEN is designed to run on single-channel or multi-channel EEG data and supports both research and clinical workflows.

---

## 🚀 Features

- **No prior sleep staging required** — works independently of hypnograms.
- **Versatile data format input** — Accepts EEG data in EDF or in SLF format.  
- **Single-channel compatible** — CAP detection optimized for two channels.  
- **Automated batch (group) analysis** for multi-subject datasets.  
- **Clear event outputs** saved per subject with timestamps.

---

## 📥 Download & Installation

### ✔ Option 1: Run directly from MATLAB  
1. Download all repository files into the **same folder**.  
2. Open MATLAB.  
3. Run the application by opening: RAVEN_desk.mlapp

This is the **fastest** way to use RAVEN because it runs directly in MATLAB without external runtimes.

---

### ✔ Option 2: Use the Installer  
1. Download the repository and the installer package.  
2. Install RAVEN using the provided installer.  

⚠ **Note:**  
The installed version loads MATLAB Runtime from the server, which makes startup **slightly slower** than running the `.mlapp` directly in MATLAB.

---

## 📁 Supported EEG File Formats

RAVEN processes EEG data directly from:

- **EDF** (European Data Format)  
- **SLF** (SleepLab Format)

When using **Group Analysis**, ensure all EDF/SLF files are placed in the **same folder** before starting the analysis.

---

## 📦 Output Files

RAVEN creates:

- `.mat` files for each event type  
- One folder per subject  
- Automatic timestamps to prevent overwriting  

At the end, RAVEN provides a summary including:
- Processed subjects  
- Skipped subjects  
- Any errors encountered


## 👁️ ODIN: Visualization & Post‑Processing

For visual review and post‑processing of RAVEN results, users can optionally use **ODIN**, a companion application available from the UEF-SmartSleepLab GitHub main page.

ODIN allows you to:
- Visualize detected events (spindles, K-complexes, Delta waves, CAP) on the EEG signal  
- Inspect, validate, and manually edit events  
- Perform additional post‑processing and quality control  

### 🤝 RAVEN + ODIN Workflow
1. Run RAVEN to generate detection result files (.mat).
2. Open ODIN and load these files.
3. Review, visualize, and optionally refine the detections.
