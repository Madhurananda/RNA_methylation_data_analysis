NanoCNN – A CNN-Based Tool to Detect RNA m⁶A Methylations Using Oxford Nanopore Sequencing

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.6.0-red.svg)](https://pytorch.org/)
[![TensorFlow](https://img.shields.io/badge/TensorFlow-2.x-orange.svg)](https://tensorflow.org/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

**NanoCNN** is an end-to-end deep learning pipeline for detecting RNA modifications (e.g., m⁶A, m⁵C, pseudouridine, inosine) from Oxford Nanopore direct RNA sequencing data. It is designed to operate directly on raw electrical signal data and produce nucleotide-resolution modification predictions.

---

## Overview

NanoCNN provides a complete workflow from raw signal processing to model training and evaluation. The pipeline is built to handle large-scale Nanopore datasets efficiently while maintaining flexibility for different modification types and experimental setups.

Key capabilities:

* Extraction of per-event signal features, normalisation parameters, and alignment data from FAST5 files processed with Tombo
* Signal discretisation into fixed-length histogram representations (default: 100 bins)
* Efficient dataset generation using TFRecord format
* Self-supervised pre-training via masked signal prediction
* Supervised fine-tuning using a bidirectional LSTM architecture
* Modification-specific classification at user-defined nucleotide targets (A, C, G, U)
* Generation of per-site predictions and ROC-based performance evaluation

---

## Installation

### Clone the repository

```bash
git clone https://github.com/Madhurananda/RNA_methylation_data_analysis.git
cd RNA_methylation_data_analysis
```

### External tools required

Ensure the following tools are installed and accessible in your environment:

* **Tombo** – for signal alignment and normalisation
* **minimap2** – for generating BAM alignment files

---

## Dependencies

Core Python libraries:

* **PyTorch** – model architecture and training
* **TensorFlow** – TFRecord handling
* **h5py** – FAST5 file access
* **pysam** – BAM file parsing
* **NumPy / pandas** – data processing
* **scikit-learn** – evaluation metrics and ROC analysis
* **matplotlib** – visualization
* **tqdm, natsort** – utility functions

---

## Pipeline Steps

### 1. FAST5 → TFRecord Generation

Convert FAST5 files and Tombo outputs into TFRecord shards.

Basic usage:

```bash
python generate_tfrecord_bp.py /path/to/fast5_folder /output/tfrecords_dir 8 1 /path/to/alignment.bam 200 20 RNA
```

Parallel processing across directories:

```bash
python do_generate_tfrec_BP_FT_MP.py /basecall_folder ALL_IN_DIR /analysis_folder 8 1 200 20 BP
```

---

### 2. Pre-training (Self-Supervised)

Train a bidirectional LSTM using masked signal reconstruction:

```bash
python TS_preTrainLSTM_bp.py MSE 3 512 0 1500 512
```

Output:

```
Bilst_bp.layers3.hs512.lr1500.b512.GPU.ep*.pt
```

---

### 3. Fine-tuning (Supervised)

Train the classifier for a specific nucleotide modification:

```bash
python TS_finetune_run_MP.py /data/p5/ /ref/IVT_seq.fa model_output 5 3 512 0 F 1500 1024 1 A A1
```

**METHYL_TYPE options:** `A`, `C`, `G`, `U`

---

### 4. Prediction & Evaluation

Run predictions on new data:

```bash
python TS_finetune_pred_run_MP.py /basefolder /pred_output_dir 3 512 0 512 saved_model.pt /ref_genome /test_data_folder A
```

Generate ROC curves:

```bash
python plot_ROC_allTools.py
```

---

## Configuration

Key parameters are defined in `TS_global.py`:

| Variable              | Description                 | Default      |
| --------------------- | --------------------------- | ------------ |
| `dis_size`            | Number of histogram bins    | 100          |
| `input_len`           | Sequence window length      | 31           |
| `m_md_norm_limits`    | Signal normalization range  | [-5, 5, 0.1] |
| `tf_record_max_event` | Maximum events per TFRecord | 500000       |
| `g_dropout_rate`      | Dropout rate                | 0.2          |

Adjust these parameters according to your dataset characteristics before running the pipeline.

---

## Outputs

The pipeline produces several outputs across different stages:

* **TFRecord files** – compressed `.gz` shards containing processed features
* **Pre-trained models** – `.pt` checkpoint files
* **Fine-tuned models** – trained classification models
* **Prediction results** – `.predres` files with per-site predictions
* **ROC curves** – `.png` performance comparison plots

---

## Citation

If you use NanoCNN in your research, please cite:

Pahar, M. and Liu, Q., 2024, November. *NanoCNN: a CNN-based tool to detect RNA m⁶A methylations using Oxford Nanopore sequencing*. In **2024 4th International Conference on Electrical, Computer, Communications and Mechatronics Engineering (ICECCME)** (pp. 1–6). IEEE.

```bibtex id="xv7k3p"
@inproceedings{pahar2024nanocnn,
  title={Nano{CNN}: a {CNN}-based tool to detect {RNA} m$^{\text{6}}$A methylations using {O}xford {N}anopore sequencing},
  author={Pahar, Madhurananda and Liu, Qian},
  booktitle={2024 4th International Conference on Electrical, Computer, Communications and Mechatronics Engineering (ICECCME)},
  pages={1--6},
  year={2024},
  organization={IEEE}
}
```

---

## Notes

* This pipeline assumes familiarity with Nanopore direct RNA sequencing workflows.
* Performance depends on data quality, alignment accuracy, and preprocessing consistency.
* Pre-training is optional but strongly recommended for improved downstream accuracy.

---
