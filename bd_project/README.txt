BD project: Protein Function Prediction

Project Overview
----------------
This project focuses on protein function prediction using multi-label classification techniques. Three different ontologies are considered:
- CC (Cellular Component)
- MF (Molecular Function)
- BP (Biological Process)

Each ontology is modeled separately, with the final submission combining their predictions.

File Structure
--------------
- datasets_creation_CC+MF+BP: Script for dataset preprocessing.
- model_bp.pth, model_cc.pth, model_mf.pth: Best trained models for BP, CC, and MF ontologies.
- mlb_bp.pkl, mlb_cc.pkl, mlb_mf.pkl: Multi-label binarizers for encoding/decoding labels.
- BP_curves, CC_curves, MF_curves: Evaluation results and performance curves.
- plot_bp, plot_cc, plot_mf, plot_zoom: Visualization files for model performance.
- output: Script which uses the best models to produce submission
- submission: Final submission file with predictions.
- BD_report: Project report detailing methodology and results.

Usage
-----
1. Preprocessing:
Run `datasets_creation_CC+MF+BP` to generate the dataset and repeat our analysis and trials.

2. Model Inference:
Load a trained model and perform predictions. Ensure to use the corresponding multi-label binarizer for decoding outputs.
The file output handles the model trained.

3. Submission:
The final predictions are stored in `submission` in bot .txt as well as .tsv extension to better comparison with the baselines.

4. Note that `output` takes on `test_ids.txt` stored in the path `../data/test`.

Results
-------
Model performance is evaluated using precision, recall, and F1, however the highest attention is on (weighted) PR Curves and F1-scores. Detailed results and analysis are available in `BD_report`.

Authors
-------
- Valeria Ferlin
- Mattia Piazza
- Elia Tiso