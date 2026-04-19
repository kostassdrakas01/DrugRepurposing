# GTEx V8 Tissue IDs Reference

Use the following IDs with the `--tissue` flag. You can specify multiple tissues separated by spaces.

## 🧠 Brain & Central Nervous System
- `Brain_Amygdala`
- `Brain_Anterior_cingulate_cortex_BA24`
- `Brain_Caudate_basal_ganglia`
- `Brain_Cerebellar_Hemisphere`
- `Brain_Cerebellum`
- `Brain_Cortex`
- `Brain_Frontal_Cortex_BA9`
- `Brain_Hippocampus`
- `Brain_Hypothalamus`
- `Brain_Nucleus_accumbens_basal_ganglia`
- `Brain_Putamen_basal_ganglia`
- `Brain_Spinal_cord_cervical_c-1`
- `Brain_Substantia_nigra`
- `Nerve_Tibial`

## 🫀 Cardiovascular & Respiratory
- `Artery_Aorta`
- `Artery_Coronary`
- `Artery_Tibial`
- `Heart_Atrial_Appendage`
- `Heart_Left_Ventricle`
- `Lung`

## 🍕 Digestive & Metabolic
- `Adipose_Subcutaneous`
- `Adipose_Visceral_Omental`
- `Colon_Sigmoid`
- `Colon_Transverse`
- `Esophagus_Gastroesophageal_Junction`
- `Esophagus_Mucosa`
- `Esophagus_Muscularis`
- `Liver`
- `Minor_Salivary_Gland`
- `Pancreas`
- `Small_Intestine_Terminal_Ileum`
- `Stomach`

## 🚻 Reproductive & Urinary
- `Adrenal_Gland`
- `Bladder`
- `Cervix_Ectocervix`
- `Cervix_Endocervix`
- `Fallopian_Tube`
- `Kidney_Cortex`
- `Ovary`
- `Prostate`
- `Testis`
- `Uterus`
- `Vagina`

## 🩸 Blood & Immune
- `Spleen`
- `Whole_Blood`
- `Cells_EBV-transformed_lymphocytes`

## 🦵 Musculoskeletal & Skin
- `Muscle_Skeletal`
- `Skin_Not_Sun_Exposed_Suprapubic`
- `Skin_Sun_Exposed_Lower_leg`
- `Cells_Cultured_fibroblasts`

## 🧬 Endocrine & Other
- `Pituitary`
- `Thyroid`
- `Breast_Mammary_Tissue`

---

**Tip**: To run an analysis for multiple tissues (e.g., Liver and Lung), use:
`python main.py "Ibuprofen" --tissue Liver Lung
