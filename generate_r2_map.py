import SimpleITK as sitk 
import pandas as pd 
import json
import numpy as np

def main(mask_path:str, r2_path:str, output_path:str, label_dict:dict):
    """
    r2.csv file assumed to have column names "label" and "r2"
    """
    r2_col = "prot_m_vol"
    r2_df = pd.read_csv(r2_path)
    mask_img = sitk.ReadImage(mask_path)

    mask_arr = sitk.GetArrayFromImage(mask_img)
    r2_arr = np.zeros(mask_arr.shape, dtype=np.float32)

    # insert r2 values
    for label, r2 in zip(r2_df["feature"], r2_df[r2_col]):
        # Check if label exists in mapping to avoid KeyError
        if label in label_dict["name_to_id"]:
            label_id = label_dict["name_to_id"][label]
            r2_arr[mask_arr == label_id] = r2
    
    r2_img = sitk.GetImageFromArray(r2_arr)

    # Copy metadata from original
    r2_img.CopyInformation(mask_img)
    
    sitk.WriteImage(r2_img, output_path)

    
if __name__ == "__main__": 
    r2_path = "C:/Users/matti/Documents/Vault/Exjobb/results for heatmap.csv"
    segm_mask_path = "C:/Users/matti/Documents/Vault/Exjobb/Data/2369783_processed_mask.nii.gz"

    output_path = segm_mask_path.replace(".nii.gz", "_r2_map.nii.gz")
    label_dict_path = "VIBE_labels.json"
    with open(label_dict_path) as f:
        label_dict = json.load(f)


    main(segm_mask_path, r2_path, output_path, label_dict)


