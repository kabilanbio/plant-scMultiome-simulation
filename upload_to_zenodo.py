"""
upload_to_zenodo.py
Automated script to upload prepared deposition archives to Zenodo Record 23068394.
Usage:
    python upload_to_zenodo.py --token YOUR_ZENODO_PERSONAL_ACCESS_TOKEN
"""

import os
import sys
import argparse
import requests

RECORD_ID = "23068394"
BASE_URL = f"https://zenodo.org/api/deposit/depositions/{RECORD_ID}"

def parse_args():
    parser = argparse.ArgumentParser(description="Upload dataset bundle to Zenodo record.")
    parser.add_argument("--token", type=str, default=os.environ.get("ZENODO_TOKEN"),
                        help="Zenodo Personal Access Token (get one from https://zenodo.org/account/settings/applications/tokens/new/)")
    parser.add_argument("--bundle_dir", type=str, default=r"G:\PhD\data_simulation_comparison\scMultiSim\zenodo_deposition_bundle",
                        help="Directory containing the prepared zip archives and manifest")
    return parser.parse_args()

def main():
    args = parse_args()
    if not args.token:
        print("[ERROR] No Zenodo token provided.")
        print("Please provide via --token or set the ZENODO_TOKEN environment variable.")
        print("You can generate an upload token at: https://zenodo.org/account/settings/applications/tokens/new/")
        print("\nAlternatively, you can manually upload the prepared files via the Zenodo Web UI:")
        print(f"1. Navigate to: https://zenodo.org/records/{RECORD_ID}")
        print("2. Click 'Edit' (or 'New Version')")
        print(f"3. Upload the 3 zip archives and manifest from: {args.bundle_dir}")
        print("4. Save and Publish.")
        sys.exit(1)

    headers = {"Authorization": f"Bearer {args.token}"}
    
    # Check deposition status
    res = requests.get(BASE_URL, headers=headers)
    if res.status_code == 404:
        print(f"[INFO] Accessing record {RECORD_ID} as a published record...")
        pub_url = f"https://zenodo.org/api/records/{RECORD_ID}"
        res_pub = requests.get(pub_url, headers=headers)
        if res_pub.status_code == 200:
            print("[INFO] Record is published. Creating a new draft version...")
            res_new = requests.post(f"{pub_url}/actions/newversion", headers=headers)
            if res_new.status_code == 201:
                new_draft_url = res_new.json()["links"]["latest_draft"]
                print(f"[SUCCESS] New version draft created: {new_draft_url}")
                deposition_id = new_draft_url.split("/")[-1]
                upload_files_to_deposition(deposition_id, args.bundle_dir, headers)
                return
            else:
                print(f"[ERROR] Failed to create new version draft: {res_new.status_code} {res_new.text}")
                sys.exit(1)
    elif res.status_code == 200:
        upload_files_to_deposition(RECORD_ID, args.bundle_dir, headers)
    else:
        print(f"[ERROR] Could not query Zenodo: {res.status_code} {res.text}")
        sys.exit(1)

def upload_files_to_deposition(dep_id, bundle_dir, headers):
    bucket_url = f"https://zenodo.org/api/deposit/depositions/{dep_id}"
    dep_res = requests.get(bucket_url, headers=headers).json()
    bucket = dep_res.get("links", {}).get("bucket")
    
    files_to_upload = [
        "paga_lineage_trees_and_metadata.zip",
        "paga_processed_anndata_h5ad.zip",
        "paga_publication_figures_600dpi.zip",
        "ZENODO_DATASET_MANIFEST.md"
    ]
    
    print(f"\n[INFO] Uploading files to Zenodo deposition {dep_id}...")
    for fn in files_to_upload:
        fp = os.path.join(bundle_dir, fn)
        if not os.path.exists(fp):
            print(f"[WARNING] File not found: {fp}, skipping...")
            continue
        print(f"Uploading {fn} ({os.path.getsize(fp) / (1024*1024):.2f} MB)...")
        if bucket:
            upload_url = f"{bucket}/{fn}"
            with open(fp, "rb") as f:
                r = requests.put(upload_url, data=f, headers=headers)
        else:
            upload_url = f"https://zenodo.org/api/deposit/depositions/{dep_id}/files"
            with open(fp, "rb") as f:
                r = requests.post(upload_url, data={"name": fn}, files={"file": f}, headers=headers)
                
        if r.status_code in [200, 201]:
            print(f"  --> Successfully uploaded {fn}")
        else:
            print(f"  [ERROR] Upload failed for {fn}: {r.status_code} {r.text}")

    print("\n[SUCCESS] All files uploaded to Zenodo! Log into https://zenodo.org to review and publish the draft.")

if __name__ == "__main__":
    main()
