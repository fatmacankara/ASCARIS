import ssl
import requests
import time
import json

def get_pdb_ids(protein_id):
    POLLING_INTERVAL = 3
    API_URL = "https://rest.uniprot.org"

    def submit_id_mapping(fromDB, toDB, ids):
        r = requests.post(
            f"{API_URL}/idmapping/run", data={"from": fromDB, "to": toDB, "ids": ids},
        )
        r.raise_for_status()
        return r.json()["jobId"]

    def get_id_mapping_results(job_id):
        while True:
            r = requests.get(f"{API_URL}/idmapping/status/{job_id}")
            r.raise_for_status()
            job = r.json()
            if "jobStatus" in job:
                if job["jobStatus"] == "RUNNING":
                    print(f"Retrying in {POLLING_INTERVAL}s")
                    time.sleep(POLLING_INTERVAL)
                else:
                    raise Exception(job["jobStatus"])
            else:
                return job

    job_id = submit_id_mapping(
        fromDB="UniProtKB_AC-ID", toDB="PDB", ids=protein_id
    )
    results = get_id_mapping_results(job_id)
    pdbs = {}
    for i in results['results']:
        if i['from'] not in pdbs.keys():
            pdbs[i['from']] = []
            pdbs[i['from']].append(i['to'])
        else:
            pdbs[i['from']].append(i['to'])
    return pdbs
