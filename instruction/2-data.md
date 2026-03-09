# Part 2 - Data

## Datasets
This part uses 2025 `ParkingDoubleMuonLowMass` prompt-reco MINIAOD samples:

| Era | Dataset | Recorded by CMS [$\mathrm{fb}^{-1}$] | Run range |
| --- | --- | --- | --- |
| 2025F | `/ParkingDoubleMuonLowMass*/Run2025F-PromptReco-v1/MINIAOD` | 27.3 | 396598-397853 |
| 2025G | `/ParkingDoubleMuonLowMass*/Run2025G-PromptReco-v1/MINIAOD` | 23.0 | 397854-398903 |

## CMSSW
```bash
cd /path/to/CMSDAS/data

cmsrel CMSSW_15_0_18
cd CMSSW_15_0_18/src
cmsenv

git clone https://github.com/yiyangzha/Onia2MuMu.git
scram b -j 16

cd Onia2MuMu/Analyzers/MuMu/test
```

- `cmsrel/cmsenv` provides the runtime expected by `cmsRun`.
- `scram b` builds local analyzer plugins and configuration dependencies.

### Local Test
First query available files:
```bash
dasgoclient --query="file dataset=/ParkingDoubleMuonLowMass0/Run2025E-PromptReco-v1/MINIAOD" | head -n 5
```

> #### **Task**
> Replace the placeholder with a real file path from your DAS query and run:
> ```bash
> cmsRun run_upsilon.py inputFiles=<one_MINIAOD_file_from_your_DAS_query>
> ```

After production, inspect the output:
```bash
root -l rootuple.root
```
```cpp
rootuple->cd();

mm_tree->GetEntries();
mm_tree->Show(0)
mm_tree->Scan("dimuon_p4.Pt():dimuon_p4.Rapidity():trigger","","",5);

.q
```

Code-output relation:
- `run_upsilon.py` builds and filters dimuon candidates.
- `MMrootupler` stores selected event content into `rootuple.root/mm_tree`.
- The scanned branches (`dimuon_p4`, `trigger`) are analysis-driving observables used again in later parts.

> #### **Checkpoint**
> Before CRAB submission, confirm:
> - `mm_tree->GetEntries()` is non-zero.
> - `Scan` can read key branches such as `dimuon_p4`, `trigger`, and `vProb`.

### CRAB
Edit `mydata`, `myname`, and the lumi mask in `crab_upsilon.py`, then submit:
```bash
crab submit -c crab_upsilon.py
crab status -d CernJobs/crab_<your_request_name>
```

These lines control dataset input, certified-lumi filtering, and job splitting:
```python
config.Data.inputDataset = mydata
config.Data.lumiMask = 'Cert_Collisions2025_391658_398903_Muon.json'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
```

- CRAB is the scalable path for large data processing.
- `Cert_Collisions2025_391658_398903_Muon.json` is the certified good-luminosity list for this run range.
- It removes luminosity sections with known detector/DAQ/data-quality problems (for example unstable detector conditions, data-taking interruptions, or subdetector quality flags failing certification).
- This ensures that accepted events come from periods where detector performance is validated for physics analysis, so your yield and cross-section normalization remain consistent and reproducible.

## Data Distributions
Run the plotting program:
```bash
cd /path/to/CMSDAS/data
root -l -b -q plot.C
```
```

> #### **Task**
> Add analysis-level selections directly inside the event loop:
> - acceptance kinematics for reconstructed muons: $|\eta|<2.0$ and $p_T>3.1$ GeV,
> - trigger passed,
> - $vProb>0.01$.
>
> Loop location hint:
> ```cpp
> for (Long64_t i = 0; i < n_entries; ++i) {
>   chain.GetEntry(i);
> ...
> }
> ```

> #### **Question**
> 1. Compare distributions before and after your selections: what observables have changed, and why?
> 2. Is the behavior consistent with trigger and vertex-quality expectations?
