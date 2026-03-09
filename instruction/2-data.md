# Part 2 - Data

## Datasets
We will use 2025 `ParkingDoubleMuonLowMass` samples in MINIAOD format:

| Era | Dataset | Integrated Luminosity Recorded by CMS [$\mathrm{fb}^{-1}$] | Run range |
| --- | --- | --- | --- |
| 2025F | `/ParkingDoubleMuonLowMass*/Run2025F-PromptReco-v1/MINIAOD` | 27.3 | 396598-397853 |
| 2025G | `/ParkingDoubleMuonLowMass*/Run2025G-PromptReco-v1/MINIAOD` | 23.0 | 397854-398903 |

## CMSSW
```bash
cd /path/to/CMSDAS_Upsilon/data

cmsrel CMSSW_15_0_18
cd CMSSW_15_0_18/src
cmsenv

git clone https://github.com/yiyangzha/Onia2MuMu.git Analyzers/MuMu
scram b -j 16

cd Analyzers/MuMu/test
```

- `cmsrel/cmsenv` provides the runtime expected by `cmsRun`.
- `scram b` builds local analyzer plugins and configuration dependencies.

### Local Test
First query available files:
```bash
dasgoclient --query="file dataset=/ParkingDoubleMuonLowMass0/Run2025E-PromptReco-v1/MINIAOD" | head -n 5
```

> #### **Task**
> Replace the placeholder with a real file path from above DAS query output and run:
> ```bash
> cmsRun run_upsilon.py inputFiles=<one_MINIAOD_file_from_your_DAS_query>
> ```

After production, inspect the output:
```bash
root -l rootuple.root
...
root [1] rootuple->cd()
root [2] mm_tree->GetEntries()
root [3] mm_tree->Show(0)
root [4] mm_tree->Scan("dimuon_p4.Pt():dimuon_p4.Rapidity():trigger","","",5)
root [5] .q
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
In `crab_upsilon.py`, the following lines control request name (`myname`), dataset input (`mydata`), certified-lumi filtering (`lumiMask`), and output directory (`outLFNDirBase`):
```python
myname='CMSDAS_Upsilon_2025E_ParkingDoubleMuonLowMass0_v1'
mydata='/ParkingDoubleMuonLowMass0/Run2025E-PromptReco-v1/MINIAOD'
...
config.General.requestName = myname
...
config.Data.inputDataset = mydata
config.Data.lumiMask = 'Cert_Collisions2025_391658_398903_Muon.json'

config.Data.outLFNDirBase = '/store/user/<user_name>/CMSDAS/'
```

> #### **Task**
> Edit the `outLFNDirBase` in `crab_upsilon.py`.

Then submit:
```bash
crab submit crab_upsilon.py
crab status -d CernJobs/crab_CMSDAS_Upsilon_2025E_ParkingDoubleMuonLowMass0_v1
```

- CRAB is the scalable path for large data processing.
- `Cert_Collisions2025_391658_398903_Muon.json` is the certified good-luminosity list for this run range.
    - It removes luminosity sections with known detector/DAQ/data-quality problems (for example unstable detector conditions, data-taking interruptions, or subdetector quality flags failing certification).
    - This ensures that accepted events come from periods where detector performance is validated for physics analysis, so your yield and cross-section normalization remain consistent and reproducible.

## Check Data Distributions
Run the plotting program:
```bash
cd /path/to/CMSDAS_Upsilon/data
root -l plot.C
```

Inspect the output varible distributions in `/path/to/CMSDAS_Upsilon/data/results`.

> #### **Task**
> The current program draws variables distributions without applying any selections. To obtain variable distributions of events that are actually used in the analysis, the following selections mentioned earlier should be added inside the event loop:
> - acceptance kinematics for muons: $|\eta|<2.0$ and $p_\mathrm{T}>3.1$ GeV,
> - trigger (`HLT_Dimuon0_Upsilon`) passed,
> - $vProb>0.01$.
>
> **Hints**:
> 1. Event loop location:
> ```cpp
> for (Long64_t i = 0; i < n_entries; ++i) {
>   chain.GetEntry(i);
> ...
> }
> ```
> 2. Besides adding selections, remember to also edit the save path of figures.

Inspect the output varible distributions with selections in the save path you set.

> #### **Question**
> 1. Compare distributions before and after your selections: what observables have changed, and why?
> 2. Is the behavior consistent with trigger and vertex-quality expectations?
