# Part 2 - Data

## Datasets
We will use 2025 `ParkingDoubleMuonLowMass` samples in MINIAOD format:

| Era | Dataset | Integrated luminosity recorded by CMS [$\mathrm{fb}^{-1}$] | Run range |
| --- | --- | --- | --- |
| 2025E | `/ParkingDoubleMuonLowMass*/Run2025E-PromptReco-v1/MINIAOD` | 14.2 | 395968-396597 |
| 2025F | `/ParkingDoubleMuonLowMass*/Run2025F-PromptReco-v1/MINIAOD` | 27.3 | 396598-397853 |
| 2025G | `/ParkingDoubleMuonLowMass*/Run2025G-PromptReco-v1/MINIAOD` | 23.0 | 397854-398903 |

In this section, we will use programs caller `ntuplizer` to process these raw datasets in order to:
- Select events containing two muons;
- Convert the file format to `ntuple`, which is more convenient for direct use.

## CMSSW
CMSSW is the software framework used by CMS for simulation, reconstruction, and data analysis. It provides the standard environment and tools to read CMS data formats, run analysis code, and process events within the CMS software ecosystem. For CMS analyses, it is the basic platform for most official workflows.
- [CMS Offline WorkBook (official CMS software introduction and workflow reference)](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBook)
- [CMSSW setup and first steps ](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter1)

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

After production, let's explore the file. The file is stored in `.root` format.
- [ROOT documentation](https://root.cern.ch/root/htmldoc/guides/users-guide/ROOTUsersGuide.html)
- [ROOT tutorials and examples](https://root.cern/doc/master/group__Tutorials.html)


```bash
root -l rootuple.root
...
root [1] rootuple->cd()
root [2] mm_tree->GetEntries()
root [3] mm_tree->Print()
root [4] mm_tree->Show(0)
root [5] mm_tree->Scan("dimuon_p4.Pt():dimuon_p4.Rapidity():trigger","","",5)
root [6] .q
```

The output of `mm_tree->Print()` lists all the available variables in our sample file. The meanings of variables that we will use are as follows:

| Variable | Format | Meaning  |
| --- | --- | --- |
| `run` | `unsigned int` | run number |
| `event` | `unsigned long` | event number |
| `lumiblock` | `unsigned int` | luminosity section number |
| `nonia` | `unsigned int` | index of the onia-candidate within the event |
| `trigger` | `unsigned int` | trigger flag for the selected $\Upsilon$ candidate, `1` means the event passed the `HLT_Dimuon0_Upsilon` trigger |
| `charge` | `int` | sum of the two muon charges |
| `dimuon_p4` | `TLorentzVector` | reconstructed 4-momentum of the $\Upsilon$ candidate, can be used to compute $p_\mathrm{T}$, $y$, mass, $\eta$, $\phi$, etc. |
| `muonP_p4` | `TLorentzVector` | reconstructed 4-momentum of the $\mu^{+}$ |
| `muonM_p4` | `TLorentzVector` | reconstructed 4-momentum of the $\mu^{-}$ |
| `vProb` | `float` | vertex-fit probability of the dimuon candidate, larger values indicate a better fit quality, thus more likely represent real $\Upsilon$ rather than background |


Code-output relation:
- `run_upsilon.py` builds and filters dimuon candidates.
- `MMrootupler` stores selected event content into `rootuple.root/mm_tree`.
- The scanned branches (`dimuon_p4`, `trigger`) are analysis-driving observables used again in later parts.

> #### **Checkpoint**
> Before CRAB submission, confirm:
> - `mm_tree->GetEntries()` is non-zero.
> - `Scan` can read key branches such as `dimuon_p4`, `trigger`, and `vProb`.

### CRAB
CRAB, CMS Remote Analysis Builder, is the standard CMS tool for submitting analysis jobs to distributed computing resources. It is commonly used to process large datasets or Monte Carlo samples on the Grid without running locally. CRAB handles job splitting, submission, monitoring, output collection, and interaction with CMS computing infrastructure. 

We will use the `CRAB` system to process these datasets. In `crab_upsilon.py`, the following lines control request name (`myname`), dataset input (`mydata`), certified-lumi filtering (`lumiMask`), and output directory (`outLFNDirBase`):

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
```

- CRAB is the scalable path for large data processing.
- `Cert_Collisions2025_391658_398903_Muon.json` is the certified good-luminosity list for this run range.
    - It removes luminosity sections with known detector/DAQ/data-quality problems (for example unstable detector conditions, data-taking interruptions, or subdetector quality flags failing certification).
    - This ensures that accepted events come from periods where detector performance is validated for physics analysis, so your yield and cross-section normalization remain consistent and reproducible.

You can use `crab status` to check the status of the jobs. 
```bash
crab status -d CernJobs/crab_CMSDAS_Upsilon_2025E_ParkingDoubleMuonLowMass0_v1
```

**Note**: there is no need to wait for the jobs to finish running completely. Due to time constraints, the subsequent steps will only use the data samples that have already been prepared.

## Check Data Distributions
**Note: DO NOT MODIFY the `kInputFile` path in `plot.C`. The existing path can be accessed.**
Run the plotting program:
```bash
cd /path/to/CMSDAS_Upsilon/data
root -l plot.C
```

Inspect the output variable distributions in `/path/to/CMSDAS_Upsilon/data/results`.

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
>   if (dimuon_p4 == nullptr || muonP_p4 == nullptr || muonM_p4 == nullptr) {
>   continue;
>   }
>   
>   /*EDIT HERE*/
>   ...
> }
> ```
> 2. Besides adding selections, remember to also edit the save path of figures.

Inspect the output variable distributions with selections in the save path you set.

> #### **Question**
> 1. Compare distributions before and after your selections: what observable distributions have changed, and why? Examine `dimuon_mass.pdf`. After applying these selections, has the background level under the $\Upsilon$ signal peak decreased? Based on this observation, explain the purpose of applying these filters.
> #### **Task**
> Try applying a stricter selection criterion, such as $vProb>0.1$. Observe how the observable distributions change. 
> #### **Question**
> 2. Examine `dimuon_mass.pdf`. Has the background level beneath the $\Upsilon$ signal peak decreased further? 
> 3. Consider why we did not use this stricter selection criterion to further reduce the background.
