# Part 8 - Event Display (Optional)

This part shows how to visualize event topology with cmsShow.

## Python Setup
**Note**: open a new terminal or reconnect to the `lxplus8` server before continuing.

```bash
cd /path/to/CMSDAS_Upsilon/data/CMSSW_15_0_18/src
cmsenv

cd /path/to/CMSDAS_Upsilon/event_display

python3 -c "import sys; print(sys.version)"
python3 -c "import uproot,numpy; print('event display python setup ok')"
```

> #### **Checkpoint**
> Continue only after both `python3 -c` tests complete successfully.

## Pick Events

Run the script (**do not** edit the input ntuple and dataset):

```bash
python3 event_display.py \
  --dataset "/ParkingDoubleMuonLowMass0/Run2025G-PromptReco-v1/MINIAOD" \
  --ntuple "/eos/home-y/yiyangz/public/CMSDAS/data/2025G_Parking0.root" \
  --max-events 10 \
  --output-dir results
```

Outputs:
- `/path/to/CMSDAS_Upsilon/event_display/results/event_display.txt`: event IDs in the format required by `edmPickEvents.py`.
- `/path/to/CMSDAS_Upsilon/event_display/results/selected_events.csv`: tabulated metadata (`run`, `lumi`, `event`).
- `/path/to/CMSDAS_Upsilon/event_display/results/pickevents.root`: picked event content used for visual inspection.

> #### **Task**
> Run the script on your ntuple and verify that 10 events were selected and written.
> #### **Checkpoint**
> Confirm that `pickevents.root` exists and is non-empty before moving to cmsShow.

## CERNBox Preparation
Before opening cmsShow, upload `pickevents.root` to CERNBox (`/eos/home-<u>/<user_name>/*`) and grant `viewer` permission to `cms-vis-access`.

Steps:
1. Open [CERNBox](https://cernbox.cern.ch).
2. Create a folder. Upload `pickevents.root` to this folder.
3. As shown below, share the folder with `cms-vis-access` using `invite` with `viewer` permission.

<p align="center">
  <img src="figures/access.png" alt="CERNBox sharing setup" width="30%">
</p>

## Open Events in cmsShow
Open the cmsShow web interface:
- [https://fireworks.cern.ch](https://fireworks.cern.ch)

Then:
1. Load `pickevents.root` from your CERNBox location using `Load File EOS`.
2. Open one event and inspect global topology.
3. Add/remove collections as needed from the left panel.

![cmsShow event view](figures/event_display_web.png)

> #### **Task**
> Inspect all 10 selected events and answer the following questions.
> 
> Pick the one you think looks best to use as the background picture for your final presentation slides!
> #### **Question**
> 1. Which reconstructed features are easiest to interpret visually in cmsShow for your selected events?
> 2. Which details are difficult to evaluate from event display alone and still require quantitative analysis plots?
