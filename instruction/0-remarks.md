# Part 0 - Remarks

## Miscellaneous Notes
The core tutorial should be completed in order. `Part 8` and `Part 9` are optional extensions.

All sections use the same three callout styles:
> #### **Task**
> An operation you should complete yourself (for example: edit code, fill a command placeholder, verify an output).

> #### **Question**
> A physics or methodology question for your discussion.

> #### **Checkpoint**
> A required status check before moving forward.

Code excerpts show only key lines.
If intermediate lines are omitted, they are marked with `...`.

## Course Structure
- `Part 0`: conventions and environment initialization.
- `Part 1`: physics context, measurement formula, trigger strategy, and analysis scope.
- `Part 2`: data production workflow and first distribution-level validation.
- `Part 3`: signal-yield extraction from dimuon-mass fits in $(p_\mathrm{T},|y|)$ bins.
- `Part 4`: acceptance $A$ from generated-level kinematic phase space.
- `Part 5`: efficiency $\epsilon$ from reconstruction-level event selection in MC.
- `Part 6`: systematic-uncertainty strategy and acceptance-related example.
- `Part 7`: luminosity, cross-section computation, plotting, and validation.
- `Part 8` (optional): event-display workflow with `cmsShow`.
- `Part 9` (optional): polarization impact on acceptance.

## Environment Setup
Each time you start working in your project area on `lxplus8`, initialize the environment first:

```bash
ssh <cern_username>@lxplus8.cern.ch

cd /path/to/CMSDAS_Upsilon

source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init -voms cms
voms-proxy-info --timeleft
```

> #### **Checkpoint**
> Continue only after `voms-proxy-info --timeleft` returns a non-zero remaining lifetime.
