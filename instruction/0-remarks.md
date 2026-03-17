# Part 0 - Remarks

## Facilitators

* [Yiyang Zhao (Tsinghua University)](mailto:yiyang.zhao@cern.ch)
* [Yuxiao Wang (Tsinghua University)](mailto:yuxiao.wang@cern.ch)
* [Junkai Qin  (Tsinghua University)](mailto:junkai.qin@cern.ch)

## Miscellaneous Notes
The core tutorial should be completed in order. `Part 8` and `Part 9` are optional extensions.

All sections use the same three callout styles:
> #### **Task**
> An operation you should _**complete by yourself**_ (for example: edit code, fill a command placeholder, verify an output).

> #### **Question**
> A physics or methodology question for you to _**think and discuss**_.

> #### **Checkpoint**
> A _**required status check**_ before moving forward.

Code excerpts show **only key lines**.
If intermediate lines are omitted, they are marked with `...`.

## Course Structure
- `Part 0`: conventions and environment initialization.
- `Part 1`: physics context, cross section formula, trigger strategy, and analysis scope.
- `Part 2`: data production and first distribution-level validation.
- `Part 3`: signal-yield extraction from dimuon-mass fits in each $(p_\mathrm{T},|y|)$ bin.
- `Part 4`: acceptance $A$ calculation with generated-level MC.
- `Part 5`: efficiency $\epsilon$ calculation from reconstruction-level MC.
- `Part 6`: systematic uncertainty estimation.
- `Part 7`: cross-section calculation, plotting, and validation.
- `Part 8` (optional): event-display with `cmsShow`.
- `Part 9` (optional): polarization effects on acceptance.

## Environment Setup
Each time you start working in your project area on `lxplus8`, initialize the environment first:

```bash
ssh <cern_username>@lxplus8.cern.ch

cd /path/to/CMSDAS_Upsilon

voms-proxy-init -voms cms
voms-proxy-info --timeleft
```

> #### **Checkpoint**
> Continue only after `voms-proxy-info --timeleft` returns a non-zero remaining lifetime.
