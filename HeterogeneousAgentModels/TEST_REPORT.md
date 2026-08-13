# Heterogeneous-agent model test report

Tested on 13 August 2026 with MATLAB R2026a Update 4 and one available GPU. The examples repository was at `fbecce303f6a95b7a03aa3ace918cd2c29031877`; VFI Toolkit was at `eb44124b0d68a5d6f2155480fdf4377aec1fa64b`.

Each run used a fresh MATLAB process and `addpath(genpath('../VFIToolkit-matlab'))`. Tests ran from temporary directories so `GL2017A.mat` and `GL2017B.mat` were not written to the repository. Both modes used `vfoptions.ngridinterp=15`, with matching `simoptions` fields. No model grids, parameters, tolerances, or other options were changed.

## Results

| Script | `gridinterplayer` | Runtime | Convergence evidence | Result |
|---|---:|---:|---|---|
| `Aiyagari1994.m` | 0 | 15.4 s | Stationary GE returned finite results; capital-market residual `-8.52e-5` | Pass |
| `Aiyagari1994.m` | 1 | 26.4 s | Stationary GE returned finite results; capital-market residual `-1.32e-4` | Pass |
| `Aiyagari1994TransitionPath.m` | 0 | 89.7 s | Transition converged in 29 iterations; final convergence ratio `0.87` | Pass |
| `Aiyagari1994TransitionPath.m` | 1 | 195.9 s | Transition converged; final convergence ratio `0.88` | Pass |
| `GuerrieriLorenzoni2017_Example.m` | 0 | 2,274.3 s | Flexible and sticky transitions converged in 310 and 30 iterations; both final ratios were `1.00` | Pass, but runtime needs improvement |
| `GuerrieriLorenzoni2017_Example.m` | 1 | 3,210.1 s before error | Flexible transition returned after 2,587.8 s; Figure 6 then failed before the sticky transition | Fail; runtime also needs improvement |

The mode-0 Guerrieri--Lorenzoni full run took 37.9 minutes. The test did not isolate the two transition runtimes, but the flexible transition's 310 iterations make it the likely dominant cost. In mode 1, the flexible transition alone took 43.1 minutes. Both observations exceed the requested ten-minute threshold and should be improved.

## Interpolation failure

With `gridinterplayer=1`, the script fails at the first Figure 6 call to `SimPanelValues_TransPath_InfHorz`:

```text
Index exceeds matrix dimension.
Error in EvalFnOnSimPanelIndex (line 174)
Error in SimPanelValues_TransPath_InfHorz (line 332)
```

The failure is caused by the read-only toolkit routine `SimPanelValues_TransPath_InfHorz`. An interpolated policy contains two additional index rows, but that routine computes the number of actual decision/next-state values by subtracting only one row. It therefore tells `EvalFnOnSimPanelIndex` to read a nonexistent third column from `daprime_val`.

Comparable toolkit routines use the appropriate calculation:

```matlab
l_daprime = size(PolicyPath,1) - 2*simoptions.gridinterplayer;
```

The suggested fix is to make the same change near lines 191--195 of the toolkit's `SimPanelValues_TransPath_InfHorz.m`. This was not applied because the toolkit repository was explicitly out of scope. An example-side workaround would require copying or shadowing a substantial toolkit routine, which would be fragile and is not recommended.

The saved mode-1 state was used to begin a follow-up test of the sticky transition while skipping only the known failing Figure 6 block. That run was stopped at the user's request because it was taking too long, so the mode-1 sticky transition remains unverified. No Aiyagari source correction was needed, and no example source files were changed.
