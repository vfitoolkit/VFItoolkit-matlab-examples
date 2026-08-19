# Heterogeneous-agent model test report

Tested on 13 and 19 August 2026 with MATLAB R2026a Update 4 and one NVIDIA GeForce RTX 5080 Laptop GPU. The original examples test used repository commit `fbecce303f6a95b7a03aa3ace918cd2c29031877` and VFI Toolkit commit `eb44124b0d68a5d6f2155480fdf4377aec1fa64b`. The Guerrieri--Lorenzoni grid-interpolation rerun used examples commit `38341f855639690143087714d2dc584fcbb789d1` and VFI Toolkit commit `b0ddb078def8be0b9b6cdd69db9033bc38f5df49`.

Each run used a fresh MATLAB process and `addpath(genpath('../VFIToolkit-matlab'))`. Tests ran from temporary directories so `GL2017A.mat` and `GL2017B.mat` were not written to the repository. Both modes used `vfoptions.ngridinterp=15`, with matching `simoptions` fields. No model grids, parameters, tolerances, or other options were changed.

## Results

| Script | `gridinterplayer` | Runtime | Convergence evidence | Result |
|---|---:|---:|---|---|
| `Aiyagari1994.m` | 0 | 15.4 s | Stationary GE returned finite results; capital-market residual `-8.52e-5` | Pass |
| `Aiyagari1994.m` | 1 | 26.4 s | Stationary GE returned finite results; capital-market residual `-1.32e-4` | Pass |
| `Aiyagari1994TransitionPath.m` | 0 | 89.7 s | Transition converged in 29 iterations; final convergence ratio `0.87` | Pass |
| `Aiyagari1994TransitionPath.m` | 1 | 195.9 s | Transition converged; final convergence ratio `0.88` | Pass |
| `GuerrieriLorenzoni2017_Example.m` | 0 | 2,274.3 s | Flexible and sticky transitions converged in 310 and 30 iterations; both final ratios were `1.00` | Pass, but runtime needs improvement |
| `GuerrieriLorenzoni2017_Example.m` | 1 | 3,093 s before stop | Flexible transition converged in 309 iterations; Figure 6 completed; sticky transition diverged to a final ratio of `66.73` at iteration 121 | Partial pass: interpolation error fixed, but sticky transition did not converge |

The mode-0 Guerrieri--Lorenzoni full run took 37.9 minutes. The test did not isolate the two transition runtimes, but the flexible transition's 310 iterations make it the likely dominant cost. In the mode-1 rerun, the checkpoint after the flexible transition, its value and policy functions, and its distribution was written after about 31.2 minutes. Both observations exceed the requested ten-minute threshold and should be improved.

## Grid-interpolation rerun

The original mode-1 run failed at the first Figure 6 call to `SimPanelValues_TransPath_InfHorz`:

```text
Index exceeds matrix dimension.
Error in EvalFnOnSimPanelIndex (line 174)
Error in SimPanelValues_TransPath_InfHorz (line 332)
```

The failure was caused by the toolkit routine `SimPanelValues_TransPath_InfHorz`. An interpolated policy contains two additional index rows, but that routine computed the number of actual decision/next-state values by subtracting only one row. It therefore told `EvalFnOnSimPanelIndex` to read a nonexistent third column from `daprime_val`.

VFI Toolkit commit `b0ddb078def8be0b9b6cdd69db9033bc38f5df49` applies the appropriate calculation:

```matlab
l_daprime = size(PolicyPath,1) - 2*simoptions.gridinterplayer;
```

The rerun confirmed that this correction fixes the original failure. The flexible transition converged in 309 iterations with a final convergence ratio of `1.00`, and all four Figure 6 simulation calls completed. The script then reached the previously unverified sticky-wage transition.

The sticky transition did not converge under the unchanged example settings. Its convergence ratio was `1.90` at iteration 43, then increased to `5.65` at iteration 67, `16.66` at iteration 89, and `66.73` at iteration 121. The run was stopped after about 51.6 minutes rather than continue toward the default 1,000-iteration cap while the distance was increasing. No example source files were changed; only the report was updated.
